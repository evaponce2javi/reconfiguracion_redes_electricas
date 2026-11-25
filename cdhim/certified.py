"""Certified brackets of the exact range of the voltage magnitudes over an uncertainty set.

Inside: a point of the set whose power-flow solution is verified (cdhim.verified.verified_point) proves that the
exact maximum of |V_i| is at least the lower end of that point's enclosure, and the exact minimum at most its
upper end. Outside: the union of the verified enclosures of a partition of the box contains the exact range, and
it shrinks to it as the partition is refined (the excess of every enclosure is of third order in the box size).
Together they bracket the exact extremes:

    vmin in [vmin_lo, vmin_hi],     vmax in [vmax_lo, vmax_hi].

All four are floating-point rigorous for the equations and the set defined by the stored data (Ybus, the centre
injections and the injection matrix); the rounding of every point and sub-box centre computed from them is
accounted for (``s_err`` of the verifier).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import scipy.sparse as sp

from types import SimpleNamespace

from .discs import disc_bound
from .krawczyk import Krawczyk
from .network import Grid
from .powerflow import newton_pf
from .second_order import SecondOrder
from .study import Study, custom_study
from .verified import abs_up, add_up, err_const, matmul_up, up, verified_enclosure, verified_point


def matrix_study(grid: Grid, sbus: np.ndarray, M: np.ndarray, **kw) -> Study:
    """Study with the injections sbus + M eps, eps in [-1, 1]^K, for a complex matrix M (n_bus, K)."""
    M = np.asarray(M, dtype=complex)
    bp = np.flatnonzero(np.abs(M.real).sum(axis=1) > 0)
    bq = np.flatnonzero(np.abs(M.imag).sum(axis=1) > 0)
    Z = sp.csr_matrix(np.vstack([M.real[bp], M.imag[bq]]))
    return custom_study(grid, sbus, np.r_[bp, bq], np.r_[np.zeros(len(bp), bool), np.ones(len(bq), bool)],
                        np.ones(len(bp) + len(bq)), Z=Z, **kw)


def injection_matrix(st: Study, delta: float = 1.0) -> np.ndarray:
    """Complex change of the specified injection of every bus per unit of every noise symbol, (n_bus, K): the
    same numbers as cdhim.enclosure.injection_changes, from the study alone."""
    src, n = st.src, st.grid.n_bus
    M = sp.csr_matrix((np.where(src.is_q, 1j, 1.0) * src.unit, (src.bus, np.arange(src.n_injections))), shape=(n, src.n_injections))
    if src.Z is not None:
        M = M @ src.Z
    if src.G is not None:
        M = M + src.G
    return np.asarray((M * delta).todense())


def point_error(grid: Grid, sbus: np.ndarray, M: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Bound of |fl(sbus + M eps) - (sbus + M eps)| at every PQ bus."""
    Mn = abs_up(np.asarray(M)[grid.pq])
    return up(err_const(M.shape[1] + 1) * add_up(abs_up(np.asarray(sbus, dtype=complex)[grid.pq]), matmul_up(Mn, np.abs(eps), M.shape[1])))


def certify_points(grid: Grid, sbus: np.ndarray, M: np.ndarray, eps: np.ndarray, tol: float = 1e-9):
    """Verified (lo, hi) of |V_i| at every PQ bus for each point eps[:, k] of the box: (n_pq, n_points) each.
    NaN columns where the power flow did not converge or the verification failed."""
    n = eps.shape[1]
    lo, hi = np.full((len(grid.pq), n), np.nan), np.full((len(grid.pq), n), np.nan)
    V = grid.V0
    cache = {}                                           # one approximate inverse for all the points
    for k in range(n):
        s = sbus + M @ eps[:, k]
        Vk, ok, _, _ = newton_pf(grid.Ybus, s, V, grid.pv, grid.pq, tol=tol, max_iter=40)
        if not ok:
            Vk, ok, _, _ = newton_pf(grid.Ybus, s, grid.V0, grid.pv, grid.pq, tol=tol, max_iter=40)
        if ok:
            lo[:, k], hi[:, k], _ = verified_point(grid, s, Vk, point_error(grid, sbus, M, eps[:, k]), cache)
    return lo, hi


def in_half_disks(eps: np.ndarray) -> bool:
    """Exact test (rational arithmetic) that eps = (eps_p, eps_q) lies in the product of half disks whose bounding
    boxes are the box: u = (1 + eps_p)/2 >= 0 and u^2 + eps_q^2 <= 1 for every source."""
    n = len(eps) // 2
    for a, b in zip(eps[:n], eps[n:]):
        u, w = (1 + Fraction(float(a))) / 2, Fraction(float(b))
        if u < 0 or u * u + w * w > 1:
            return False
    return True


def half_disk_eps(points: np.ndarray) -> np.ndarray:
    """Points u + j w of the unit half disk, (n_points, n_sources), as noise symbols (2 n_sources, n_points),
    pulled inwards until the exact membership test passes."""
    out = np.empty((2 * points.shape[1], points.shape[0]))
    for k, pt in enumerate(points):
        f = 1.0
        while True:
            e = np.r_[2.0 * np.maximum(pt.real * f, 0.0) - 1.0, pt.imag * f]
            if in_half_disks(e):
                break
            f *= 1.0 - 1e-12
        out[:, k] = e
    return out


@dataclass
class Bracket:
    """vmin in [vmin_lo, vmin_hi] and vmax in [vmax_lo, vmax_hi] for every PQ bus."""

    vmin_lo: np.ndarray
    vmin_hi: np.ndarray
    vmax_lo: np.ndarray
    vmax_hi: np.ndarray
    leaves: int                 # boxes of the final partition
    evaluated: int              # verified enclosures computed
    unverified: int             # leaves without a verified enclosure (then the outer ends are not valid: NaN)
    gap: float                  # largest excess of a leaf over the attained values, relative to the attained width
    root_lo: np.ndarray         # the enclosure of the whole box, for reference
    root_hi: np.ndarray


def box_bracket(grid: Grid, sbus: np.ndarray, M: np.ndarray, tol: float = 1e-3, max_boxes: int = 2000, node_budget: int = 0,
                pf_tol: float = 1e-9, min_width: float = 1e-9, log=None, leaf: str = "verified") -> Bracket:
    """Bracket of the exact extremes of |V_i| over {sbus + M eps, eps in [-1, 1]^K} by adaptive bisection.

    Every leaf (centre c, half-widths h, both dyadic) gets the verified enclosure of the study centred at
    sbus + M c with matrix M diag(h). The leaf with the largest excess over the best attained values is bisected
    along the symbol with the largest first-order effect, until the excess is below ``tol`` times the attained
    width or ``max_boxes`` enclosures have been computed. Attained values: the centres of the boxes and, at the
    root, the vertices that the first-order model points to; the best point per bus is verified at the end.

    ``leaf="discs"`` uses, in every leaf, the per-bus disc of Wang et al. (cdhim.discs) around the solution at
    the leaf centre instead of the verified enclosure: the subdivision of the uncertainty set is how that
    method is made tight in its own paper. Its outer ends are valid in exact arithmetic only (no rounding
    control), and the leaf is bisected along the symbol with the largest injection.
    ``leaf="krawczyk"`` uses the Krawczyk operator in current-injection form with slopes and rotated axes
    (cdhim.krawczyk), also without rounding control.
    """
    M = np.asarray(M, dtype=complex)
    sbus = np.asarray(sbus, dtype=complex)
    K, pq = M.shape[1], grid.pq
    n_pq = len(pq)

    def evaluate(c, h):
        try:
            st = matrix_study(grid, sbus + M @ c, M * h, tol=pf_tol)
        except RuntimeError:
            return None
        so = SecondOrder(st, 1.0)
        enc = verified_enclosure(so, node_budget=node_budget, s_err=point_error(grid, sbus, M, c))
        return dict(c=c, h=h, lo=enc.lo, hi=enc.hi, ok=enc.verified, vm=np.abs(st.V_star[pq]), infl=np.abs(so.c_lin).sum(axis=0), st=st)

    def evaluate_discs(c, h):
        s_c = sbus + M @ c
        Vc, ok, _, _ = newton_pf(grid.Ybus, s_c, grid.V0, grid.pv, grid.pq, tol=pf_tol, max_iter=40)
        if not ok:
            return None
        db = disc_bound(grid, s_c, Vc, np.abs(M) @ h)
        vm = np.abs(Vc[pq])
        rad = db.radius_bus if db.ok else np.full(n_pq, np.nan)
        return dict(c=c, h=h, lo=vm - rad, hi=vm + rad, ok=bool(db.ok), vm=vm, infl=np.abs(M[pq]).sum(axis=0) * h, st=None)

    def evaluate_krawczyk(c, h):
        try:
            st = matrix_study(grid, sbus + M @ c, M * h, tol=pf_tol)
        except RuntimeError:
            return None
        kr = Krawczyk(SimpleNamespace(st=st, delta=1.0), "current", True, True)
        res = kr.solve()
        return dict(c=c, h=h, lo=res.vm_lo, hi=res.vm_hi, ok=bool(res.closes), vm=np.abs(st.V_star[pq]), infl=np.abs(kr.L).sum(axis=0), st=st)

    if leaf not in ("verified", "discs", "krawczyk"):
        raise ValueError(leaf)
    if leaf == "discs":
        evaluate = evaluate_discs
    elif leaf == "krawczyk":
        evaluate = evaluate_krawczyk
    root = evaluate(np.zeros(K), np.ones(K))
    if root is None:
        raise RuntimeError("no AC solution at the centre of the box")
    st0 = root.pop("st")
    if st0 is None:
        st0 = matrix_study(grid, sbus, M, tol=pf_tol)
    best_max, best_min = root["vm"].copy(), root["vm"].copy()
    arg_max, arg_min = np.zeros((n_pq, K)), np.zeros((n_pq, K))
    # vertices suggested by the first-order model of the root
    c1 = SecondOrder(st0, 1.0).c_lin
    sg = np.where(c1 >= 0, 1.0, -1.0)
    cand = np.unique(np.vstack([sg, -sg]), axis=0)
    V, conv, _ = st0.solver.solve(st0.injections(cand.T, 1.0))
    vm = np.abs(V[pq])
    for k in np.flatnonzero(conv):
        up_, dn_ = vm[:, k] > best_max, vm[:, k] < best_min
        best_max[up_], arg_max[up_] = vm[up_, k], cand[k]
        best_min[dn_], arg_min[dn_] = vm[dn_, k], cand[k]

    leaves, evaluated = [root], 1
    while True:
        width = best_max - best_min
        act = width > min_width
        ex = np.array([np.inf if not lf["ok"] else
                       (np.maximum(lf["hi"] - best_max, best_min - lf["lo"])[act] / width[act]).max() if act.any() else 0.0 for lf in leaves])
        b = int(np.argmax(ex))
        if ex[b] <= tol or evaluated + 2 > max_boxes:
            break
        lf = leaves.pop(b)
        d = int(np.argmax(lf["infl"]))                   # first-order width of every symbol inside this leaf
        for sgn in (-1.0, 1.0):
            c, h = lf["c"].copy(), lf["h"].copy()
            h[d] *= 0.5
            c[d] += sgn * h[d]
            child = evaluate(c, h)
            evaluated += 1
            if child is None:
                child = dict(c=c, h=h, lo=np.full(n_pq, np.nan), hi=np.full(n_pq, np.nan), ok=False, vm=None, infl=lf["infl"])
            else:
                child.pop("st")
                up_, dn_ = child["vm"] > best_max, child["vm"] < best_min
                best_max[up_], arg_max[up_] = child["vm"][up_], c
                best_min[dn_], arg_min[dn_] = child["vm"][dn_], c
            leaves.append(child)
        if log is not None and evaluated % 200 < 2:
            log(f"{evaluated} enclosures, {len(leaves)} leaves, largest excess {ex[b]:.2e}")

    unverified = sum(not lf["ok"] for lf in leaves)
    if unverified:
        out_lo, out_hi = np.full(n_pq, np.nan), np.full(n_pq, np.nan)
    else:
        out_lo, out_hi = np.min([lf["lo"] for lf in leaves], axis=0), np.max([lf["hi"] for lf in leaves], axis=0)
    # verify the best attained point of every bus
    pts, inv = np.unique(np.vstack([arg_max, arg_min]), axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    lo, hi = certify_points(grid, sbus, M, pts.T, pf_tol)
    j = np.arange(n_pq)
    return Bracket(out_lo, hi[j, inv[n_pq:]], lo[j, inv[:n_pq]], out_hi, len(leaves), evaluated, unverified, float(ex[b]), root["lo"], root["hi"])


def half_disk_inner(st: Study, points_max: np.ndarray, points_min: np.ndarray, tol: float = 1e-9):
    """Certified inner bounds over the product of half disks of the study ``st`` (cdhim.halfdisk.half_disk_study).

    ``points_max[i]``, ``points_min[i]``: the point (u + j w per source) at which the largest and the smallest
    |V_i| were found. Returns (vmin_hi, vmax_lo): the exact minimum of bus i is at most vmin_hi[i] and the
    exact maximum at least vmax_lo[i]; NaN where the point could not be verified.
    """
    g = st.grid
    n_pq = len(g.pq)
    M = injection_matrix(st)
    pts, inv = np.unique(np.vstack([points_max, points_min]), axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    lo, hi = certify_points(g, st.sbus, M, half_disk_eps(pts), tol)
    j = np.arange(n_pq)
    return hi[j, inv[n_pq:]], lo[j, inv[:n_pq]]
