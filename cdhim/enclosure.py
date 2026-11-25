"""Fixed-point enclosure of the power-flow solution set under box uncertainty (prototype, not verified).

In rectangular coordinates the state deviation D = x - x* solves, exactly,

    D = C eps - J^-1 q(D),        q(a + r) = q(a) + dJ(a) r + q(r),

where q is the quadratic part of the equations and dJ(a) = J(x* + a) - J(x*) is LINEAR in a.
Write D = C eps + r. The remainder r is a fixed point of

    T_eps(r) = -J^-1 [ q(C eps) + dJ(C eps) r + q(r) ] = r2(eps) - sum_k eps_k P_k r - J^-1 q(r),

with r2(eps) = -J^-1 q(C eps) (a quadratic form per state component) and P_k = J^-1 dJ(C_k).
For |r| <= rho componentwise and every eps in [-1, 1]^K:

    |T_eps(r)|  <=  b2 + G rho + N(rho),
    b2  >= |r2(eps)|                   (bound of the state-level quadratic forms),
    G   = sum_k |P_k|                  (entrywise; the cancellations inside J^-1 dJ are kept),
    N(rho) = 1/2 sum_j rho_j |J^-1 dJ(e_j)| rho     (q(r) = 1/2 dJ(r) r, so the same cancellations are kept);
             for large systems the cruder |J^-1| qbar(rho), qbar_m = rho_m (|Y| rho)_m, is used instead.

If rho* >= b2 + G rho* + N(rho*) the box [-rho*, rho*] is mapped into itself for every eps, and by
Brouwer's theorem the power-flow equations have a solution with |D - C eps| <= rho* for every eps of
the box. Then s = r - r2(eps) satisfies |s| <= sigma = G rho* + N(rho*), which is of third order.

For the magnitude of a PQ bus, with w = (C eps + r2)_i, p and q its components along and across V*_i:

    |V*_i + w| - |V*_i| = p + q^2 / (2 (v + p)) - t q^4 / (8 (v + p)^3),   0 <= t <= 1,

so d|V_i| = c_i' eps + eps' H_i eps + rem with |rem| <= e3 + |s_i|, where e3 collects the terms of the
expansion of the modulus above the second order, bounded with intervals (see ``modulus_allowance``).

Everything is valid in exact arithmetic. Nothing is rounded outward, and uniqueness is not shown.

``InjectionEnclosure`` is a second bound of the same remainder for networks whose non-slack buses are all PQ.
It writes the equations in current-injection form, where the non-linearity is local to every bus, so it does
not degrade with very short lines (large entries of Y) as the bound above does.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from scipy.sparse.linalg import LinearOperator, eigs

from .knomial import knomial_bound
from .powerflow import polar_jacobian, solve_blocked
from .second_order import SecondOrder
from .study import drive_matrix

ASSEMBLED_QUADRATIC_UP_TO = 1500      # state dimension up to which N(rho) keeps the cancellations of J^-1 dJ(e_j)
INFLATION = 1e-9                      # relative inflation of the computed fixed point before the self-mapping check


def rect_jacobian(grid, V):
    """Jacobian of the rectangular power-flow equations at complex voltages ``V``; linear in V.

    Rows: P at pv + pq, Q at pq, |V|^2 at pv. Columns: e, then f, of the non-slack buses."""
    Y = grid.Ybus.tocsr()
    ns, pq, pv = grid.pvpq, grid.pq, grid.pv
    m = len(ns)
    col = np.full(grid.n_bus, -1)
    col[ns] = np.arange(m)
    cur = Y @ V
    dS_de = (sp.diags(V) @ Y.conj() + sp.diags(cur.conj())).tocsr()
    dS_df = (1j * (sp.diags(cur.conj()) - sp.diags(V) @ Y.conj())).tocsr()
    blocks = [sp.hstack([dS_de[ns][:, ns].real, dS_df[ns][:, ns].real]), sp.hstack([dS_de[pq][:, ns].imag, dS_df[pq][:, ns].imag])]
    if len(pv):
        k = np.arange(len(pv))
        blocks.append(sp.hstack([sp.csr_matrix((2 * V[pv].real, (k, col[pv])), shape=(len(k), m)),
                                 sp.csr_matrix((2 * V[pv].imag, (k, col[pv])), shape=(len(k), m))]))
    return sp.vstack(blocks).tocsc()


def state_forms(so: SecondOrder) -> np.ndarray:
    """r2 = eps' Hs eps for every state component: (2m, K, K). Liao et al. (2017), Eq. 25, halved."""
    g = so.st.grid
    E, F = so.E, so.F
    k = E.shape[1]
    z = 0.5 * (E[:, :, None] * np.conj(F)[:, None, :] + E[:, None, :] * np.conj(F)[:, :, None])
    rows = [z[g.pvpq].real.reshape(-1, k * k), z[g.pq].imag.reshape(-1, k * k)]
    if so.has_pv:
        rows.append((np.conj(E[g.pv])[:, :, None] * E[g.pv][:, None, :]).real.reshape(-1, k * k))
    return -solve_blocked(so.lu, np.vstack(rows)).reshape(2 * so.m, k, k)


def monomial_bounds(H: np.ndarray):
    """Term-by-term enclosure of eps' H eps on [-1, 1]^K for a stack of forms (n, K, K)."""
    d = np.einsum("ikk->ik", H)
    off = np.abs(H).sum(axis=(1, 2)) - np.abs(d).sum(axis=1)
    return -np.maximum(-d, 0.0).sum(axis=1) - off, np.maximum(d, 0.0).sum(axis=1) + off


@dataclass
class FixedPoint:
    spectral_radius: float      # of G; the linear part contracts only below 1
    converged: bool
    iterations: int
    rho: np.ndarray             # |D - C eps| <= rho for every eps of the box (per state component)
    sigma: np.ndarray           # |D - C eps - r2(eps)| <= sigma: third order and above
    b2: np.ndarray              # bound of the second-order state term


class Enclosure:
    """Operators of the fixed-point test at delta = so.delta; ``solve(scale)`` rescales to another half-width."""

    def __init__(self, so: SecondOrder):
        g = so.st.grid
        self.so, self.g = so, g
        n2 = 2 * so.m
        self.Hs = state_forms(so)
        lo, hi = monomial_bounds(self.Hs)
        self.b2 = np.maximum(-lo, hi)
        G = np.zeros((n2, n2))
        for k in range(so.E.shape[1]):                               # G = sum_k |J^-1 dJ(C_k)|
            G += np.abs(solve_blocked(so.lu, rect_jacobian(g, so.E[:, k]).toarray()))
        self.G = G
        self.Jinv_abs = np.abs(solve_blocked(so.lu, np.eye(n2)))
        self.Yabs = abs(g.Ybus).tocsr()
        self.quad = None
        if n2 <= ASSEMBLED_QUADRATIC_UP_TO:                           # |J^-1 dJ(e_j)|, stored by its non-zero columns
            self.quad = []
            for j in range(n2):
                unit = np.zeros(g.n_bus, complex)
                unit[g.pvpq[j % so.m]] = 1.0 if j < so.m else 1j
                A = rect_jacobian(g, unit)
                cols = np.flatnonzero(np.diff(A.indptr))
                self.quad.append((cols, np.abs(so.lu.solve(A[:, cols].toarray()))))
        x = np.ones(n2)
        for _ in range(300):                                         # Perron root of the non-negative G
            y = G @ x
            lam = float(np.linalg.norm(y))
            x = y / lam
        self.spectral_radius = lam

    def qbar(self, rho):
        """Interval bound of every row of the quadratic part for |r| <= rho."""
        g, m = self.g, self.so.m
        mag = np.zeros(g.n_bus)
        mag[g.pvpq] = np.hypot(rho[:m], rho[m:])
        prod = mag * (self.Yabs @ mag)                               # |dV_m| |(Y dV)_m|
        return np.concatenate([prod[g.pvpq], prod[g.pq], mag[g.pv] ** 2])

    def N(self, rho):
        """Bound of |J^-1 q(r)| for |r| <= rho."""
        if self.quad is None:
            return self.Jinv_abs @ self.qbar(rho)
        out = np.zeros(len(rho))
        for j, (cols, P) in enumerate(self.quad):
            if rho[j] > 0.0:
                out += rho[j] * (P @ rho[cols])
        return 0.5 * out

    def true_radius(self, n_vertices: int = 64) -> float:
        """Largest spectral radius of J^-1 dJ(C eps) over a sample of vertices: what G over-estimates."""
        so, g = self.so, self.g
        n2, K = 2 * so.m, so.E.shape[1]
        rng = np.random.default_rng(0)
        worst = 0.0
        for t in range(n_vertices):
            eps = np.ones(K) if t == 0 else np.sign(rng.standard_normal(K))
            A = rect_jacobian(g, so.E @ eps)
            op = LinearOperator((n2, n2), matvec=lambda x, A=A: so.lu.solve(A @ x), dtype=float)
            try:
                lam = eigs(op, k=1, which="LM", tol=1e-6, maxiter=2000, return_eigenvectors=False)
                worst = max(worst, float(np.abs(lam[0])))
            except Exception:                                         # no convergence: skip the vertex
                continue
        return worst

    def solve(self, scale: float = 1.0, max_iter: int = 2000, tol: float = 1e-13) -> FixedPoint:
        """Smallest rho with rho >= b2 + G rho + N(rho) at half-width scale * delta (b2 ~ scale^2, G ~ scale)."""
        b2, G = self.b2 * scale * scale, self.G * scale
        rho = b2.copy()
        ok = False
        for it in range(1, max_iter + 1):
            new = b2 + G @ rho + self.N(rho)
            if not np.all(np.isfinite(new)) or new.max() > 1e3:
                break
            done = np.abs(new - rho).max() <= tol * max(new.max(), 1e-300)
            rho = new
            if done:
                rho = rho * (1.0 + INFLATION)                        # the iterates approach the fixed point from below
                ok = bool(np.all(b2 + G @ rho + self.N(rho) <= rho))  # the box is mapped into itself
                break
        sigma = G @ rho + self.N(rho)
        return FixedPoint(self.spectral_radius * scale, ok, it, rho, sigma, b2)

    def modulus_allowance(self, fp: FixedPoint, scale: float = 1.0):
        """(e3 + |s_i|, e3, |s_i| bound) for every PQ bus; inf where the bound of the modulus breaks down."""
        so, g = self.so, self.g
        e3 = modulus_higher_order(so, self.Hs, scale)
        s_mag = np.hypot(fp.sigma[so.col[g.pq]], fp.sigma[so.m + so.col[g.pq]])
        return e3 + s_mag, e3, s_mag


def modulus_higher_order(so: SecondOrder, Hs: np.ndarray, scale: float = 1.0) -> np.ndarray:
    """e3: bound of the terms above second order in the expansion of |V*_i + w|, w the first- plus second-order
    state deviation, for every PQ bus (see the module docstring); inf where the bound breaks down."""
    g = so.st.grid
    ie, jf = so.col[g.pq], so.m + so.col[g.pq]
    v = so.vm
    er, ei = (so.st.V_star[g.pq].real / v)[:, None, None], (so.st.V_star[g.pq].imag / v)[:, None, None]
    s2 = scale * scale
    p2 = np.maximum(*[np.abs(b) for b in monomial_bounds(er * Hs[ie] + ei * Hs[jf])]) * s2     # along V*
    q2 = np.maximum(*[np.abs(b) for b in monomial_bounds(er * Hs[jf] - ei * Hs[ie])]) * s2     # across V*
    q1 = np.abs(so.c_perp).sum(axis=1) * scale
    pm = np.abs(so.c_lin).sum(axis=1) * scale + p2
    den = v - pm
    with np.errstate(divide="ignore", invalid="ignore"):
        e3 = q1 ** 2 * pm / (2.0 * v * den) + (2.0 * q1 * q2 + q2 ** 2) / (2.0 * den) + (q1 + q2) ** 4 / (8.0 * den ** 3)
    return np.where(den > 0, e3, np.inf)


def injection_changes(so: SecondOrder) -> np.ndarray:
    """Complex change of the specified injection of every bus per unit of every noise symbol: (n_bus, K)."""
    src, n = so.st.src, so.st.grid.n_bus
    M = sp.csr_matrix((np.where(src.is_q, 1j, 1.0) * src.unit, (src.bus, np.arange(src.n_injections))), shape=(n, src.n_injections))
    if src.Z is not None:
        M = M @ src.Z
    if src.G is not None:
        M = M + src.G
    return np.asarray((M * so.delta).todense())


@dataclass
class InjectionFixedPoint:
    converged: bool
    iterations: int
    s_abs: np.ndarray           # per bus (all buses, zero at the slack): |V - V* - first order - second order| <= s_abs
    rel_radius: float           # largest bound of |V_m - V*_m| / |V*_m|; the argument needs it below 1


class InjectionEnclosure:
    """Fixed-point bound of the state remainder above second order, in current-injection form.

    The exact equations are J D = dS - D conj(dI), with dI = Y D the change of the injected currents. For a
    network whose non-slack buses are all PQ the injected current of bus m is conj(S_m / V_m), so
    conj(dI) = (S* + dS)/(V* + D) - S*/V* is a function of the voltage of the same bus only. With
    d = D/V*, sig = S*/V*^2 and tau(eps) = dS(eps)/V*^2 (all per bus),

        J D = dS + V*^2 [ -tau d + (sig + tau) d^2 / (1 + d) ],

    J being the power-flow Jacobian at V* acting on complex power perturbations. The first-order terms give
    D1 = C eps and the second-order ones, V*^2 (-tau d1 + sig d1^2), give D2 = r2(eps) (the forms of
    ``state_forms``). What is left, s = D - D1 - D2, is a fixed point of

        s = J^-1 [ V*^2 c(s) ],
        c = -tau (d - d1) + sig (d^2 - d1^2) + tau d^2 - (sig + tau) d^3 / (1 + d),    d = d1 + d2 + s/V*.

    With |d1| <= a1, |d2| <= a2, |tau| <= t, |s/V*| <= z and a = a1 + a2 + z < 1, per bus,

        |c| <= t (a2 + z) + |sig| (a2 + z)(2 a1 + a2 + z) + t a^2 + (|sig| + t) a^3 / (1 - a) =: beta(z),

    and |s_m| <= sum_n T_mn |V*_n|^2 beta_n, T_mn the spectral norm of the 2x2 real block of J^-1 that maps
    the complex power perturbation of bus n to the complex voltage of bus m. If some z with a < 1 satisfies
    T (|V*|^2 beta(z)) / |V*| <= z, the product of discs |s_m| <= |V*_m| z_m is mapped into itself for every
    eps of the box and Brouwer's theorem gives a solution inside it. ``solve`` iterates from zero, inflates
    the limit slightly and checks that inequality.

    Buses without injection have sig = tau = 0 and do not contribute; Y does not appear in the bound.
    Valid in exact arithmetic; nothing is rounded outward and uniqueness is not shown.
    """

    def __init__(self, so: SecondOrder, Hs: np.ndarray | None = None):
        g, st = so.st.grid, so.st
        if len(g.pv):
            raise ValueError("the current-injection bound needs all non-slack buses to be PQ")
        self.so, self.g = so, g
        m, ns = so.m, g.pvpq
        self.Hs = state_forms(so) if Hs is None else Hs
        V = st.V_star[ns]
        self.v_abs = np.abs(V)
        self.sig = np.abs(st.sbus[ns]) / self.v_abs ** 2
        self.t = np.abs(injection_changes(so)[ns]).sum(axis=1) / self.v_abs ** 2
        self.a1 = np.abs(so.E[ns]).sum(axis=1) / self.v_abs
        lo_e, hi_e = monomial_bounds(self.Hs[:m])
        lo_f, hi_f = monomial_bounds(self.Hs[m:])
        b2 = np.hypot(np.maximum(-lo_e, hi_e), np.maximum(-lo_f, hi_f))
        self.a2 = np.minimum(b2, np.hypot(self.Hs[:m], self.Hs[m:]).sum(axis=(1, 2))) / self.v_abs
        W = solve_blocked(so.lu, np.eye(2 * m))                       # J^-1; rows (e, f), columns (P, Q)
        a, b, c, d = W[:m, :m], W[:m, m:], W[m:, :m], W[m:, m:]
        s = a * a + b * b + c * c + d * d
        det = a * d - b * c
        self.T = np.sqrt(0.5 * (s + np.sqrt(np.maximum(s * s - 4.0 * det * det, 0.0))))     # 2x2 spectral norms

    def beta(self, z, scale: float = 1.0):
        t, a1, a2, sig = self.t * scale, self.a1 * scale, self.a2 * scale * scale, self.sig
        a = a1 + a2 + z
        with np.errstate(divide="ignore", invalid="ignore"):
            out = t * (a2 + z) + sig * (a2 + z) * (2.0 * a1 + a2 + z) + t * a * a + (sig + t) * a ** 3 / (1.0 - a)
        return np.where((sig + t) > 0.0, out, 0.0), a

    def _map(self, z, scale):
        b, a = self.beta(z, scale)
        if not np.all(np.isfinite(b)) or a.max() >= 1.0:
            return None, float(a.max())
        return (self.T @ (self.v_abs ** 2 * b)) / self.v_abs, float(a.max())

    def solve(self, scale: float = 1.0, max_iter: int = 2000, tol: float = 1e-13) -> InjectionFixedPoint:
        z = np.zeros(self.so.m)
        ok, a = False, float((self.a1 * scale).max())
        for it in range(1, max_iter + 1):
            new, a = self._map(z, scale)
            if new is None:
                break
            done = np.abs(new - z).max() <= tol * max(new.max(), 1e-300)
            z = new
            if done:
                z = z * (1.0 + INFLATION)                             # the iterates approach the fixed point from below
                img, a = self._map(z, scale)
                ok = img is not None and bool(np.all(img <= z))
                break
        s_abs = np.zeros(self.g.n_bus)
        s_abs[self.g.pvpq] = z * self.v_abs
        return InjectionFixedPoint(ok, it, s_abs, a)

    def modulus_allowance(self, fp: InjectionFixedPoint, scale: float = 1.0):
        """(e3 + |s_i|, e3, |s_i| bound) for every PQ bus."""
        e3 = modulus_higher_order(self.so, self.Hs, scale)
        s_mag = fp.s_abs[self.g.pq]
        return e3 + s_mag, e3, s_mag


def polar_contraction_radius(st, delta: float, h: float = 1e-3) -> float:
    """Perron root of sum_k |J_p^-1 (dJ_p/dx . X_k)| in polar coordinates (theta, |V|): the first-order
    analogue of G when the Jacobian depends on angle differences only. Diagnostic; not used in a proof."""
    g = st.grid
    n_th, n = len(g.pvpq), len(g.pvpq) + len(g.pq)
    row_p = np.full(g.n_bus, -1)
    row_p[g.pvpq] = np.arange(n_th)
    X = solve_blocked(st.lin.lu, (drive_matrix(st.src, n, row_p) * delta).toarray())
    va0, vm0 = np.angle(st.V_star), np.abs(st.V_star)
    G = np.zeros((n, n))
    for k in range(X.shape[1]):
        Jd = []
        for s in (h, -h):
            va, vm = va0.copy(), vm0.copy()
            va[g.pvpq] += s * X[:n_th, k]
            vm[g.pq] += s * X[n_th:, k]
            Jd.append(polar_jacobian(g.Ybus, vm * np.exp(1j * va), g.pvpq, g.pq))
        G += np.abs(solve_blocked(st.lin.lu, ((Jd[0] - Jd[1]) / (2.0 * h)).toarray()))
    x = np.ones(n)
    lam = 0.0
    for _ in range(300):
        y = G @ x
        lam = float(np.linalg.norm(y))
        x = y / lam
    return lam


def polar_true_radius(st, delta: float, n_vertices: int = 64) -> float:
    """Largest spectral radius of J_p^-1 (J_p(x* + X eps) - J_p(x*)) over a sample of vertices, in polar
    coordinates, with the first-order state X eps. No bound is taken: it shows how much the Jacobian
    really changes across the box, i.e. the room any contraction argument has."""
    g = st.grid
    n_th, n = len(g.pvpq), len(g.pvpq) + len(g.pq)
    row_p = np.full(g.n_bus, -1)
    row_p[g.pvpq] = np.arange(n_th)
    X = solve_blocked(st.lin.lu, (drive_matrix(st.src, n, row_p) * delta).toarray())
    va0, vm0 = np.angle(st.V_star), np.abs(st.V_star)
    J0 = polar_jacobian(g.Ybus, st.V_star, g.pvpq, g.pq)
    rng = np.random.default_rng(0)
    worst = 0.0
    for t in range(n_vertices):
        eps = np.ones(X.shape[1]) if t == 0 else np.sign(rng.standard_normal(X.shape[1]))
        va, vm = va0.copy(), vm0.copy()
        va[g.pvpq] += X[:n_th] @ eps
        vm[g.pq] += X[n_th:] @ eps
        A = (polar_jacobian(g.Ybus, vm * np.exp(1j * va), g.pvpq, g.pq) - J0).tocsr()
        op = LinearOperator((n, n), matvec=lambda x, A=A: st.lin.lu.solve(A @ x), dtype=float)
        try:
            lam = eigs(op, k=1, which="LM", tol=1e-6, maxiter=2000, return_eigenvectors=False)
            worst = max(worst, float(np.abs(lam[0])))
        except Exception:
            continue
    return worst


@dataclass
class VoltageEnclosure:
    """Bounds of |V_i| for every PQ bus over the whole uncertainty set."""

    lo: np.ndarray              # enclosure: second-order bound minus the allowance (NaN if no fixed point was found)
    hi: np.ndarray
    lo2: np.ndarray             # second-order bound alone: not an enclosure
    hi2: np.ndarray
    allowance: np.ndarray       # bound of everything above second order
    converged: bool             # a self-mapped set was found: lo, hi hold for every point of the set (exact arithmetic)
    rel_radius: float
    iterations: int


def voltage_enclosure(so: SecondOrder, first_order=None, max_k: int | None = None, node_budget: int = 0) -> VoltageEnclosure:
    """Enclosure of the voltage magnitudes of a network without PV buses over the box of ``so``.

    |V_i| - |V*_i| = c_i' eps + eps' H_i eps + rem. The quadratic form is bounded with RAA-K (blocks up to
    ``max_k`` symbols, all of them by default; every form is normalised first, see cdhim.knomial) and |rem|
    with ``InjectionEnclosure``.
    ``first_order``: (lo, hi) of c_i' eps over the uncertainty set when it is a subset of the box with a
    known support function (e.g. half disks); by default the box itself, -+ |c_i|_1.
    """
    enc = InjectionEnclosure(so)
    fp = enc.solve(1.0)
    allow, _, _ = enc.modulus_allowance(fp, 1.0)
    H = so.forms()
    k = H.shape[1] if max_k is None else max_k
    q = np.array([knomial_bound(h, k, node_budget, normalise=True)[:2] for h in H])
    if first_order is None:
        r = np.abs(so.c_lin).sum(axis=1)
        first_order = (-r, r)
    lo2, hi2 = so.vm + first_order[0] + q[:, 0], so.vm + first_order[1] + q[:, 1]
    if fp.converged:
        lo, hi = lo2 - allow, hi2 + allow
    else:
        lo, hi = np.full(len(lo2), np.nan), np.full(len(hi2), np.nan)
    return VoltageEnclosure(lo, hi, lo2, hi2, allow, fp.converged, fp.rel_radius, fp.iterations)
