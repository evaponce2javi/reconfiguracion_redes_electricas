"""Worst-case load flow for radial feeders after Sou, Khong and Sandberg (IEEE TPWRS 41(1), 2026).

Reimplemented from the paper as the published baseline for radial networks. Model: balanced radial
network, no shunt elements, DistFlow equations

    nu  = nu0 1 + 2 R p + 2 X q + M l                                   (their Eq. 2a)
    l   = nu^-1 . |D (p + j q) - (D - I) Dz l|^2                        (2b)

with nu the squared voltages, l the squared line currents, D the descendant matrix and
R = D'Dr D, X = D'Dx D, M = D'Dr (I - 2D) Dr + D'Dx (I - 2D) Dx  (entries of M are <= 0).

Uncertainty: the injection of bus n is (p0_n + dp_n, q0_n + dq_n) with (dp_n, dq_n) in the half disk
dp_n >= 0, |dp_n + j dq_n| <= C_n.

Outputs, as in the paper:
  v_min   exact minimum voltages from ONE load flow at (dp, dq) = (0, -C), assuming the voltages are
          monotone in the injections (argued in the paper for typical data, not proved in general);
  v_max_L lower bound of the maximum voltages: largest value over load flows on the boundary (Eqs. 19-21);
  v_max_U upper bound of the maximum voltages from the fixed-point iteration of Algorithm 1.

The load flows are solved with the exact AC solver of this package, which is equivalent to DistFlow on
a radial network without shunts.

Also here: the rule the paper uses to add distributed generation to its benchmarks (Sec. VI-I) and its net
injection model (Eq. 43).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .network import Grid
from .powerflow import BatchACSolver


# Sec. VI-I: (F) forward and (R) reverse power flow; the paper uses rho_U = 0.8 in F and 0.5 in R
SCENARIOS = {"F": dict(rho_D=1.0, rho_G=0.5, pf=0.8, rho_U=0.8), "R": dict(rho_D=0.5, rho_G=1.0, pf=0.8, rho_U=0.5)}


def benchmark_ratings(grid: Grid) -> np.ndarray:
    """Inverter rating per internal bus [p.u.] after Sec. VI-I: a source at bus 1 + 10 k (bus 0 = reference), all
    with the same rating, the total equal to the total apparent power load (read here as |sum of the loads|).

    The non-reference buses are taken in the order of their labels, so for the benchmarks of the paper (reference
    = first bus) this is its rule; for other networks it is every tenth non-reference bus."""
    lab = np.asarray(grid.bus_label)
    order = np.argsort(lab)
    order = order[~np.isin(order, grid.ref)]                # paper bus n >= 1 = n-th non-reference bus
    dg = order[np.arange(0, len(order), 10)]
    rating = np.zeros(grid.n_bus)
    rating[dg] = abs(grid.s_load.sum()) / len(dg)
    return rating


def net_injections(grid: Grid, rating: np.ndarray, rho_D: float, rho_G: float, pf: float, rho_U: float):
    """Eq. (43): nominal complex injection per bus and capacity of the uncertain part, (s0, cap)."""
    s0 = rho_G * (pf + 1j * np.sqrt(1.0 - pf ** 2)) * (1.0 - rho_U) * rating - rho_D * grid.s_load
    return s0, rho_U * rating


@dataclass
class Feeder:
    """Tree data of a radial grid in the notation of the paper (non-reference buses in ``nodes`` order)."""

    nodes: np.ndarray        # internal bus index of every non-reference bus
    parent: np.ndarray       # internal bus index of its parent
    r: np.ndarray            # series resistance / reactance of the line entering each node [p.u.]
    x: np.ndarray
    D: np.ndarray            # D[n, m] = 1 if m is n or a descendant of n
    nu0: float

    @property
    def R(self):
        return self.D.T @ (self.r[:, None] * self.D)

    @property
    def X(self):
        return self.D.T @ (self.x[:, None] * self.D)

    @property
    def M(self):
        A = np.eye(len(self.r)) - 2.0 * self.D
        return self.D.T @ (self.r[:, None] * A * self.r[None, :]) + self.D.T @ (self.x[:, None] * A * self.x[None, :])


def feeder(grid: Grid, shunt_tol: float = 1e-8) -> Feeder:
    """Tree structure of ``grid``; raises if it is not a radial network without shunts and PV buses."""
    n = grid.n_bus
    if len(grid.pv) or len(grid.ref) != 1:
        raise ValueError("the DistFlow baseline needs a single slack bus and no voltage-controlled buses")
    Y = grid.Ybus.tocsr()
    if np.abs(np.asarray(Y.sum(axis=1)).ravel()).max() > shunt_tol * np.abs(Y).max():
        raise ValueError("the network has shunt elements or off-nominal taps: not a DistFlow model")
    edges = {(min(a, b), max(a, b)) for a, b in zip(grid.br_f, grid.br_t) if a != b}
    if len(edges) != n - 1:
        raise ValueError(f"{len(edges)} distinct branches for {n} buses: the network is not radial")
    nbr = [[] for _ in range(n)]
    for a, b in edges:
        nbr[a].append(b)
        nbr[b].append(a)
    root = int(grid.ref[0])
    parent = np.full(n, -1)
    order, seen = [root], {root}
    for u in order:                                       # breadth-first: parents come before children
        for v in nbr[u]:
            if v not in seen:
                seen.add(v)
                parent[v] = u
                order.append(v)
    if len(order) != n:
        raise ValueError("the network is not connected")
    nodes = np.array(order[1:])
    pos = np.full(n, -1)
    pos[nodes] = np.arange(n - 1)
    z = np.array([-1.0 / Y[parent[b], b] for b in nodes])
    D = np.eye(n - 1)
    for b in nodes[::-1]:                                 # children before parents
        if parent[b] != root:
            D[pos[parent[b]]] += D[pos[b]]
    return Feeder(nodes, parent[nodes], z.real, z.imag, D, float(np.abs(grid.V0[root]) ** 2))


def _half_disk_distance(a_lo, a_hi, b_lo, b_hi, rad):
    """Distance between the rectangle [a_lo, a_hi] x [b_lo, b_hi] and the half disk x >= 0, |x + jy| <= rad (Eq. 30)."""
    a = np.where(a_hi < 0, a_hi, np.where(a_lo > 0, a_lo, 0.0))          # point of the rectangle closest to the half disk
    b = np.clip(0.0, b_lo, b_hi)
    right = np.maximum(np.hypot(a, b) - rad, 0.0)
    left = np.hypot(a, np.maximum(np.abs(b) - rad, 0.0))
    return np.where(a >= 0, right, left)


def _half_disk_reach(a, b, rad):
    """Largest |(a + jb) + (x + jy)| over the half disk x >= 0, |x + jy| <= rad (Eq. 37)."""
    return np.where(a > 0, np.hypot(a, b) + rad, np.where(b >= 0, np.hypot(a, b + rad), np.hypot(a, b - rad)))


_L_CAP = 1e100          # cap of the squared-current upper bounds during the iteration [p.u.]


def _rel_change(new, old):
    den = np.linalg.norm(new)
    return float(np.linalg.norm(new - old) / den) if den > 0.0 else float(np.linalg.norm(old))


@dataclass
class WorstCase:
    v_min: np.ndarray          # exact under monotonicity (one load flow)
    v_max_L: np.ndarray        # attained by a load flow: v_max >= v_max_L
    v_max_U: np.ndarray        # Algorithm 1: v_max <= v_max_U
    iterations: int
    converged: bool
    nodes: np.ndarray


def worst_case_load_flow(grid: Grid, s0: np.ndarray, cap: np.ndarray, n_samples: int = 8, tol: float = 1e-6,
                         max_iter: int = 500, current_start: float = 3.0, current_start_amp: float | None = None) -> WorstCase:
    """Extremal voltages for injections s0 + (dp + j dq), (dp_n, dq_n) in the half disk of radius cap_n.

    ``s0``, ``cap``: complex nominal injection and capacity of the uncertain part, per internal bus.
    Initial upper bound of the line currents: the paper uses "a large number (e.g. 5,000 A)"; pass it as
    ``current_start_amp``. Otherwise the bound starts at ``current_start`` times the largest current of the
    sampled load flows. The update of a line uses the bounds of its strict descendants only, so the iteration
    sweeps from the leaves to the root and forgets the start; on the way the bounds of the upstream lines grow
    as the square of the previous ones, and are capped at ``_L_CAP`` here to keep the arithmetic finite.
    """
    fd = feeder(grid)
    nodes = fd.nodes
    solver = BatchACSolver(grid.Ybus, grid.pv, grid.pq, grid.V0)
    k = np.arange(n_samples + 1)
    S = np.repeat(s0[:, None], n_samples + 1, axis=1) + cap[:, None] * (np.sin(k * np.pi / n_samples) + 1j * np.cos(k * np.pi / n_samples))[None, :]
    V, conv, _ = solver.solve(S)
    if not conv.all():
        raise RuntimeError("a boundary load flow did not converge")
    vm = np.abs(V[nodes])
    v_min = vm[:, n_samples]                               # (dp, dq) = (0, -C): Eq. (7)
    v_max_L = vm[:, 1:].max(axis=1)                        # Eq. (21)
    cur = np.abs((V[fd.parent] - V[nodes]) / (fd.r + 1j * fd.x)[:, None]) ** 2      # squared line currents of the samples

    D = fd.D
    DmI = D - np.eye(len(nodes))
    R, X, M = fd.R, fd.X, fd.M
    Z = np.hypot(R, X)
    p0, q0, C = s0[nodes].real, s0[nodes].imag, cap[nodes]
    Dp, Dq, DC = D @ p0, D @ q0, D @ C
    base = fd.nu0 + 2.0 * R @ p0 + 2.0 * X @ q0 + 2.0 * Z @ C
    nu_U = base.copy()                                     # Eq. (39a)
    l_L = np.zeros(len(nodes))
    if current_start_amp is None:
        l_start = (current_start ** 2) * cur.max()
    else:
        i_base = grid.base_mva * 1e3 / (np.sqrt(3.0) * float(grid.base_kv[grid.ref[0]]))        # [A]
        l_start = (current_start_amp / i_base) ** 2
    l_U = np.full(len(nodes), min(l_start, _L_CAP))
    ok = False
    for it in range(1, max_iter + 1):
        sr_L, sr_U = DmI @ (fd.r * l_L), DmI @ (fd.r * l_U)
        sx_L, sx_U = DmI @ (fd.x * l_L), DmI @ (fd.x * l_U)
        nu_new = base + M @ l_L                            # Eq. (25)
        J = _half_disk_distance(-Dp + sr_L, -Dp + sr_U, -Dq + sx_L, -Dq + sx_U, DC)                 # Eqs. (28)-(30)
        lL_new = J ** 2 / nu_U                             # Eq. (31)
        reach = np.zeros(len(nodes))
        for a in (Dp - sr_L, Dp - sr_U):                   # Eqs. (36)-(38)
            for b in (Dq - sx_L, Dq - sx_U):
                reach = np.maximum(reach, _half_disk_reach(a, b, DC))
        lU_new = np.minimum(reach ** 2 / v_min ** 2, _L_CAP)
        change = max(_rel_change(nu_new, nu_U), _rel_change(lU_new, l_U), _rel_change(lL_new, l_L))
        nu_U, l_L, l_U = nu_new, lL_new, lU_new
        if not np.isfinite(change):
            break
        if change <= tol and l_U.max() < _L_CAP:
            ok = True
            break
    return WorstCase(v_min, v_max_L, np.sqrt(nu_U), it, ok, nodes)
