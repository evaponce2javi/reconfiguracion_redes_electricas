"""A linearised uncertainty study: operating point -> sparse LU -> affine coefficients -> CD-HIM."""
from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from .core import CDHIMResult, MethodConfig, run_cdhim
from .network import Grid
from .powerflow import BatchACSolver, legacy_flat_jacobian, newton_pf, polar_jacobian, solve_blocked
from .uncertainty import SourceCovariance


@dataclass
class Linearization:
    J: object
    lu: object
    point: np.ndarray                 # complex voltages where J was evaluated
    jacobian: str                     # "exact" | "legacy_signbug"
    t_jacobian: float
    t_lu: float
    nnz_J: int
    nnz_L: int
    nnz_U: int
    permc_spec: str


def linearize(grid: Grid, V_point, jacobian: str = "exact", permc_spec: str = "COLAMD") -> Linearization:
    """Build the Jacobian at ``V_point`` and factorise it with SuperLU (``scipy.sparse.linalg.splu``)."""
    t0 = perf_counter()
    if jacobian == "exact":
        J = polar_jacobian(grid.Ybus, V_point, grid.pvpq, grid.pq)
    elif jacobian == "legacy_signbug":
        J = legacy_flat_jacobian(grid.Ybus, np.abs(V_point), grid.pvpq, grid.pq)
    else:
        raise ValueError("jacobian must be 'exact' or 'legacy_signbug'")
    t1 = perf_counter()
    lu = splu(J, permc_spec=permc_spec)
    t2 = perf_counter()
    return Linearization(J, lu, np.asarray(V_point), jacobian, t1 - t0, t2 - t1, J.nnz, lu.L.nnz, lu.U.nnz, permc_spec)


@dataclass
class Sources:
    """Uncertain injections. ``unit[k]`` is the signed injection change per unit eps at delta = 1:
    a load P_d that grows by delta*P_d*eps lowers the injection, hence unit = -P_d."""

    bus: np.ndarray        # internal bus index
    is_q: np.ndarray       # False: active power injection, True: reactive power injection
    unit: np.ndarray       # [p.u.]
    row: np.ndarray        # row of the Jacobian that the injection drives
    Z: object = None       # (n_injections, n_symbols) sparse aggregation; None: one symbol per uncertain injection
    G: object = None       # (n_bus, n_symbols) sparse: generator active-power response per unit eps [p.u.]; None: slack only

    @property
    def n(self) -> int:
        """Number of noise symbols."""
        return len(self.bus) if self.Z is None else self.Z.shape[1]

    @property
    def n_injections(self) -> int:
        return len(self.bus)

    def expand(self, eps: np.ndarray) -> np.ndarray:
        """Noise symbols (n, batch) -> relative change of every uncertain injection (n_injections, batch)."""
        return eps if self.Z is None else self.Z @ eps


def build_sources(grid: Grid, load_factor: float = 1.0, min_pu: float = 1e-6) -> Sources:
    """One P symbol per non-slack bus with demand and one independent Q symbol per PQ bus with demand."""
    n_th = len(grid.pvpq)
    row_p = np.full(grid.n_bus, -1)
    row_p[grid.pvpq] = np.arange(n_th)
    row_q = np.full(grid.n_bus, -1)
    row_q[grid.pq] = n_th + np.arange(len(grid.pq))
    pd, qd = load_factor * grid.s_load.real, load_factor * grid.s_load.imag
    bp = np.flatnonzero((np.abs(pd) > min_pu) & (row_p >= 0))
    bq = np.flatnonzero((np.abs(qd) > min_pu) & (row_q >= 0))
    return Sources(np.concatenate([bp, bq]), np.concatenate([np.zeros(len(bp), bool), np.ones(len(bq), bool)]),
                   np.concatenate([-pd[bp], -qd[bq]]), np.concatenate([row_p[bp], row_q[bq]]))


def zone_aggregation(src: Sources, zone_of_bus, couple_pq: bool = True):
    """Aggregation matrix Z: every uncertain demand of a zone moves with the same noise symbol.

    ``couple_pq=True`` gives one symbol per zone (P and Q together, constant power factor);
    ``False`` gives one symbol for the active and one for the reactive demand of each zone.
    Zones without uncertain demand get no symbol.
    """
    zone = np.asarray(zone_of_bus)[src.bus]
    key = zone if couple_pq else 2 * zone + src.is_q
    _, col = np.unique(key, return_inverse=True)
    n = len(col)
    return sp.csr_matrix((np.ones(n), (np.arange(n), col)), shape=(n, int(col.max()) + 1))


def participation_factors(grid: Grid) -> np.ndarray:
    """Share of every bus in covering an active-power imbalance: proportional to its scheduled generation
    (PV and slack buses with positive output). Sums to one."""
    alpha = np.zeros(grid.n_bus)
    gen = np.concatenate([grid.pv, grid.ref])
    alpha[gen] = np.maximum(grid.s_gen.real[gen], 0.0)
    if alpha.sum() <= 0.0:                                  # no dispatched generation data: the slack takes it all
        alpha[grid.ref] = 1.0
    return alpha / alpha.sum()


def generator_response(src: Sources, participation) -> sp.csr_matrix:
    """Affine redispatch G: the generators cover the active-power change of every symbol in fixed shares.

    A symbol that changes the total demand by t (p.u. per unit eps) raises the generation of bus g by
    alpha_g * t, so the specified injections stay balanced and the slack bus only supplies its own share
    and the change in losses.
    """
    alpha = np.asarray(participation, dtype=float)
    p_unit = np.where(src.is_q, 0.0, src.unit)             # active-power injection change of each uncertain injection
    t = p_unit if src.Z is None else src.Z.T @ p_unit      # the same, per symbol
    g = np.flatnonzero(alpha)
    k = np.flatnonzero(t)
    return sp.csr_matrix((-np.outer(alpha[g], t[k]).ravel(), (np.repeat(g, len(k)), np.tile(k, len(g)))),
                         shape=(len(alpha), len(t)))


def drive_matrix(src: Sources, n_rows: int, row_p_of_bus: np.ndarray, rows: np.ndarray | None = None):
    """(n_rows, n_symbols) right-hand side of the linearised power-flow equations per unit eps at delta = 1.

    ``rows``: equation row driven by each uncertain injection (default: ``src.row``, polar Jacobian);
    ``row_p_of_bus``: row of the active-power equation of every bus, -1 where there is none (slack).
    """
    rows = src.row if rows is None else rows
    D = sp.csr_matrix((src.unit, (rows, np.arange(src.n_injections))), shape=(n_rows, src.n_injections))
    if src.Z is not None:
        D = D @ src.Z
    if src.G is not None:
        b = np.flatnonzero(row_p_of_bus >= 0)
        D = D + sp.csr_matrix((np.ones(len(b)), (row_p_of_bus[b], b)), shape=(n_rows, len(row_p_of_bus))) @ src.G
    return D.tocsc()


def voltage_coefficients(lin: Linearization, grid: Grid, src: Sources, chunk: int = 512, threads=None) -> np.ndarray:
    """Affine coefficients at delta = 1: C1[i, k] = d|V_i| / d eps_k for every PQ bus i.

    Solves J X = E for the N_u columns that carry a source and keeps only the |V| rows; J^-1 is never
    formed. With the legacy Jacobian the unknown is d|V|/|V|, which equals d|V| at the flat start.
    """
    n = lin.J.shape[0]
    n_th = len(grid.pvpq)
    if src.Z is not None or src.G is not None:     # aggregated symbols: one solve per symbol, not per injection
        row_p = np.full(grid.n_bus, -1)
        row_p[grid.pvpq] = np.arange(n_th)
        drive = drive_matrix(src, n, row_p)
        C1 = np.empty((len(grid.pq), src.n))
        for s in range(0, src.n, chunk):
            C1[:, s:s + chunk] = solve_blocked(lin.lu, drive[:, s:s + chunk].toarray(), threads=threads)[n_th:]
        return C1
    C1 = np.empty((len(grid.pq), src.n))
    for s in range(0, src.n, chunk):
        cols = np.arange(s, min(s + chunk, src.n))
        rhs = np.zeros((n, len(cols)))
        rhs[src.row[cols], np.arange(len(cols))] = src.unit[cols]
        C1[:, cols] = solve_blocked(lin.lu, rhs, threads=threads)[n_th:]
    return C1


@dataclass
class Study:
    """Everything that depends on the operating point but not on delta."""

    grid: Grid
    load_factor: float
    V_star: np.ndarray             # converged AC solution at this load factor (centre of every interval)
    sbus: np.ndarray
    lin: Linearization
    src: Sources
    C1: np.ndarray                 # (n_pq, N_u) coefficients at delta = 1
    fi: np.ndarray
    ti: np.ndarray
    branch_idx: np.ndarray
    solver: BatchACSolver
    timings: dict = field(default_factory=dict)

    @property
    def vm_star(self) -> np.ndarray:
        return np.abs(self.V_star[self.grid.pq])

    def coefficients(self, delta: float) -> np.ndarray:
        return self.C1 * delta

    def run(self, delta: float, cfg: MethodConfig = MethodConfig(), sigma_eps: SourceCovariance | None = None) -> CDHIMResult:
        return run_cdhim(self.coefficients(delta), self.fi, self.ti, cfg, sigma_eps)

    def injections(self, eps: np.ndarray, delta: float) -> np.ndarray:
        """Specified complex injections for a batch of noise samples ``eps`` (N_u, batch)."""
        S = np.repeat(self.sbus[:, None], eps.shape[1], axis=1)
        d = delta * self.src.unit[:, None] * self.src.expand(eps)
        p, q = ~self.src.is_q, self.src.is_q
        S[self.src.bus[p]] += d[p]
        S[self.src.bus[q]] += 1j * d[q]
        if self.src.G is not None:
            S += delta * (self.src.G @ eps)
        return S


def solve_operating_point(grid: Grid, load_factor: float, V_init=None, tol: float = 1e-9, step: float = 0.1):
    """Converged AC solution at ``load_factor``. Returns (V, ok).

    With ``V_init`` a single warm-started Newton solve is made; without it the load factor is marched
    from 1.0 in steps of ``step`` (continuation), which keeps Newton on the operating branch when the
    target is far from the base case.
    """
    if V_init is not None:
        V, ok, _, _ = newton_pf(grid.Ybus, grid.sbus(load_factor), V_init, grid.pv, grid.pq, tol=tol, max_iter=40)
        return V, ok
    V, lam = grid.V0, 1.0
    n_steps = max(1, int(np.ceil(abs(load_factor - 1.0) / step - 1e-9)))
    for target in np.linspace(1.0, load_factor, n_steps + 1)[1:]:
        V, ok, _, _ = newton_pf(grid.Ybus, grid.sbus(float(target)), V, grid.pv, grid.pq, tol=tol, max_iter=40)
        if not ok:
            return V, False
    return V, True


def make_study(grid: Grid, load_factor: float = 1.0, point: str = "converged", jacobian: str = "exact",
               permc_spec: str = "COLAMD", V_star=None, tol: float = 1e-9, max_chord: int = 40, threads=None,
               min_injection_pu: float = 1e-6, zones=None, couple_pq: bool = True, participation=None) -> Study:
    """``point``: "converged" linearises at the AC solution (what the manuscript states);
    "flat" linearises at theta = 0, |V| = 1 (what the original code did).

    ``zones``: zone label of every bus. If given, all the uncertain demands of a zone share one noise
    symbol (see ``zone_aggregation``); otherwise every P_d and every Q_d has its own symbol.

    ``participation``: share of every bus in covering an active-power imbalance (sums to one, see
    ``participation_factors``). If given, generation follows the demand change of every symbol in those
    shares; otherwise the slack bus absorbs the whole imbalance, as in the manuscript."""
    timings = {}
    if V_star is None:
        if load_factor == 1.0:
            V_star = grid.V0
        else:
            V_star, ok = solve_operating_point(grid, load_factor, tol=tol)
            if not ok:
                raise RuntimeError(f"no AC solution at load factor {load_factor}")
    if point not in ("converged", "flat"):
        raise ValueError("point must be 'converged' or 'flat'")
    if jacobian == "legacy_signbug" and point != "flat":
        raise ValueError("the legacy Jacobian is only defined at the flat start (point='flat')")
    V_lin = V_star if point == "converged" else grid.flat_voltage()
    lin = linearize(grid, V_lin, jacobian, permc_spec)
    src = build_sources(grid, load_factor, min_injection_pu)
    if zones is not None:
        src.Z = zone_aggregation(src, zones, couple_pq)
    if participation is not None:
        src.G = generator_response(src, participation)
    t0 = perf_counter()
    C1 = voltage_coefficients(lin, grid, src, threads=threads)
    timings.update(jacobian=lin.t_jacobian, lu=lin.t_lu, solves=perf_counter() - t0)
    fi, ti, bidx = grid.pq_pq_branches()
    exact_at_star = lin.lu if (point == "converged" and jacobian == "exact") else None
    solver = BatchACSolver(grid.Ybus, grid.pv, grid.pq, V_star, lu=exact_at_star, tol=tol, max_chord=max_chord, threads=threads)
    return Study(grid, load_factor, V_star, grid.sbus(load_factor), lin, src, C1, fi, ti, bidx, solver, timings)


def custom_study(grid: Grid, sbus: np.ndarray, bus, is_q, unit, Z=None, permc_spec: str = "COLAMD", tol: float = 1e-9,
                 max_chord: int = 40, threads=None) -> Study:
    """Study around the AC solution for the specified injections ``sbus`` with explicitly given uncertain injections.

    Injection k is the active (``is_q[k]`` False) or reactive power injected at ``bus[k]`` and changes by
    ``unit[k]`` p.u. per unit of its noise symbol, so ``sbus`` is the centre of the uncertainty box and
    ``unit`` its half-widths. Use it for sources that are not demands with a symmetric relative range, e.g.
    distributed generation between zero and its rating. ``Study.load_factor`` is NaN for these studies.
    """
    bus, is_q, unit = np.asarray(bus, dtype=int), np.asarray(is_q, dtype=bool), np.asarray(unit, dtype=float)
    n_th = len(grid.pvpq)
    row_p = np.full(grid.n_bus, -1)
    row_p[grid.pvpq] = np.arange(n_th)
    row_q = np.full(grid.n_bus, -1)
    row_q[grid.pq] = n_th + np.arange(len(grid.pq))
    row = np.where(is_q, row_q[bus], row_p[bus])
    if (row < 0).any():
        raise ValueError("an uncertain injection sits on a bus without that power-flow equation (slack, or Q at a PV bus)")
    V_star, ok, _, _ = newton_pf(grid.Ybus, sbus, grid.V0, grid.pv, grid.pq, tol=tol, max_iter=40)
    if not ok:
        raise RuntimeError("no AC solution at the centre of the uncertainty box")
    lin = linearize(grid, V_star, "exact", permc_spec)
    src = Sources(bus, is_q, unit, row, Z)
    t0 = perf_counter()
    C1 = voltage_coefficients(lin, grid, src, threads=threads)
    timings = dict(jacobian=lin.t_jacobian, lu=lin.t_lu, solves=perf_counter() - t0)
    fi, ti, bidx = grid.pq_pq_branches()
    solver = BatchACSolver(grid.Ybus, grid.pv, grid.pq, V_star, lu=lin.lu, tol=tol, max_chord=max_chord, threads=threads)
    return Study(grid, float("nan"), V_star, np.asarray(sbus, dtype=complex), lin, src, C1, fi, ti, bidx, solver, timings)


def loadability_limit(grid: Grid, lo: float = 1.0, hi: float = 8.0, tol_lambda: float = 2e-3, tol: float = 1e-9):
    """Largest load factor with an AC solution reachable by continuation from the base case (bisection)."""
    V = grid.V0
    step = 0.1
    lam = lo
    while lam + step <= hi:                      # march up until Newton fails
        V_new, ok = solve_operating_point(grid, lam + step, V, tol)
        if not ok:
            break
        lam, V = lam + step, V_new
    a, b = lam, min(lam + step, hi)
    while b - a > tol_lambda:
        mid = 0.5 * (a + b)
        V_new, ok = solve_operating_point(grid, mid, V, tol)
        if ok:
            a, V = mid, V_new
        else:
            b = mid
    return a, V
