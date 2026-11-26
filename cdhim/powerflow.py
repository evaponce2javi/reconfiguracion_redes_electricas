"""Sparse AC power flow in polar coordinates: mismatch, Jacobian, Newton and a batched solver.

State vector  x = [theta_pvpq ; |V|_pq],  mismatch  F = [Re(S_calc - S)_pvpq ; Im(S_calc - S)_pq].
The Jacobian uses true d/d|V| (not |V| d/d|V|), so the last block of J^-1 is directly d|V|/dS.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

SOLVE_BLOCK = 32                                   # right-hand sides per SuperLU call
SOLVE_THREADS = max(1, (os.cpu_count() or 2) // 2)


def solve_blocked(lu, rhs, block: int = SOLVE_BLOCK, threads: int | None = None, trans: str = "N"):
    """``lu.solve`` for a dense (n, m) right-hand side, in blocks of ``block`` columns.

    SuperLU sweeps every supernode over all columns of the right-hand side, so one call with hundreds
    of columns thrashes the cache (measured on the 5,227 x 5,227 PEGASE-2869 Jacobian: 1.2 ms per
    column with 512 columns, 0.2 ms with 16-32). The solve releases the GIL, so blocks run in threads.
    """
    rhs = np.asarray(rhs)
    if rhs.ndim == 1 or rhs.shape[1] <= block:
        return lu.solve(np.ascontiguousarray(rhs), trans=trans)
    threads = SOLVE_THREADS if threads is None else threads
    out = np.empty(rhs.shape, dtype=np.result_type(rhs.dtype, np.float64))
    starts = range(0, rhs.shape[1], block)

    def work(s):
        out[:, s:s + block] = lu.solve(np.ascontiguousarray(rhs[:, s:s + block]), trans=trans)

    if threads <= 1:
        for s in starts:
            work(s)
    else:
        with ThreadPoolExecutor(threads) as pool:
            list(pool.map(work, starts))
    return out


def mismatch(Ybus, V, Sbus):
    """Complex power mismatch S_calc(V) - S_specified; works for 1-D V or a (n_bus, batch) array."""
    return V * np.conj(Ybus @ V) - Sbus


def stack_mismatch(mis, pvpq, pq):
    return np.concatenate([mis[pvpq].real, mis[pq].imag], axis=0)


def polar_jacobian(Ybus, V, pvpq, pq):
    """Exact power-flow Jacobian dF/dx at complex voltages ``V`` (CSC)."""
    n = V.shape[0]
    idx = np.arange(n)
    Ibus = Ybus @ V
    diag_v = sp.csr_matrix((V, (idx, idx)), shape=(n, n))
    diag_i = sp.csr_matrix((Ibus, (idx, idx)), shape=(n, n))
    diag_vn = sp.csr_matrix((V / np.abs(V), (idx, idx)), shape=(n, n))
    ds_dvm = diag_v @ (Ybus @ diag_vn).conj() + diag_i.conj() @ diag_vn
    ds_dva = 1j * (diag_v @ (diag_i - Ybus @ diag_v).conj())
    ds_dva, ds_dvm = ds_dva.tocsr(), ds_dvm.tocsr()
    a_p, m_p = ds_dva[pvpq].tocsc(), ds_dvm[pvpq].tocsc()
    a_q, m_q = ds_dva[pq].tocsc(), ds_dvm[pq].tocsc()
    return sp.bmat(
        [[a_p[:, pvpq].real, m_p[:, pq].real], [a_q[:, pvpq].imag, m_q[:, pq].imag]], format="csc"
    )


def legacy_flat_jacobian(Ybus, v_flat, pvpq, pq):
    """Jacobian exactly as built by ``finalFinalCode.py`` (lines 104-131). Kept ONLY to reproduce the
    figures of the submitted manuscript.

    It differs from the exact Jacobian in three ways: (i) the off-diagonal entries of H = dP/dtheta
    carry the wrong sign, (ii) it is evaluated at a flat start (theta = 0, |V| = 1 at PQ buses) and
    (iii) the diagonal entries drop the injection terms (-Q_i, +P_i, ...).
    """
    n = v_flat.shape[0]
    idx = np.arange(n)
    dv = sp.csr_matrix((v_flat, (idx, idx)), shape=(n, n))
    G, B = Ybus.real.tocsr(), Ybus.imag.tocsr()
    b_diag = sp.csr_matrix((B.diagonal(), (idx, idx)), shape=(n, n))
    H = (dv @ (B - 2.0 * b_diag) @ dv).tocsr()   # off-diagonal +ViVjBij  (exact: -ViVjBij)
    N = (dv @ G @ dv).tocsr()
    M = -N
    L = (-(dv @ B @ dv)).tocsr()
    return sp.bmat(
        [[H[pvpq].tocsc()[:, pvpq], N[pvpq].tocsc()[:, pq]], [M[pq].tocsc()[:, pvpq], L[pq].tocsc()[:, pq]]],
        format="csc",
    )


def newton_pf(Ybus, Sbus, V_init, pv, pq, tol=1e-9, max_iter=30):
    """Full Newton-Raphson. Returns (V, converged, iterations, final mismatch norm)."""
    pvpq = np.concatenate([pv, pq])
    n_th = len(pvpq)
    V = np.array(V_init, dtype=complex)
    va, vm = np.angle(V), np.abs(V)
    norm = np.inf
    for it in range(max_iter + 1):
        F = stack_mismatch(mismatch(Ybus, V, Sbus), pvpq, pq)
        norm = float(np.max(np.abs(F))) if F.size else 0.0
        if not np.isfinite(norm):
            return V, False, it, norm
        if norm < tol:
            return V, True, it, norm
        if it == max_iter:
            break
        try:
            dx = splu(polar_jacobian(Ybus, V, pvpq, pq)).solve(F)
        except RuntimeError:          # exactly singular Jacobian
            return V, False, it, norm
        va[pvpq] -= dx[:n_th]
        vm[pq] -= dx[n_th:]
        if np.any(vm[pq] <= 0.0) or not np.all(np.isfinite(dx)):
            return V, False, it, norm
        V = vm * np.exp(1j * va)
    return V, False, max_iter, norm


@dataclass
class BatchStats:
    n: int = 0
    converged: int = 0
    newton_fallback: int = 0
    chord_iterations: int = 0

    def add(self, other: "BatchStats") -> None:
        self.n += other.n
        self.converged += other.converged
        self.newton_fallback += other.newton_fallback
        self.chord_iterations += other.chord_iterations


class BatchACSolver:
    """Solve many AC power flows that share topology and differ only in the injections.

    Every sample is an exact AC solution (mismatch below ``tol``). The iteration is a chord
    (fixed-Jacobian) Newton run on the whole batch with one sparse LU, warm-started at the base
    solution; columns that stall or diverge are re-solved one by one with full Newton-Raphson.
    A sample is reported as non-converged only if full Newton-Raphson also fails.

    Internally the buses are reordered as [other, pv, pq] so that the state update touches contiguous
    row blocks (views) instead of fancy-indexed copies.
    """

    def __init__(self, Ybus, pv, pq, V_base, lu=None, tol=1e-9, max_chord=40, max_newton=25, threads=None):
        self.Ybus = Ybus.tocsr()
        self.pv, self.pq = np.asarray(pv), np.asarray(pq)
        self.pvpq = np.concatenate([self.pv, self.pq])
        self.n_th = len(self.pvpq)
        self.V_base = np.asarray(V_base, dtype=complex)
        self.lu = lu if lu is not None else splu(polar_jacobian(self.Ybus, self.V_base, self.pvpq, self.pq))
        self.tol, self.max_chord, self.max_newton, self.threads = tol, max_chord, max_newton, threads
        n = self.Ybus.shape[0]
        rest = np.setdiff1d(np.arange(n), self.pvpq)
        self._order = np.concatenate([rest, self.pvpq])           # internal position -> bus
        self._r0 = len(rest)                                       # rows [r0:] are pv+pq
        self._q0 = len(rest) + len(self.pv)                        # rows [q0:] are pq
        self._Yp = self.Ybus[self._order][:, self._order].tocsr()
        self._va0 = np.angle(self.V_base[self._order])
        self._vm0 = np.abs(self.V_base[self._order])

    def solve(self, Sbus):
        """``Sbus``: (n_bus, batch) complex injections. Returns (V, converged mask, BatchStats)."""
        n_b = Sbus.shape[1]
        r0, q0, n_th = self._r0, self._q0, self.n_th
        S = Sbus[self._order]                                      # permuted copy, compacted with the active set
        va = np.repeat(self._va0[:, None], n_b, axis=1)
        vm = np.repeat(self._vm0[:, None], n_b, axis=1)
        V = vm * np.exp(1j * va)
        Vout = np.empty((len(self._order), n_b), dtype=complex)    # in permuted order
        conv = np.zeros(n_b, dtype=bool)
        cols = np.arange(n_b)                                      # original column of each active column
        prev = np.full(n_b, np.inf)
        fallback: list[np.ndarray] = []
        stats = BatchStats(n=n_b)
        F = np.empty((n_th + (len(self._order) - q0), n_b))
        for it in range(self.max_chord + 1):
            mis = self._Yp @ V
            np.conj(mis, out=mis)
            mis *= V
            mis -= S
            F = np.concatenate([mis[r0:].real, mis[q0:].imag], axis=0)
            norm = np.abs(F).max(axis=0)
            done = norm < self.tol
            bad = ~np.isfinite(norm) | (norm > 1e3) | ((it >= 6) & (norm > 0.7 * prev) & ~done)
            last = it == self.max_chord
            if done.any():
                Vout[:, cols[done]] = V[:, done]
                conv[cols[done]] = True
            leave = done | bad if not last else np.ones_like(done)
            if (bad | (last & ~done)).any():
                fallback.append(cols[(bad | last) & ~done])
            if leave.any():
                keep = ~leave
                cols, S, va, vm, F, norm = cols[keep], S[:, keep], va[:, keep], vm[:, keep], F[:, keep], norm[keep]
            if cols.size == 0:
                break
            prev = norm
            dx = solve_blocked(self.lu, F, threads=self.threads)
            va[r0:] -= dx[:n_th]
            vm[q0:] -= dx[n_th:]
            V = vm * np.exp(1j * va)
            stats.chord_iterations += cols.size
        V_full = np.empty_like(Vout)
        V_full[self._order] = Vout
        if fallback:
            for k in np.concatenate(fallback):
                v_k, ok, _, _ = newton_pf(self.Ybus, Sbus[:, k], self.V_base, self.pv, self.pq, self.tol, self.max_newton)
                stats.newton_fallback += 1
                if ok:
                    V_full[:, k] = v_k
                    conv[k] = True
        if not conv.all():
            V_full[:, ~conv] = np.nan
        stats.converged = int(conv.sum())
        return V_full, conv, stats
