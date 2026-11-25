"""Uncertain injections in half disks: (dp_n, dq_n) with dp_n >= 0 and |dp_n + j dq_n| <= C_n at every source bus.

This is the set of Sou, Khong and Sandberg (2026) for inverter-based generation. The affine machinery of this
package works on boxes, so every half disk is covered by its bounding box,

    dp_n = C_n/2 (1 + eps_p),   dq_n = C_n eps_q,     eps in [-1, 1]^2,

with the P symbols of all sources first and their Q symbols after them. The first-order term is then
maximised over the half disks themselves (closed form); the higher orders are bounded over the boxes.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import splu

from .network import Grid
from .powerflow import newton_pf, polar_jacobian
from .study import Study, custom_study


def half_disk_study(grid: Grid, s0: np.ndarray, cap: np.ndarray, **kw) -> Study:
    """Study centred on the bounding boxes of the half disks of radius ``cap`` added to the injections ``s0``."""
    dg = np.flatnonzero(cap > 0)
    n = len(dg)
    centre = np.asarray(s0, dtype=complex).copy()
    centre[dg] += cap[dg] / 2.0
    return custom_study(grid, centre, np.concatenate([dg, dg]), np.r_[np.zeros(n, bool), np.ones(n, bool)],
                        np.r_[cap[dg] / 2.0, cap[dg]], **kw)


def first_order_range(C: np.ndarray):
    """(lo, hi) of c' eps over the product of half disks, for coefficients ``C`` (rows; P symbols then Q symbols).

    With u = dp/C, w = dq/C in the unit half disk, c_p eps_p + c_q eps_q = -c_p + 2 c_p u + c_q w, and
    max (a u + b w) = hypot(max(a, 0), b) over the unit half disk.
    """
    n = C.shape[1] // 2
    cp, cq = C[:, :n], C[:, n:]
    hi = (-cp + np.hypot(np.maximum(2.0 * cp, 0.0), cq)).sum(axis=1)
    lo = (-cp - np.hypot(np.maximum(-2.0 * cp, 0.0), cq)).sum(axis=1)
    return lo, hi


def sample_half_disks(rng: np.random.Generator, n_sources: int, n: int, on_arc: float = 0.5) -> np.ndarray:
    """(n_sources, n) complex points u + j w of the unit half disk; the first ``on_arc`` share lies on the arc."""
    ang = rng.uniform(-np.pi / 2, np.pi / 2, (n_sources, n))
    rad = np.sqrt(rng.uniform(0.0, 1.0, (n_sources, n)))
    rad[:, :int(on_arc * n)] = 1.0
    return rad * (np.cos(ang) + 1j * np.sin(ang))


def attained_extreme(grid: Grid, s0: np.ndarray, cap: np.ndarray, sign: int, n_iter: int = 15, tol: float = 1e-9) -> np.ndarray:
    """Best value of |V_i| found for every PQ bus over the product of half disks (see ``attained_points``)."""
    return attained_points(grid, s0, cap, sign, n_iter, tol)[0]


def attained_points(grid: Grid, s0: np.ndarray, cap: np.ndarray, sign: int, n_iter: int = 15, tol: float = 1e-9):
    """Best value of |V_i| found for every PQ bus over the product of half disks: maximum (sign=+1) or minimum (-1),
    and the point that attains it, as u + j w of the unit half disk of every source: (n_pq,), (n_pq, n_sources).

    Repeated linear maximisation: solve the power flow, read the sensitivities of |V_i| to the injections of
    the sources at that point (one transposed solve), move every source to the point of its half disk that
    is best for those sensitivities, and repeat while the value improves. Every value is attained by a feasible
    injection, so it bounds the true extreme from inside; it is a local search, not a proof of optimality.
    """
    dg = np.flatnonzero(cap > 0)
    n_th = len(grid.pvpq)
    row_p = np.full(grid.n_bus, -1)
    row_p[grid.pvpq] = np.arange(n_th)
    row_q = np.full(grid.n_bus, -1)
    row_q[grid.pq] = n_th + np.arange(len(grid.pq))
    centre = np.asarray(s0, dtype=complex).copy()
    centre[dg] += cap[dg] / 2.0
    V0, ok, _, _ = newton_pf(grid.Ybus, centre, grid.V0, grid.pv, grid.pq, tol=tol, max_iter=40)
    if not ok:
        raise RuntimeError("no AC solution at the centre of the uncertainty set")
    lu0 = splu(polar_jacobian(grid.Ybus, V0, grid.pvpq, grid.pq).tocsc())
    best = np.full(len(grid.pq), -np.inf)
    point = np.zeros((len(grid.pq), len(dg)), complex)
    for i in range(len(grid.pq)):
        e = np.zeros(n_th + len(grid.pq))
        e[n_th + i] = 1.0
        y = lu0.solve(e, trans="T")
        V = V0
        for _ in range(n_iter):
            a, b = sign * y[row_p[dg]], sign * y[row_q[dg]]
            ap = np.maximum(a, 0.0)
            nrm = np.hypot(ap, b)
            u = np.divide(ap, nrm, out=np.zeros_like(ap), where=nrm > 0)
            w = np.divide(b, nrm, out=np.zeros_like(b), where=nrm > 0)
            s = np.asarray(s0, dtype=complex).copy()
            s[dg] += cap[dg] * (u + 1j * w)
            V, ok, _, _ = newton_pf(grid.Ybus, s, V, grid.pv, grid.pq, tol=tol, max_iter=40)
            if not ok:
                break
            val = sign * abs(V[grid.pq[i]])
            if val <= best[i] + 1e-13:
                break
            best[i] = val
            point[i] = u + 1j * w
            y = splu(polar_jacobian(grid.Ybus, V, grid.pvpq, grid.pq).tocsc()).solve(e, trans="T")
    return sign * best, point
