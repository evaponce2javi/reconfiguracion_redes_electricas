"""Published disc bounds of the load-flow solution for networks with PQ buses only (any topology, shunts allowed).

Wang, Bernstein, Le Boudec, Paolone, "Explicit conditions on existence and uniqueness of load-flow solutions in
distribution networks", IEEE Trans. Smart Grid 9(2), 2018, Theorem 1. With w the zero-load voltages,
W = diag(w), xi(s) = || W^-1 Y_LL^-1 conj(W)^-1 diag(conj(s)) ||_inf and u_min = min_j |v^_j / w_j| for a known
solution (v^, s^): if xi(s^) < u_min^2 and

    Delta = (u_min - xi(s^)/u_min)^2 - 4 xi(s - s^) > 0,

the load flow for s has a unique solution in { v : |v_j - v^_j| <= rho |w_j| },
rho = ((u_min - xi(s^)/u_min) - sqrt(Delta)) / 2.

xi is monotone in |s - s^|, so evaluating it with the largest |s_k - s^_k| of an uncertainty set gives a disc
that holds for the whole set. This is the baseline for meshed networks, where the radial method does not apply.

``per_bus`` refines the uniform radius bus by bus: from v - v^ = Y_LL^-1 (conj(s/v) - conj(s^/v^)) and the disc,

    |v_j - v^_j| <= sum_l |Z_jl| (u_min |s_l - s^_l| + rho |s^_l|) / (|w_l| u_min (u_min - rho)).

It is the single-phase reading of Proposition 3.4 of C. Wang's thesis (EPFL 9456, 2019), re-derived here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import splu

from .network import Grid


def disc_radius(u_min: float, xi_hat: float, xi_delta: float):
    """rho of Theorem 1, or None if its conditions (9)-(10) fail."""
    if not xi_hat < u_min ** 2:
        return None
    a = u_min - xi_hat / u_min
    delta = a * a - 4.0 * xi_delta
    if not delta > 0.0:
        return None
    return 0.5 * (a - np.sqrt(delta))


@dataclass
class DiscBound:
    ok: bool                  # the conditions of the theorem hold for the whole set
    rho: float                # |v_j - v^_j| <= rho |w_j|
    radius: np.ndarray        # per PQ bus, uniform bound rho |w_j|            (NaN if not ok)
    radius_bus: np.ndarray    # per PQ bus, refined bound
    u_min: float
    xi_hat: float
    xi_delta: float


def disc_bound(grid: Grid, s_hat: np.ndarray, v_hat: np.ndarray, ds_abs: np.ndarray) -> DiscBound:
    """Disc around the solution ``v_hat`` for ``s_hat`` that contains the solution for every injection s with
    |s_k - s_hat_k| <= ds_abs_k. All arrays per internal bus; the result follows the order of ``grid.pq``."""
    if len(grid.pv):
        raise ValueError("the disc bound needs PQ buses only")
    L, ref = grid.pq, grid.ref
    Y = grid.Ybus.tocsc()
    lu = splu(Y[L][:, L].tocsc())
    Z = lu.solve(np.eye(len(L), dtype=complex))
    w = -Z @ (Y[L][:, ref] @ grid.V0[ref])                  # zero-load voltages (the paper has one slack bus)
    A = np.abs(Z) / (np.abs(w)[:, None] * np.abs(w)[None, :])
    xi_hat = float((A @ np.abs(s_hat[L])).max())
    xi_delta = float((A @ ds_abs[L]).max())
    u_min = float(np.abs(v_hat[L] / w).min())
    rho = disc_radius(u_min, xi_hat, xi_delta)
    nan = np.full(len(L), np.nan)
    if rho is None or not rho < u_min:
        return DiscBound(False, np.nan, nan, nan, u_min, xi_hat, xi_delta)
    per_bus = (np.abs(Z) / np.abs(w)[None, :]) @ (u_min * ds_abs[L] + rho * np.abs(s_hat[L])) / (u_min * (u_min - rho))
    return DiscBound(True, float(rho), rho * np.abs(w), np.minimum(per_bus, rho * np.abs(w)), u_min, xi_hat, xi_delta)
