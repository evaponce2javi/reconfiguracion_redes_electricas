"""B9 -- Branch current magnitudes |I| (both ends of every in-service branch) with the affine-forms method.

In transmission, what an operator watches on a branch is its loading, i.e. the current. The current at the
from end is I_f = Yf V. Its first-order affine form follows from the complex voltage change dV = E eps:
    d|I_f| = Re( conj(I_f0)/|I_f0| * (Yf E) ) eps.
The nonlinear remainder is sampled on the second-order state dV1 + dV2 (SecondOrder.state):
    r(eps) = |I_f0 + Yf (dV1 + dV2)| - |I_f0| - (linear part),
and its range is added to the radius, as for the voltages (B4, B5). Same for the to end with Yt.

Methods per branch end (interval = quantiles alpha/2, 1 - alpha/2): Bernstein + range, 3 sigma + range, the cumulant
method with Cornish-Fisher (linear) and Hong's PEM with cumulants on the exact AC power flow (as in B8).
Checks: AC Monte Carlo per branch end and for the whole system. Branch ends with |I_0| below 1e-6 p.u. are left out.

Cells: design point of the four systems and the two stressed cells of B8.

    python experiments/b9_branch_currents.py
Writes results/b9_branch_currents.csv
"""
from time import perf_counter

import numpy as np
import scipy.sparse as sp

from common import CFG, grid, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, loadability_limit, sample_eps, save_table
from cdhim.affine_forms import affine_radius, bernstein_radius, sigma
from cdhim.validation import VIOLATION_ATOL
from b8_ppf_baselines import CASES, K4_U, KURT_U, VAR_U, interval_from_moments

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
N_SUR = 1000
I_MIN = 1e-6


def current_model(st, so):
    g = st.grid
    Yb = sp.vstack([g.Yf, g.Yt]).tocsr()
    I0 = Yb @ st.V_star
    keep = np.abs(I0) > I_MIN
    Yb, I0 = Yb[keep], I0[keep]
    u = np.conj(I0) / np.abs(I0)
    C = np.real(u[:, None] * (Yb @ so.E))                  # (n_ends, N_u) first-order coefficients of |I|
    return Yb, I0, C, keep


def remainder_range(st, so, Yb, I0, C, n, seed, batch=250):
    rng = np.random.default_rng(seed)
    lo, hi = np.full(C.shape[0], np.inf), np.full(C.shape[0], -np.inf)
    for s in range(0, n, batch):
        eps = sample_eps(rng, st.src.n, min(batch, n - s))
        r = np.abs(I0[:, None] + Yb @ so.state(eps)) - np.abs(I0)[:, None] - C @ eps
        np.minimum(lo, r.min(axis=1), out=lo)
        np.maximum(hi, r.max(axis=1), out=hi)
    return lo, hi


def pem_currents(st, delta, Yb, I0, batch=250):
    """Hong 2m+1 PEM with cumulants (as b8_ppf_baselines.pem) for the current magnitudes."""
    m = st.src.n
    e = np.sqrt(KURT_U * VAR_U)
    p = 1.0 / (2.0 * KURT_U)
    K = np.zeros((Yb.shape[0], 4))
    for s in range(0, m, batch):
        cols = np.arange(s, min(s + batch, m))
        X = np.zeros((m, 2 * cols.size))
        X[cols, np.arange(cols.size)] = e
        X[cols, cols.size + np.arange(cols.size)] = -e
        V, conv, _ = st.solver.solve(st.injections(X, delta))
        y = np.abs(Yb @ V) - np.abs(I0)[:, None]
        y[:, ~conv] = 0.0
        a, b = y[:, :cols.size], y[:, cols.size:]
        m1, m2, m3, m4 = p * (a + b), p * (a ** 2 + b ** 2), p * (a ** 3 + b ** 3), p * (a ** 4 + b ** 4)
        K[:, 0] += m1.sum(axis=1)
        K[:, 1] += (m2 - m1 ** 2).sum(axis=1)
        K[:, 2] += (m3 - 3 * m2 * m1 + 2 * m1 ** 3).sum(axis=1)
        K[:, 3] += (m4 - 4 * m3 * m1 - 3 * m2 ** 2 + 12 * m2 * m1 ** 2 - 6 * m1 ** 4).sum(axis=1)
    var = K[:, 1]
    ok = var > 0
    g1 = np.divide(K[:, 2], var ** 1.5, out=np.zeros_like(var), where=ok)
    g2 = np.divide(K[:, 3], var ** 2, out=np.zeros_like(var), where=ok)
    return interval_from_moments(K[:, 0], var, g1, g2)


def run(task):
    case, lf, tag, delta, n = task
    t0 = perf_counter()
    seed = seed_for("b9", case, lf, delta)
    st = study(case, load_factor=lf, threads=worker_threads())
    t1 = perf_counter()
    so = SecondOrder(st, delta)
    Yb, I0, C, keep = current_model(st, so)
    lo_r, hi_r = remainder_range(st, so, Yb, I0, C, N_SUR, seed + 1)
    t_ours = perf_counter() - t1
    rb, s3 = bernstein_radius(C, ALPHA), 3.0 * sigma(C)
    k2, k4 = VAR_U * (C ** 2).sum(axis=1), K4_U * (C ** 4).sum(axis=1)
    g2 = np.divide(k4, k2 ** 2, out=np.zeros_like(k2), where=k2 > 0)
    t1 = perf_counter()
    pem = pem_currents(st, delta, Yb, I0)
    t_pem = perf_counter() - t1
    bounds = {"bern": (-rb + np.minimum(0, lo_r), rb + np.maximum(0, hi_r)),
              "sig3": (-s3 + np.minimum(0, lo_r), s3 + np.maximum(0, hi_r)),
              "affine": (-affine_radius(C) + np.minimum(0, lo_r), affine_radius(C) + np.maximum(0, hi_r)),
              "cumulant_cf": interval_from_moments(np.zeros_like(k2), k2, np.zeros_like(k2), g2),
              "pem_cf": pem,
              "bern_linear_only": (-rb, rb)}
    counts = {k: np.zeros(C.shape[0], np.int64) for k in bounds}
    anyc = {k: 0 for k in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    err_lin = err_so = 0.0
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        y = np.abs(Yb @ V[:, conv]) - np.abs(I0)[:, None]
        n_ac += int(conv.sum())
        if s == 0:                                         # accuracy of the two models on the first batch
            e1 = eps[:, conv]
            err_lin = float(np.abs(y - C @ e1).max())
            err_so = float(np.abs(y - (np.abs(I0[:, None] + Yb @ so.state(e1)) - np.abs(I0)[:, None])).max())
        for k, (lo, hi) in bounds.items():
            bad = (y < lo[:, None] - VIOLATION_ATOL) | (y > hi[:, None] + VIOLATION_ATOL)
            counts[k] += bad.sum(axis=1)
            anyc[k] += int(bad.any(axis=0).sum())
    rows = []
    live = affine_radius(C) > 1e-10
    wb = bounds["bern"][1] - bounds["bern"][0]
    for k, c in counts.items():
        lo, hi = bounds[k]
        rate = c / max(n_ac, 1)
        rows.append({"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "n_samples": n, "n_ac_not_converged": n - n_ac,
                     "seed": seed, "n_branch_ends": int(C.shape[0]), "method": k, "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()),
                     "n_over_alpha": int((rate > ALPHA).sum()), "any_end_pct": pct(anyc[k] / max(n_ac, 1)),
                     "width_over_bern_median": float(np.median((hi - lo)[live] / wb[live])),
                     "width_over_affine_median": float(np.median((hi - lo)[live] / (2 * affine_radius(C)[live]))),
                     "rel_halfwidth_median_pct": pct(np.median((0.5 * (hi - lo))[live] / np.abs(I0)[live])),
                     "max_err_linear_pu": err_lin, "max_err_second_order_pu": err_so,
                     "t_method_s": round(t_pem if k == "pem_cf" else (t_ours if k in ("bern", "sig3", "affine") else 0.0), 3)})
    print(f"[b9] {case:26s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s (ours {t_ours:.1f}s, PEM {t_pem:.1f}s)", flush=True)
    return rows


def main():
    tasks = []
    for case, d in CASES.items():
        tasks.append((case, 1.0, "design", d, EXP["n_ac"]))
        tasks.append((case, 1.2, "stress_delta", 0.30, EXP["n_ac"]))
        lam, _ = loadability_limit(grid(case), tol=CFG["powerflow"]["tol_pu"])
        tasks.append((case, round(EXP["near_limit_fraction"] * lam, 4), "near_limit", 0.15, EXP["n_ac"]))
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b9_branch_currents", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR,
                                                  "note": "targets: |I| at both ends of every in-service branch with |I_0| > 1e-6 p.u.; "
                                                          "bern_linear_only: without the remainder range"}))


if __name__ == "__main__":
    main()
