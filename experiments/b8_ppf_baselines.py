"""B8 -- Published probabilistic power flow methods against the affine-forms method.

Baselines (per target: bus |V_i| and branch drop |V_i| - |V_j|; interval = quantiles alpha/2 and 1 - alpha/2):
  cumulant_cf   Cumulant method on the linear model (Zhang & Lee 2004 and successors): output cumulants
                kappa_n = sum_k c_k^n kappa_n(eps), U[-1,1] inputs (kappa_2 = 1/3, kappa_4 = -2/15), quantiles by the
                Cornish-Fisher expansion to fourth order. Linear, no distribution assumption beyond the inputs' cumulants.
  pem_cf        Hong's 2m+1 point estimate method (Hong 1998) on the exact AC power flow: two AC solutions per
                source at eps_k = xi_{k,1,2} sigma_eps; the cumulants of each source's contribution are added
                (independent sources) and the quantiles come from Cornish-Fisher. Captures the nonlinearity source by
                source, not the interactions between sources.
Ours (B5): bern (Bernstein radius + sampled range of the quadratic term), sig3 (3 sigma + range).

Checks: AC Monte Carlo per target and for the whole system, widths, and cost (PEM: 2m AC solutions).
Design point of the four systems and two stressed cells of each (load 1.2 with delta 30 %, and 95 % of the
loadability limit with delta 15 %).

    python experiments/b8_ppf_baselines.py
Writes results/b8_ppf_baselines.csv
"""
from time import perf_counter

import numpy as np
from scipy.special import ndtri

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, loadability_limit, sample_eps, save_table
from cdhim.affine_forms import affine_radius, bernstein_radius, bus_targets, drop_targets, quadratic_ranges, sigma
from cdhim.validation import VIOLATION_ATOL

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
N_SUR = 1000
VAR_U, K4_U = 1.0 / 3.0, -2.0 / 15.0             # variance and fourth cumulant of U[-1, 1]
KURT_U = 1.8                                     # its kurtosis (lambda_4); skewness 0
CASES = {"case33bw": 0.15, "pglib_opf_case2869_pegase": 0.05, "case118": 0.05, "case1354pegase": 0.05}


def cornish_fisher(p, g1, g2):
    """Standardised quantile at probability p for skewness g1 and excess kurtosis g2 (fourth-order CF)."""
    z = ndtri(p)
    return z + (z * z - 1) * g1 / 6 + (z ** 3 - 3 * z) * g2 / 24 - (2 * z ** 3 - 5 * z) * g1 ** 2 / 36


def interval_from_moments(mean, var, g1, g2):
    s = np.sqrt(np.maximum(var, 0.0))
    return mean + s * cornish_fisher(ALPHA / 2, g1, g2), mean + s * cornish_fisher(1 - ALPHA / 2, g1, g2)


def cumulant_cf(Ct):
    k2 = VAR_U * (Ct ** 2).sum(axis=1)
    k4 = K4_U * (Ct ** 4).sum(axis=1)
    g2 = np.divide(k4, k2 ** 2, out=np.zeros_like(k2), where=k2 > 0)
    return interval_from_moments(np.zeros_like(k2), k2, np.zeros_like(k2), g2)


def pem(st, delta, Ws, batch=250):
    """Hong's 2m+1 point estimate method on the AC power flow for U[-1,1] symbols, combined with cumulants.

    Each source k is represented by three points of its own law, eps_k in {+e, 0, -e} with e = sqrt(lambda_4 var),
    probabilities 1/(2 lambda_4), 1 - 1/lambda_4, 1/(2 lambda_4) (they match its first five moments), the other
    sources at their mean. The output change h_k at those points gives the cumulants of h_k; for independent
    sources the cumulants of the sum add up. (Summing raw moments instead, as in the plain scheme, misses the
    cross terms of the fourth moment when there are thousands of sources and wrecks the Cornish-Fisher
    quantile.) Returns {target kind: (lo, hi)}, number of AC solves, solves not converged."""
    m = st.src.n
    e = np.sqrt(KURT_U * VAR_U)                   # 0.775, inside [-1, 1]
    p = 1.0 / (2.0 * KURT_U)
    kap = {k: np.zeros((W.shape[1], 4)) for k, W in Ws.items()}
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    n_bad = 0
    for s in range(0, m, batch):
        cols = np.arange(s, min(s + batch, m))
        X = np.zeros((m, 2 * cols.size))
        X[cols, np.arange(cols.size)] = e
        X[cols, cols.size + np.arange(cols.size)] = -e
        V, conv, _ = st.solver.solve(st.injections(X, delta))
        n_bad += int((~conv).sum())
        dv = np.abs(V[st.grid.pq]) - st.vm_star[:, None]
        dv[:, ~conv] = 0.0
        for k in Ws:
            y = WT[k] @ dv
            a, b = y[:, :cols.size], y[:, cols.size:]
            m1 = p * (a + b)
            m2 = p * (a ** 2 + b ** 2)
            m3 = p * (a ** 3 + b ** 3)
            m4 = p * (a ** 4 + b ** 4)
            kap[k][:, 0] += m1.sum(axis=1)
            kap[k][:, 1] += (m2 - m1 ** 2).sum(axis=1)
            kap[k][:, 2] += (m3 - 3 * m2 * m1 + 2 * m1 ** 3).sum(axis=1)
            kap[k][:, 3] += (m4 - 4 * m3 * m1 - 3 * m2 ** 2 + 12 * m2 * m1 ** 2 - 6 * m1 ** 4).sum(axis=1)
    out = {}
    for k, K in kap.items():
        var = K[:, 1]
        ok = var > 0
        g1 = np.divide(K[:, 2], var ** 1.5, out=np.zeros_like(var), where=ok)
        g2 = np.divide(K[:, 3], var ** 2, out=np.zeros_like(var), where=ok)
        out[k] = interval_from_moments(K[:, 0], var, g1, g2)
    return out, 2 * m, n_bad


def run(task):
    case, lf, tag, delta, n = task
    t0 = perf_counter()
    seed = seed_for("b8", case, lf, delta)
    st = study(case, load_factor=lf, threads=worker_threads())
    res = st.run(delta, method())
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    Ct = {k: np.asarray(W.T @ (st.C1 * delta)) for k, W in Ws.items()}
    t1 = perf_counter()
    so = SecondOrder(st, delta)
    rq = quadratic_ranges(so, Ws, N_SUR, seed + 1)
    t_ours = perf_counter() - t1
    bounds = {}
    for k, c in Ct.items():
        rb, s3 = bernstein_radius(c, ALPHA), 3.0 * sigma(c)
        bounds[(k, "bern")] = (-rb + np.minimum(0.0, rq[k][0]), rb + np.maximum(0.0, rq[k][1]))
        bounds[(k, "sig3")] = (-s3 + np.minimum(0.0, rq[k][0]), s3 + np.maximum(0.0, rq[k][1]))
        bounds[(k, "cumulant_cf")] = cumulant_cf(c)
    t1 = perf_counter()
    pe, n_pem, n_pem_bad = pem(st, delta, Ws)
    t_pem = perf_counter() - t1
    for k in Ws:
        bounds[(k, "pem_cf")] = pe[k]
    counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    t1 = perf_counter()
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ ac for k in Ws}
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    t_mc = perf_counter() - t1
    rows = []
    for (k, name), c in counts.items():
        lo, hi = bounds[(k, name)]
        w = hi - lo
        live = affine_radius(Ct[k]) > 1e-10
        wb = bounds[(k, "bern")][1] - bounds[(k, "bern")][0]
        rate = c / max(n_ac, 1)
        rows.append({"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "n_samples": n, "n_ac_not_converged": n - n_ac,
                     "seed": seed, "target": k, "method": name, "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()),
                     "n_over_alpha": int((rate > ALPHA).sum()), "any_target_pct": pct(anyc[(k, name)] / max(n_ac, 1)),
                     "width_over_bern_median": float(np.median(w[live] / wb[live])),
                     "width_over_affine_median": float(np.median(w[live] / (2 * affine_radius(Ct[k])[live]))),
                     "n_ac_solves": n_pem if name == "pem_cf" else 0, "pem_not_converged": n_pem_bad if name == "pem_cf" else 0,
                     "t_method_s": round(t_pem if name == "pem_cf" else (t_ours if name in ("bern", "sig3") else 0.0), 3),
                     "t_mc_s": round(t_mc, 2)})
    print(f"[b8] {case:26s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s (PEM {t_pem:.1f}s, MC {t_mc:.0f}s)", flush=True)
    return rows


def main():
    tasks = []
    for case, d in CASES.items():
        tasks.append((case, 1.0, "design", d, EXP["n_ac"]))
        tasks.append((case, 1.2, "stress_delta", 0.30, EXP["n_ac"]))
        lam, _ = loadability_limit(grid(case), tol=CFG["powerflow"]["tol_pu"])
        tasks.append((case, round(EXP["near_limit_fraction"] * lam, 4), "near_limit", 0.15, EXP["n_ac"]))
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b8_ppf_baselines", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR,
                                                "note": "cumulant_cf: linear model, Cornish-Fisher; pem_cf: Hong 2m+1 on AC, Cornish-Fisher; "
                                                        "t_method_s for pem_cf is its 2m AC solutions; for bern/sig3 it is the second-order model and the range"}))


if __name__ == "__main__":
    main()
