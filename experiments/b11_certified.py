"""B11 -- Certified level: Bernstein radius at alpha - 2/(n+1) plus the range of n samples of the quadratic term.

By the union bound and order statistics (cdhim.affine_forms.certified_split), for the second-order model
    P(target outside) <= alpha_L + 2/(n + 1) = alpha.
This experiment measures what the certificate costs in width against the uncertified bound of B5 to B9
(Bernstein at alpha with the range of 1.000 samples), and checks both against AC Monte Carlo.
Targets: bus |V|, branch drops |V_i| - |V_j| and branch-end currents |I|. Cells: those of B8 and B9.

    python experiments/b11_certified.py
Writes results/b11_certified.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, loadability_limit, sample_eps, save_table
from cdhim.affine_forms import affine_radius, bernstein_radius, bus_targets, certified_split, drop_targets
from cdhim.validation import VIOLATION_ATOL
from b8_ppf_baselines import CASES
from b9_branch_currents import current_model

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
SIZES = (1000, 4000, 10000)
N_REF = 1000                                         # the uncertified bound uses the first 1.000 samples


def nonlinear_samples(st, so, Ws, Yb, I0, Ci, n, seed, batch=250):
    """n samples of the nonlinear part of every target: quadratic term for |V| targets, remainder of |I|."""
    rng = np.random.default_rng(seed)
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    out = {k: [] for k in Ws}
    out["current"] = []
    for s in range(0, n, batch):
        eps = sample_eps(rng, st.src.n, min(batch, n - s))
        q = so.quadratic(eps)
        for k in Ws:
            out[k].append(WT[k] @ q)
        out["current"].append(np.abs(I0[:, None] + Yb @ so.state(eps)) - np.abs(I0)[:, None] - Ci @ eps)
    return {k: np.hstack(v) for k, v in out.items()}


def run(task):
    case, lf, tag, delta, n_mc = task
    t0 = perf_counter()
    seed = seed_for("b11", case, lf, delta)
    st = study(case, load_factor=lf, threads=worker_threads())
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    Yb, I0, Ci, _ = current_model(st, so)
    C = {k: np.asarray(W.T @ (st.C1 * delta)) for k, W in Ws.items()}
    C["current"] = Ci
    t1 = perf_counter()
    Q = nonlinear_samples(st, so, Ws, Yb, I0, Ci, max(SIZES), seed + 1)
    t_q = perf_counter() - t1
    bounds = {}
    for k, c in C.items():
        rb = bernstein_radius(c, ALPHA)
        q = Q[k][:, :N_REF]
        bounds[(k, "uncertified", N_REF)] = (-rb + np.minimum(0, q.min(axis=1)), rb + np.maximum(0, q.max(axis=1)))
        for m in SIZES:
            a = certified_split(ALPHA, m)
            r = bernstein_radius(c, a)
            q = Q[k][:, :m]
            bounds[(k, "certified", m)] = (-r + np.minimum(0, q.min(axis=1)), r + np.maximum(0, q.max(axis=1)))
    counts = {b: np.zeros(C[b[0]].shape[0], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n_mc, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n_mc - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        dv = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ dv for k in Ws}
        tv["current"] = np.abs(Yb @ V[:, conv]) - np.abs(I0)[:, None]
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    rows = []
    for (k, kind, m), c in counts.items():
        lo, hi = bounds[(k, kind, m)]
        lo0, hi0 = bounds[(k, "uncertified", N_REF)]
        live = affine_radius(C[k]) > 1e-10
        rate = c / max(n_ac, 1)
        a_l = ALPHA if kind == "uncertified" else certified_split(ALPHA, m)
        rows.append({"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "n_samples": n_mc,
                     "n_ac_not_converged": n_mc - n_ac, "seed": seed, "target": k, "kind": kind, "n_surrogate": m,
                     "alpha_linear": a_l, "certified_level": a_l + 2.0 / (m + 1) if kind == "certified" else np.nan,
                     "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "any_target_pct": pct(anyc[(k, kind, m)] / max(n_ac, 1)),
                     "width_over_uncertified_median": float(np.median((hi - lo)[live] / (hi0 - lo0)[live])),
                     "width_over_uncertified_max": float(np.max((hi - lo)[live] / (hi0 - lo0)[live])),
                     "width_over_affine_median": float(np.median((hi - lo)[live] / (2 * affine_radius(C[k])[live]))),
                     "t_samples_s": round(t_q * m / max(SIZES), 2)})
    print(f"[b11] {case:26s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s (samples {t_q:.1f}s)", flush=True)
    return rows


def main():
    tasks = []
    for case, d in CASES.items():
        tasks.append((case, 1.0, "design", d, EXP["n_ac"]))
        tasks.append((case, 1.2, "stress_delta", 0.30, EXP["n_ac"]))
        lam, _ = loadability_limit(grid(case), tol=CFG["powerflow"]["tol_pu"])
        tasks.append((case, round(EXP["near_limit_fraction"] * lam, 4), "near_limit", 0.15, EXP["n_ac"]))
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b11_certified", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "sizes": list(SIZES),
                                             "note": "certified: Bernstein at alpha - 2/(n+1) plus the range of n samples of the "
                                                     "nonlinear part; uncertified: Bernstein at alpha plus the range of 1.000 samples; "
                                                     "t_samples_s prorated from the 10.000-sample draw"}))


if __name__ == "__main__":
    main()
