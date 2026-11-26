"""B10 -- Non-uniform source distributions in the affine-forms method (design point).

The twelve marginals of C3. Each radius rests on different assumptions about the symbols:
  hoeffding        zero mean, support [-1, 1]                            (any such marginal)
  bern_var13       zero mean, support [-1, 1], variance <= 1/3           (the default radius)
  bern_true_var    zero mean, support [-1, 1], the marginal's variance   (known marginal)
  sig3_true_var    3 sigma with the true variance                        (normality)
The marginals outside the model (unbounded, or rescaled to leave [-1, 1]) show where the assumptions break.
All radii carry the range of the quadratic term sampled under the true marginal (1.000 samples).

    python experiments/b10_marginals.py
Writes results/b10_marginals.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, marginal_catalogue, sample_eps, save_table
from cdhim.affine_forms import affine_radius, bernstein_radius, bus_targets, drop_targets, quadratic_ranges, sigma
from cdhim.validation import VIOLATION_ATOL

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
N_SUR = 1000
CASES = {"case33bw": 0.15, "case118": 0.05, "pglib_opf_case2869_pegase": 0.05}


def run(task):
    case, delta, name = task
    t0 = perf_counter()
    marg = marginal_catalogue()[name]
    seed = seed_for("b10", case, name)
    st = study(case, threads=worker_threads())
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    sampler = lambda rng, m: sample_eps(rng, st.src.n, m, marg)
    rq = quadratic_ranges(so, Ws, N_SUR, seed + 1, sampler=sampler)
    var = marg.variance
    bounds = {}
    for k, W in Ws.items():
        Ct = np.asarray(W.T @ (st.C1 * delta))
        L = np.log(2.0 / ALPHA)
        radii = {"hoeffding": np.minimum(np.sqrt(2.0 * L * (Ct ** 2).sum(axis=1)), affine_radius(Ct)),
                 "bern_var13": bernstein_radius(Ct, ALPHA),
                 "bern_true_var": bernstein_radius(Ct, ALPHA, var_eps=var),
                 "sig3_true_var": 3.0 * np.sqrt(var) * np.sqrt((Ct ** 2).sum(axis=1)),
                 "affine": affine_radius(Ct)}
        for r, rr in radii.items():
            bounds[(k, r)] = (-rr + np.minimum(0.0, rq[k][0]), rr + np.maximum(0.0, rq[k][1]))
        bounds[(k, "_ref")] = (-affine_radius(Ct), affine_radius(Ct))
    counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac, n = 0, EXP["n_ac"]
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sampler(rng, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ ac for k in Ws}
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    rows = []
    for (k, r), c in counts.items():
        if r == "_ref":
            continue
        lo, hi = bounds[(k, r)]
        a = bounds[(k, "_ref")][1]
        live = a > 1e-10
        rate = c / max(n_ac, 1)
        lo_s, hi_s = marg.support
        rows.append({"case": case, "delta": delta, "marginal": name, "in_model": marg.in_model, "eps_mean": marg.mean,
                     "eps_variance": var, "support_lo": lo_s, "support_hi": hi_s, "n_samples": n, "n_ac_not_converged": n - n_ac,
                     "seed": seed, "target": k, "radius": r, "max_pct": pct(rate.max()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "any_target_pct": pct(anyc[(k, r)] / max(n_ac, 1)),
                     "width_over_affine_median": float(np.median((0.5 * (hi - lo))[live] / a[live]))})
    print(f"[b10] {case:26s} {name:18s} {perf_counter() - t0:4.0f}s", flush=True)
    return rows


def main():
    tasks = [(c, d, m) for c, d in CASES.items() for m in marginal_catalogue()]
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b10_marginals", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR,
                                             "note": "marginals of C3 (cdhim.marginal_catalogue); ranges of q sampled under the true marginal"}))


if __name__ == "__main__":
    main()
