"""B12 -- A cheaper, narrower variant without guarantee: Gaussian radius plus the range of the second-order term.

Variant: r = min(z * sigma, worst case) with z = sqrt(chi2_2(1 - alpha)) = 3.44 (the normal quantile at 99.94 %, i.e. a
3 sigma radius with a safety factor; it is what made the hybrid interval of the submitted manuscript cover once
corrected, B4), plus the sampled range of the second-order term with n = 1000 samples (no certificate is claimed, so
the 2/(n+1) term is not paid). Compared on the same AC samples with:
  bern_cert   the certified bound (Bernstein at alpha - 2/(n+1), n = 4000)
  sig3        3 sigma plus the range (n = 1000), which fails on the PEGASE systems (B8)
Targets: buses, branch drops, branch-end currents. Cells: those of B8, B9 and B11.

    python experiments/b12_gaussian_variant.py
    python experiments/b12_gaussian_variant.py timing
Writes results/b12_gaussian_variant.csv (and results/b12_timing.csv)
"""
import sys
from time import perf_counter

import numpy as np
from scipy.stats import chi2

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, loadability_limit, sample_eps, save_table
from cdhim.affine_forms import (BranchCurrents, NonlinearSampler, affine_radius, bus_targets, certified_bernstein,
                                drop_targets, sigma)
from cdhim.validation import VIOLATION_ATOL
from b8_ppf_baselines import CASES

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
Z_GAUSS = float(np.sqrt(chi2.ppf(1 - ALPHA, 2)))       # 3.44
N_VARIANT, N_CERT = 1000, 4000


def setup(st, delta):
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    bc = BranchCurrents(st, so)
    C = {k: np.asarray(W.T @ (st.C1 * delta)) for k, W in Ws.items()}
    C["current"] = bc.C
    return so, Ws, bc, C


def run(task):
    case, lf, tag, delta, n = task
    t0 = perf_counter()
    seed = seed_for("b12", case, lf, delta)
    st = study(case, load_factor=lf, threads=worker_threads())
    so, Ws, bc, C = setup(st, delta)
    ns = NonlinearSampler(so, Ws, bc)
    r_var = ns.ranges(N_VARIANT, seed + 1)
    r_cert = ns.ranges(N_CERT, seed + 2)
    bounds = {}
    for k, c in C.items():
        aff, sg = affine_radius(c), sigma(c)
        for name, r, rq in (("gauss344", np.minimum(Z_GAUSS * sg, aff), r_var[k]),
                            ("sig3", np.minimum(3.0 * sg, aff), r_var[k]),
                            ("bern_cert", certified_bernstein(c, ALPHA, N_CERT), r_cert[k])):
            bounds[(k, name)] = (-r + np.minimum(0.0, rq[0]), r + np.maximum(0.0, rq[1]))
    counts = {b: np.zeros(C[b[0]].shape[0], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        dv = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ dv for k in Ws}
        tv["current"] = np.abs(bc.Yb @ V[:, conv]) - np.abs(bc.I0)[:, None]
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    rows = []
    for (k, name), c in counts.items():
        lo, hi = bounds[(k, name)]
        lo0, hi0 = bounds[(k, "bern_cert")]
        live = affine_radius(C[k]) > 1e-10
        rate = c / max(n_ac, 1)
        rows.append({"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "n_samples": n, "n_ac_not_converged": n - n_ac,
                     "seed": seed, "target": k, "bound": name, "z": Z_GAUSS if name == "gauss344" else (3.0 if name == "sig3" else np.nan),
                     "n_surrogate": N_CERT if name == "bern_cert" else N_VARIANT,
                     "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "any_target_pct": pct(anyc[(k, name)] / max(n_ac, 1)),
                     "width_over_certified_median": float(np.median((hi - lo)[live] / (hi0 - lo0)[live])),
                     "width_over_affine_median": float(np.median((hi - lo)[live] / (2 * affine_radius(C[k])[live])))})
    print(f"[b12] {case:26s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s", flush=True)
    return rows


def run_timing(case, delta, repeats=3):
    rows = []
    for rep in range(repeats):
        st = study(case, threads=1)
        T = {"linear_study_s": st.timings["jacobian"] + st.timings["lu"] + st.timings["solves"]}
        t = perf_counter(); so, Ws, bc, C = setup(st, delta); T["forms_and_second_order_s"] = perf_counter() - t
        for n in (N_VARIANT, N_CERT):
            t = perf_counter(); NonlinearSampler(so, Ws).ranges(n, rep); T[f"samples_{n}_voltages_s"] = perf_counter() - t
            t = perf_counter(); NonlinearSampler(so, Ws, bc).ranges(n, rep); T[f"samples_{n}_all_s"] = perf_counter() - t
        base = T["linear_study_s"] + T["forms_and_second_order_s"]
        rows.append({"case": case, "delta": delta, "repeat": rep, **{k: round(x, 3) for k, x in T.items()},
                     "total_variant_voltages_s": round(base + T[f"samples_{N_VARIANT}_voltages_s"], 3),
                     "total_variant_all_s": round(base + T[f"samples_{N_VARIANT}_all_s"], 3),
                     "total_certified_voltages_s": round(base + T[f"samples_{N_CERT}_voltages_s"], 3),
                     "total_certified_all_s": round(base + T[f"samples_{N_CERT}_all_s"], 3)})
        print(f"[b12] timing {case} {rows[-1]}", flush=True)
    return rows


def main():
    if "timing" in sys.argv[1:]:
        rows = [r for c, d in CASES.items() for r in run_timing(c, d)]
        print(save_table("b12_timing", rows, {"note": "one process, run alone; forms_and_second_order_s includes the current forms"}))
        return
    tasks = []
    for case, d in CASES.items():
        tasks.append((case, 1.0, "design", d, EXP["n_ac"]))
        tasks.append((case, 1.2, "stress_delta", 0.30, EXP["n_ac"]))
        lam, _ = loadability_limit(grid(case), tol=CFG["powerflow"]["tol_pu"])
        tasks.append((case, round(EXP["near_limit_fraction"] * lam, 4), "near_limit", 0.15, EXP["n_ac"]))
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b12_gaussian_variant", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "z": Z_GAUSS,
                                                    "note": "gauss344: sqrt(chi2_2(1-alpha)) sigma plus the range of 1000 samples, no "
                                                            "guarantee; bern_cert: certified bound with 4000 samples; sig3: 3 sigma "
                                                            "plus the range of 1000 samples"}))


if __name__ == "__main__":
    main()
