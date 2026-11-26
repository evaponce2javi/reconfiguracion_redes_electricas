"""B7 -- Simultaneous coverage and two more systems for the affine-forms method.

(1) Simultaneous coverage. The per-target radii of B5 hold each bus (or drop) at level alpha; the probability that
SOME target is outside is much larger (C5). With a union bound each target gets alpha / m. A concentration radius
grows only like sqrt(ln(m / alpha)) with m, a Gaussian one like Phi^-1(1 - alpha / 2m), Chebyshev like sqrt(m).
Measured: the share of AC samples with at least one bus (or one drop, or one of either) outside, and the widths.

(2) Generality: case118 and case1354pegase (both with PV buses), over a load x delta grid, with the same bounds
as B5 (per target and simultaneous).

Every bound is centred on the converged AC solution and carries the sampled range of the quadratic term
(1.000 samples), as in B5. Independent U[-1, 1] symbols, slack balancing.

    python experiments/b7_generality.py
Writes results/b7_generality.csv
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
NEW = {"case118": "IEEE-118", "case1354pegase": "PEGASE-1354"}
GRID_LOADS = (0.8, 1.0, 1.2)            # plus 95 % of the loadability limit
GRID_DELTAS = (0.05, 0.15, 0.30)


def radii(Ct, a):
    sg, r_aff = sigma(Ct), affine_radius(Ct)
    return {"bernstein": bernstein_radius(Ct, a), "gauss": np.minimum(sg * ndtri(1.0 - a / 2.0), r_aff),
            "chebyshev": np.minimum(sg / np.sqrt(a), r_aff), "affine": r_aff}


def run_cell(task):
    case, lf, tag, delta, n = task
    t0 = perf_counter()
    seed = seed_for("b7", case, lf, delta)
    try:
        st = study(case, load_factor=lf, threads=worker_threads())
    except RuntimeError as exc:
        return [{"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "status": f"blocked: {exc}"}]
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    Ct = {k: np.asarray(W.T @ (st.C1 * delta)) for k, W in Ws.items()}
    live = {k: affine_radius(c) > 1e-10 for k, c in Ct.items()}          # targets with any uncertainty
    m = {k: int(x.sum()) for k, x in live.items()}
    m["both"] = m["bus"] + m["drop"]
    rq = quadratic_ranges(so, Ws, N_SUR, seed + 1)
    bounds = {}
    for k, c in Ct.items():
        for level, a in (("per_target", ALPHA), ("bonferroni_own", ALPHA / max(m[k], 1)), ("bonferroni_both", ALPHA / max(m["both"], 1))):
            for name, r in radii(c, a).items():
                if name == "affine" and level != "per_target":
                    continue
                bounds[(k, level, name)] = (-r + np.minimum(0.0, rq[k][0]), r + np.maximum(0.0, rq[k][1]))
    counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
    any_own = {b: 0 for b in bounds}
    any_both = {}                                                       # (level, name): any bus or drop outside
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ ac for k in Ws}
        bad_any = {}
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            ba = bad.any(axis=0)
            any_own[bk] += int(ba.sum())
            key = (bk[1], bk[2])
            bad_any[key] = ba if key not in bad_any else (bad_any[key] | ba)
        for key, ba in bad_any.items():
            any_both[key] = any_both.get(key, 0) + int(ba.sum())
    rows = []
    base = {"case": case, "load_factor": lf, "load_level": tag, "delta": delta, "status": "ok", "n_samples": n,
            "n_ac_not_converged": n - n_ac, "seed": seed, "n_surrogate": N_SUR, "n_bus_targets": m["bus"], "n_drop_targets": m["drop"]}
    ref = {k: radii(c, ALPHA) for k, c in Ct.items()}
    for (k, level, name), c in counts.items():
        lo, hi = bounds[(k, level, name)]
        half = 0.5 * (hi - lo)
        ok = live[k]
        rate = c / max(n_ac, 1)
        a = ALPHA if level == "per_target" else ALPHA / max(m[k if level == "bonferroni_own" else "both"], 1)
        capped = radii(Ct[k], a)[name][ok] >= ref[k]["affine"][ok] * (1 - 1e-12)
        rows.append({**base, "target": k, "level": level, "radius": name, "alpha_per_target": a,
                     "max_pct": pct(rate.max()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "any_target_pct": pct(any_own[(k, level, name)] / max(n_ac, 1)),
                     "any_bus_or_drop_pct": pct(any_both[(level, name)] / max(n_ac, 1)),
                     "width_over_affine_median": float(np.median(half[ok] / ref[k]["affine"][ok])),
                     "width_over_bern_per_target_median": float(np.median(half[ok] / ref[k]["bernstein"][ok])),
                     "capped_by_affine_pct": pct(np.mean(capped))})
    print(f"[b7] {case:15s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s not converged {n - n_ac}", flush=True)
    return rows


def main():
    tasks = []
    for key, sysc in CFG["systems"].items():                            # design points of the two original systems
        tasks.append((sysc["case"], 1.0, "design", sysc["delta"], EXP["n_ac"]))
    for case in NEW:
        lam, _ = loadability_limit(grid(case), tol=CFG["powerflow"]["tol_pu"])
        print(f"[b7] {case}: loadability limit {lam:.3f}", flush=True)
        levels = [(lf, "fixed") for lf in GRID_LOADS] + [(round(EXP["near_limit_fraction"] * lam, 4), "near_limit")]
        tasks += [(case, lf, tag, d, EXP["n_sweep"]) for lf, tag in levels for d in GRID_DELTAS]
    rows = [r for o in parallel(run_cell, tasks, workers=3) for r in o]
    print(save_table("b7_generality", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR,
                                             "note": "bonferroni_own: alpha / (number of targets of that kind); bonferroni_both: alpha / (buses + drops); "
                                                     "any_target_pct: samples with some target of that kind outside; any_bus_or_drop_pct: some bus or drop "
                                                     "outside, both kinds at the same level"}))


if __name__ == "__main__":
    main()
