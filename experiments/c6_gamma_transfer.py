"""C6 -- Is gamma_th transferable between systems?

The threshold sweep of the manuscript (5-25 %) is repeated on five systems and for the three boxes that
can sit inside G. For each (system, box, gamma): share of Ellipse branches, CRI, the manuscript's
selection criterion (CRI over normalised run time) and the violation of the output intervals.

Writes results/c6_gamma_transfer.csv and results/c6_optimum.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import geometry as geo
from cdhim import monte_carlo, save_table
from cdhim.core import run_cdhim

EXP = CFG["experiments"]


def run(case):
    delta = CFG["transfer_delta"].get(case, CFG["transfer_delta"]["default"])
    st = study(case, threads=worker_threads())
    C = st.coefficients(delta)
    seed = seed_for("c6", case)
    results, bounds, times = {}, {}, {}
    for box in geo.GAIN_BOXES:
        for gamma in EXP["gammas"]:
            cfg = method(gain_box=box, gamma_th=gamma)
            t = []
            for _ in range(EXP["timing_repeats"]):
                t0 = perf_counter()
                res = run_cdhim(C, st.fi, st.ti, cfg)
                t.append(perf_counter() - t0)
            results[(box, gamma)], times[(box, gamma)] = res, float(np.median(t))
            bounds[f"{box}|{gamma}"] = res.r_out
    ref = results[("covariance_1sigma", 0.15)]
    cs = monte_carlo(st, delta, ref, bounds, EXP["n_transfer"], seed, batch=CFG["powerflow"]["batch"])
    rows = []
    for (box, gamma), res in results.items():
        t_ref = times[(box, 0.15)]
        ac, lin = cs.bus_summary(f"{box}|{gamma}", "ac"), cs.bus_summary(f"{box}|{gamma}", "lin")
        cri = res.cri_area()
        rows.append({"case": case, "n_bus": st.grid.n_bus, "delta": delta, "gain_box": box, "gamma_th": gamma,
                     "n_branches": res.n_branches, "ellipse_pct": pct(res.n_ellipse / max(res.n_branches, 1)),
                     "cri_area_pct": pct(cri) if res.n_ellipse else np.nan,
                     "width_red_vs_affine_pct": pct(res.width_reduction("affine")) if res.n_ellipse else np.nan,
                     "width_red_vs_affine_all_buses_pct": pct(res.width_reduction("affine", "all")),
                     "time_branch_stage_ms": 1e3 * times[(box, gamma)], "time_normalised": times[(box, gamma)] / t_ref,
                     "criterion_cri_over_time": (cri / (times[(box, gamma)] / t_ref)) if res.n_ellipse else np.nan,
                     "bus_ac_mean_pct": pct(ac["mean"]), "bus_ac_max_pct": pct(ac["max"]), "bus_ac_n_over": ac["n_over"],
                     "bus_ac_system_pct": pct(ac["system"]), "bus_lin_mean_pct": pct(lin["mean"]), "bus_lin_max_pct": pct(lin["max"]),
                     "n_ac": cs.n, "n_ac_not_converged": cs.n - cs.n_ac, "seed": seed})
    print(f"[c6] {case} done", flush=True)
    return rows


def optimum(rows):
    out = []
    for case in dict.fromkeys(r["case"] for r in rows):
        for box in geo.GAIN_BOXES:
            sub = [r for r in rows if r["case"] == case and r["gain_box"] == box]
            sel = sorted({r["ellipse_pct"] for r in sub})
            best_t = max(sub, key=lambda r: -1 if np.isnan(r["criterion_cri_over_time"]) else r["criterion_cri_over_time"])
            best_c = max(sub, key=lambda r: -1 if np.isnan(r["cri_area_pct"]) else r["cri_area_pct"])
            out.append({"case": case, "gain_box": box, "distinct_selections_over_sweep": len(sel),
                        "ellipse_pct_min": sel[0], "ellipse_pct_max": sel[-1],
                        "gamma_opt_cri_over_time": best_t["gamma_th"], "gamma_opt_cri_only": best_c["gamma_th"],
                        "time_spread_over_sweep_pct": pct(max(r["time_normalised"] for r in sub) - min(r["time_normalised"] for r in sub)),
                        "bus_ac_max_pct_range": f"{min(r['bus_ac_max_pct'] for r in sub):.3f}-{max(r['bus_ac_max_pct'] for r in sub):.3f}"})
    return out


def main():
    rows = [r for out in parallel(run, CFG["transfer_cases"], workers=5) for r in out]
    meta = {"n_ac": EXP["n_transfer"], "timing_repeats": EXP["timing_repeats"],
            "note": "time = median of the branch stage (covariances, eigenvalues, decision, radii); the decision needs the "
                    "spectra of every branch, so run time cannot depend on gamma_th beyond noise"}
    print(save_table("c6_gamma_transfer", rows, meta))
    print(save_table("c6_optimum", optimum(rows), meta))


if __name__ == "__main__":
    main()
