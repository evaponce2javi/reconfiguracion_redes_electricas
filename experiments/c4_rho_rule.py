"""C4 -- Gain-based switching rule against a direct correlation threshold.

For the two covariance boxes kappa is a constant, so the gain rule IS a correlation rule with
rho_th = sqrt(1 - ((1 - gamma)/kappa)^2); the script verifies that both rules pick the same branches.
For the affine box kappa changes from branch to branch; the closest correlation rule is found by
scanning rho_th and the branches that change are reported with their effect.

Writes results/c4_rho_rule.csv and results/c4_kappa.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, method, pct, seed_for, study
from cdhim import geometry as geo
from cdhim import monte_carlo, save_table
from cdhim.core import run_cdhim

N_MC = CFG["experiments"]["n_rule_mc"]


def timed(fn, repeats=20):
    t = []
    for _ in range(repeats):
        t0 = perf_counter()
        out = fn()
        t.append(perf_counter() - t0)
    return out, float(np.median(t))


def main():
    rows, kap = [], []
    gamma = CFG["method"]["gamma_th"]
    for key, sysc in CFG["systems"].items():
        delta = sysc["delta"]
        for model, legacy in (("legacy", True), ("corrected", False)):
            st = study(sysc["case"], legacy=legacy)
            C = st.coefficients(delta)
            for box in geo.GAIN_BOXES:
                res, t_gain = timed(lambda: run_cdhim(C, st.fi, st.ti, method(legacy=legacy, gain_box=box)))
                ok = res.valid & (np.abs(res.rho) < 1 - res.cfg.delta_stab)
                k = res.kappa[ok]
                const = geo.kappa_constant(box, res.chi2)
                kap.append({"system": key, "model": model, "gain_box": box, "kappa_min": float(np.nanmin(k)),
                            "kappa_median": float(np.nanmedian(k)), "kappa_max": float(np.nanmax(k)),
                            "kappa_closed_form": const if const is not None else "branch dependent",
                            "min_gain_pct": pct(res.gain[ok].min()), "gain_at_rho_0_pct": pct(1 - const) if const else ""})
                if const is not None:
                    rho_th = geo.rho_threshold(gamma, const)
                else:                                              # no exact equivalent: best matching threshold
                    cand = np.unique(np.concatenate([[0.0], np.abs(res.rho[ok])]))
                    diff = [(np.count_nonzero((np.abs(res.rho[ok]) >= c) != res.ellipse[ok]), c) for c in cand]
                    rho_th = float(min(diff)[1])
                res_r, t_rho = timed(lambda: run_cdhim(C, st.fi, st.ti, method(legacy=legacy, gain_box=box,
                                                                                selection_rule="rho", rho_th=rho_th)))
                changed = res.ellipse != res_r.ellipse
                row = {"system": key, "model": model, "gain_box": box, "gamma_th": gamma, "rho_th_equivalent": rho_th,
                       "exact_equivalence_expected": const is not None, "n_branches": res.n_branches,
                       "n_ellipse_gain_rule": res.n_ellipse, "n_ellipse_rho_rule": res_r.n_ellipse,
                       "n_branches_changed": int(changed.sum()), "same_selection": bool(not changed.any()),
                       "cri_area_gain_pct": pct(res.cri_area()), "cri_area_rho_pct": pct(res_r.cri_area()),
                       "width_red_vs_affine_gain_pct": pct(res.width_reduction("affine")),
                       "width_red_vs_affine_rho_pct": pct(res_r.width_reduction("affine")),
                       "time_gain_rule_ms": 1e3 * t_gain, "time_rho_rule_ms": 1e3 * t_rho}
                if changed.any() and not legacy:                   # effect on violations (linear and AC)
                    seed = seed_for("c4", key, box)
                    for tag, r in (("gain", res), ("rho", res_r)):
                        cs = monte_carlo(st, delta, r, {"cdhim_lmax": r.r_out}, N_MC, seed, batch=CFG["powerflow"]["batch"])
                        row.update({f"bus_ac_mean_{tag}_pct": pct(cs.bus_summary("cdhim_lmax", "ac")["mean"]),
                                    f"bus_ac_max_{tag}_pct": pct(cs.bus_summary("cdhim_lmax", "ac")["max"]),
                                    f"ellipse_lin_mean_{tag}_pct": pct(cs.branch_summary("lin")["mean"]),
                                    f"ellipse_ac_mean_{tag}_pct": pct(cs.branch_summary("ac")["mean"]),
                                    "mc_seed": seed, "mc_n": N_MC})
                rows.append(row)
                print(f"[c4] {key} {model} {box}: changed={int(changed.sum())}", flush=True)
    meta = {"gamma_th": gamma, "note": "times are medians of 20 runs of the branch stage only"}
    print(save_table("c4_kappa", kap, meta))
    print(save_table("c4_rho_rule", rows, meta))


if __name__ == "__main__":
    main()
