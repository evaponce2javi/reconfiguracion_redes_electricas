"""C5 -- System-level violation and the Bonferroni correction.

A 99.73 % set per branch does not give 99.73 % for the whole system. This script reports, against the
linear model and against nonlinear AC, the fraction of samples in which at least one bus (or one
ellipse) is violated and the expected number of violations per sample; then repeats with the chi2
quantile at level 1 - alpha/m (m = number of Ellipse-mode branches).

Writes results/c5_system_level.csv
"""
import numpy as np

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import monte_carlo, save_table, standard_bounds


def run(task):
    key, variant = task
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], CFG["experiments"]["n_ac"]
    seed = seed_for("ac_design_point", key)                       # same AC samples as B1
    st = study(sysc["case"], threads=worker_threads())
    base = st.run(delta, method())
    m = base.n_ellipse if variant == "bonferroni" else None
    rows = []
    for affine_gain in (False, True):                             # second pass only to get the area CRI vs the affine box
        kw = {"gain_box": "affine"} if affine_gain else {}
        res = st.run(delta, method(bonferroni_m=m, **kw))
        if affine_gain:
            rows[-1]["cri_area_vs_affine_box_pct"] = pct(np.mean(res.gain[base.ellipse]))
            continue
        res_c = st.run(delta, method(bonferroni_m=m, radius="coordinate"))
        cs = monte_carlo(st, delta, res, standard_bounds(res, res_c), n, seed, batch=CFG["powerflow"]["batch"])
        e = res.bus_ellipse
        row = {"system": key, "variant": variant, "delta": delta, "n_samples": n, "seed": seed,
               "n_ac_not_converged": cs.n - cs.n_ac, "m_ellipse_branches": res.n_ellipse, "n_pq_buses": len(e),
               "chi2_quantile": res.chi2, "radius_factor_vs_uncorrected": float(np.sqrt(res.chi2 / base.chi2)),
               "cri_area_cov_box_pct": pct(res.cri_area()),
               "width_red_lmax_vs_affine_pct": pct(res.width_reduction("affine")),
               "width_red_coord_vs_affine_pct": pct(res_c.width_reduction("affine")),
               "n_buses_lmax_wider_than_affine": int((res.r_out[e] > res.r_affine[e]).sum()),
               "n_buses_coord_wider_than_affine": int((res_c.r_out[e] > res.r_affine[e]).sum())}
        for ref in ("lin", "ac"):
            for bound, short in (("cdhim_lmax", "lmax"), ("cdhim_coord", "coord"), ("cltm_3sigma", "cltm"), ("affine_box", "affine")):
                s = cs.bus_summary(bound, ref)
                row[f"any_bus_{short}_{ref}_pct"] = pct(s["system"])
                row[f"expected_bus_violations_{short}_{ref}"] = s["expected_per_sample"]
            b = cs.branch_summary(ref)
            row[f"any_ellipse_{ref}_pct"] = pct(b["system"])
            row[f"expected_ellipse_violations_{ref}"] = b["expected_per_sample"]
            row[f"ellipse_{ref}_mean_pct"] = pct(b["mean"])
        rows.append(row)
    print(f"[c5] {key} {variant} done", flush=True)
    return rows


def main():
    rows = [r for out in parallel(run, [(k, v) for k in CFG["systems"] for v in ("per_branch_99.73", "bonferroni")], workers=4) for r in out]
    print(save_table("c5_system_level", rows, {
        "n_samples": CFG["experiments"]["n_ac"],
        "note": "any_* = share of samples with at least one violation; Bonferroni uses chi2_2 at 1 - 0.0027/m; "
                "cri_area_cov_box is scale free (unchanged by the quantile), cri_area_vs_affine_box is not"}))


if __name__ == "__main__":
    main()
