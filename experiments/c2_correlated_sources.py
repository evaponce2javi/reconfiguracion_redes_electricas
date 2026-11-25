"""C2 -- Correlated input uncertainty: Q = C Sigma_eps C^T against the mis-specified Q = (1/3) C C^T.

Sources are coupled by a one-factor Gaussian copula with parameter r and keep U[-1,1] marginals, so the
true covariance is Sigma_eps = (1/3)[(1 - rho_u) I + rho_u 11^T] with rho_u = (6/pi) asin(r/2).
Both covariance models are tested on the same samples (linear model and nonlinear AC).

Writes results/c2_correlated_sources.csv
"""
from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SourceCovariance, monte_carlo, save_table, standard_bounds
from cdhim.uncertainty import UNIFORM_VARIANCE, copula_to_pearson_uniform

EXP = CFG["experiments"]


def run(task):
    key, r = task
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    seed = seed_for("c2", key, r)
    st = study(sysc["case"], threads=worker_threads())
    rho_u = float(copula_to_pearson_uniform(r))
    models = {"misspecified_I/3": None, "true_Sigma_eps": SourceCovariance(var=UNIFORM_VARIANCE, rho_c=rho_u)}
    rows = []
    for name, sigma in models.items():
        res = st.run(delta, method(), sigma)
        res_c = st.run(delta, method(radius="coordinate"), sigma)
        cs = monte_carlo(st, delta, res, standard_bounds(res, res_c), EXP["n_corr"], seed, copula_r=r, batch=CFG["powerflow"]["batch"])
        row = {"system": key, "copula_r": r, "pearson_rho_eps": rho_u, "covariance_model": name, "delta": delta,
               "n_samples": cs.n, "seed": seed, "n_ac_not_converged": cs.n - cs.n_ac,
               "ellipse_branches_pct": pct(res.n_ellipse / res.n_branches), "n_ellipse_branches": res.n_ellipse,
               "cri_area_pct": pct(res.cri_area()), "width_red_lmax_vs_affine_pct": pct(res.width_reduction("affine")),
               "width_red_coord_vs_affine_pct": pct(res_c.width_reduction("affine")),
               "n_buses_lmax_wider_than_affine": int((res.r_out[res.bus_ellipse] > res.r_affine[res.bus_ellipse]).sum())}
        for ref in ("lin", "ac"):
            b = cs.branch_summary(ref)
            row.update({f"ellipse_{ref}_mean_pct": pct(b["mean"]), f"ellipse_{ref}_max_pct": pct(b["max"]),
                        f"ellipse_{ref}_n_over": b["n_over"], f"ellipse_{ref}_system_pct": pct(b["system"])})
            for bound, short in (("cdhim_lmax", "lmax"), ("cdhim_coord", "coord"), ("cltm_3sigma", "cltm"), ("affine_box", "affine")):
                s = cs.bus_summary(bound, ref)
                row.update({f"bus_{short}_{ref}_mean_pct": pct(s["mean"]), f"bus_{short}_{ref}_max_pct": pct(s["max"]),
                            f"bus_{short}_{ref}_n_over": s["n_over"], f"bus_{short}_{ref}_system_pct": pct(s["system"])})
        rows.append(row)
    print(f"[c2] {key} r={r} done", flush=True)
    return rows


def main():
    rows = [r for out in parallel(run, [(k, r) for k in CFG["systems"] for r in EXP["copula_r"]], workers=4) for r in out]
    print(save_table("c2_correlated_sources", rows, {
        "n_samples": EXP["n_corr"], "copula": "one-factor Gaussian, all P and Q symbols equicorrelated, U[-1,1] marginals",
        "note": "true_Sigma_eps uses Sigma_eps = (1/3)[(1-rho)I + rho 11^T]; r = 0 reproduces the default exactly"}))


if __name__ == "__main__":
    main()
