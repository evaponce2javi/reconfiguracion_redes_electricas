"""C3 -- Non-Gaussian input uncertainty, redone.

Every marginal is zero-mean. For each one the covariance is built twice:
  true_variance  Q = Var(eps) C C^T      (second moment correctly specified)
  uniform_1/3    Q = (1/3)   C C^T      (the manuscript's assumption, whatever the real marginal is)
Violation is reported for ALL Ellipse-mode branches (linear model, N = n_linear) and, with the
true-variance covariance, against nonlinear AC (N = n_c3_ac). Marginals whose support leaves [-1, 1]
are outside the affine model and are flagged as such. One skewed marginal is also combined with
correlated sources (Gaussian copula).

Writes results/c3_marginals.csv, results/c3_nongaussian.csv and results/c3_skewed_correlated.csv
"""
import numpy as np

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SourceCovariance, marginal_catalogue, monte_carlo, save_table, standard_bounds
from cdhim.uncertainty import UNIFORM_VARIANCE, copula_to_pearson

EXP = CFG["experiments"]
COMBO = ("beta25_support", 0.6)


def summarise(row, cs, res, refs, scale=None, suffix=""):
    for ref in refs:
        b = cs.branch_summary(ref, scale)
        row.update({f"ellipse_{ref}_mean_pct": pct(b["mean"]), f"ellipse_{ref}_p99_pct": pct(b["p99"]),
                    f"ellipse_{ref}_max_pct": pct(b["max"]), f"ellipse_{ref}_n_over": b["n_over"],
                    f"ellipse_{ref}_system_pct": pct(b["system"])})
        for bound, short in (("cdhim_lmax", "lmax"), ("cdhim_coord", "coord")):
            s = cs.bus_summary(bound + suffix, ref)
            row.update({f"bus_{short}_{ref}_mean_pct": pct(s["mean"]), f"bus_{short}_{ref}_max_pct": pct(s["max"]),
                        f"bus_{short}_{ref}_n_over": s["n_over"]})
    row["n_ellipse_branches"] = res.n_ellipse
    return row


def run(task):
    key, name = task
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    m = marginal_catalogue()[name]
    st = study(sysc["case"], threads=worker_threads())
    seed = seed_for("c3", key, name)
    # i.i.d. symbols: Q = Var(eps) C C^T. rho and the branch selection do not depend on Var(eps), so the
    # uniform assumption is the same set of ellipses scaled by k2 = (1/3)/Var(eps): one pass serves both.
    k2 = UNIFORM_VARIANCE / m.variance
    res = st.run(delta, method(), SourceCovariance(var=m.variance))
    res_c = st.run(delta, method(radius="coordinate"), SourceCovariance(var=m.variance))
    res_u = st.run(delta, method())
    res_cu = st.run(delta, method(radius="coordinate"))
    assert np.array_equal(res.ellipse, res_u.ellipse)
    bounds = {"cdhim_lmax": res.r_out, "cdhim_coord": res_c.r_out, "cdhim_lmax@u": res_u.r_out, "cdhim_coord@u": res_cu.r_out}
    base = {"system": key, "marginal": name, "in_model_support": m.in_model, "variance_eps": m.variance,
            "fraction_outside_[-1,1]_pct": pct(m.fraction_outside_unit), "delta": delta, "seed": seed}
    kw = dict(marginal=m, batch=CFG["powerflow"]["batch"], ellipse_scales={"u": k2})
    cs = monte_carlo(st, delta, res, bounds, EXP["n_linear"], seed, ac=False, **kw)
    cs_ac = monte_carlo(st, delta, res, bounds, EXP["n_nongaussian_ac"], seed + 1, ac=True, **kw)
    rows = []
    for label, scale, suffix in (("true_variance", None, ""), ("uniform_1/3", "u", "@u")):
        row = summarise({**base, "covariance_model": label, "n_linear": cs.n}, cs, res, ("lin",), scale, suffix)
        row.update({"n_ac": cs_ac.n, "n_ac_not_converged": cs_ac.n - cs_ac.n_ac})
        row.update({k: v for k, v in summarise({}, cs_ac, res, ("ac",), scale, suffix).items() if "_ac_" in k})
        rows.append(row)
    print(f"[c3] {key} {name} done", flush=True)
    return rows


def run_combo(key):
    """Skewed marginal + correlated sources: true covariance vs. the uniform, independent assumption."""
    name, r = COMBO
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    m = marginal_catalogue()[name]
    st = study(sysc["case"], threads=worker_threads())
    seed = seed_for("c3_combo", key)
    rho = copula_to_pearson(m, r)
    rows = []
    for label, sigma in (("true (variance and correlation)", SourceCovariance(var=m.variance, rho_c=rho)),
                         ("true variance, independent", SourceCovariance(var=m.variance)),
                         ("uniform 1/3, independent", None)):
        res = st.run(delta, method(), sigma)
        res_c = st.run(delta, method(radius="coordinate"), sigma)
        cs = monte_carlo(st, delta, res, standard_bounds(res, res_c), EXP["n_corr"], seed, marginal=m, copula_r=r,
                         batch=CFG["powerflow"]["batch"])
        row = {"system": key, "marginal": name, "copula_r": r, "pearson_rho_eps": rho, "covariance_model": label,
               "delta": delta, "seed": seed, "n_samples": cs.n, "n_ac_not_converged": cs.n - cs.n_ac,
               "cri_area_pct": pct(res.cri_area()), "width_red_lmax_vs_affine_pct": pct(res.width_reduction("affine"))}
        rows.append(summarise(row, cs, res, ("lin", "ac")))
    print(f"[c3] {key} combo done", flush=True)
    return rows


def main():
    cat = marginal_catalogue()
    marg = [{"marginal": n, "support_lo": m.support[0], "support_hi": m.support[1], "mean": m.mean, "variance": m.variance,
             "inside_[-1,1]": m.in_model, "fraction_outside_[-1,1]_pct": pct(m.fraction_outside_unit), "definition": m.note}
            for n, m in cat.items()]
    print(save_table("c3_marginals", marg, {"note": "exact moments of each marginal; checked by tests/test_uncertainty.py"}))
    rows = [r for out in parallel(run, [(k, n) for k in CFG["systems"] for n in cat], workers=6) for r in out]
    print(save_table("c3_nongaussian", rows, {"n_linear": EXP["n_linear"], "n_ac": EXP["n_nongaussian_ac"],
                     "note": "all Ellipse-mode branches; independent sources; AC samples use seed + 1"}))
    combo = [r for out in parallel(run_combo, list(CFG["systems"]), workers=2) for r in out]
    print(save_table("c3_skewed_correlated", combo, {"n_samples": EXP["n_corr"], "marginal": COMBO[0], "copula_r": COMBO[1]}))


if __name__ == "__main__":
    main()
