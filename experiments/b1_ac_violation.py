"""B1 -- Violation rate of the CD-HIM bounds against nonlinear AC power-flow samples (design point).

Three linearisations are evaluated with the same AC samples:
  legacy     flat start + H sign error (what produced the manuscript's figures)
  flat_exact flat start, exact Jacobian (sign fixed only)
  corrected  exact Jacobian at the converged AC solution (what the manuscript describes)

Per bus: AC sample outside the output interval. Per branch: AC sample outside the covariance ellipse.
Writes results/b1_bus_violation.csv and results/b1_branch_violation.csv
"""
from time import perf_counter

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import monte_carlo, save_table, standard_bounds

MODELS = {"legacy": dict(legacy=True), "flat_exact": dict(point="flat", jacobian="exact"), "corrected": {}}


def run(task):
    key, model = task
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], CFG["experiments"]["n_ac"]
    seed = seed_for("ac_design_point", key)                      # same AC samples for every model
    legacy = model == "legacy"
    st = study(sysc["case"], threads=worker_threads(), **MODELS[model])
    res = st.run(delta, method(legacy=legacy))
    res_c = st.run(delta, method(legacy=legacy, radius="coordinate"))
    t0 = perf_counter()
    cs = monte_carlo(st, delta, res, standard_bounds(res, res_c), n, seed, batch=CFG["powerflow"]["batch"])
    t_mc = perf_counter() - t0
    base = {"system": key, "model": model, "delta": delta, "n_samples": n, "seed": seed, "n_ac_converged": cs.n_ac,
            "n_ac_not_converged": cs.n - cs.n_ac, "newton_fallback": cs.solver.newton_fallback}
    bus_rows, br_rows = [], []
    sets = {"all_pq": None, "ellipse_buses": res.bus_ellipse, "box_buses": ~res.bus_ellipse}
    for bound in ("cdhim_lmax", "cdhim_coord", "affine_box", "cltm_3sigma"):
        for set_name, mask in sets.items():
            for ref in ("ac", "lin"):
                s = cs.bus_summary(bound, ref, mask)
                bus_rows.append({**base, "bound": bound, "bus_set": set_name, "reference": ref, "n_buses": s["n"],
                                 "mean_pct": pct(s["mean"]), "p99_pct": pct(s["p99"]), "max_pct": pct(s["max"]),
                                 "n_buses_over_0.27pct": s["n_over"],
                                 "system_level_pct": pct(s["system"]) if mask is None else "",
                                 "expected_violations_per_sample": s["expected_per_sample"] if mask is None else ""})
    for ref in ("ac", "lin"):
        s = cs.branch_summary(ref)
        br_rows.append({**base, "reference": ref, "n_ellipse_branches": s["n"], "mean_pct": pct(s["mean"]),
                        "p99_pct": pct(s["p99"]), "max_pct": pct(s["max"]), "n_branches_over_0.27pct": s["n_over"],
                        "system_level_pct": pct(s["system"]), "expected_violations_per_sample": s["expected_per_sample"]})
    print(f"[b1] {key:11s} {model:10s} done in {t_mc:.0f}s", flush=True)
    return bus_rows, br_rows


def main():
    out = parallel(run, [(k, m) for k in CFG["systems"] for m in MODELS], workers=3)
    meta = {"n_samples": CFG["experiments"]["n_ac"], "violation_atol_pu": 1e-8,
            "note": "non-converged AC samples are excluded from the denominator and counted in n_ac_not_converged; "
                    "bus radius = smallest radius assigned by Algorithm 1 (bus_combine=min), centre = converged AC solution"}
    print(save_table("b1_bus_violation", [r for b, _ in out for r in b], meta))
    print(save_table("b1_branch_violation", [r for _, b in out for r in b], meta))


if __name__ == "__main__":
    main()
