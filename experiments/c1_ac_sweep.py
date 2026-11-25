"""C1 -- Validation against nonlinear AC over loading level x uncertainty level.

Loading: demand and scheduled generation scaled by the load factor (slack covers losses); the fifth
level is ``near_limit_fraction`` of the loadability limit found by continuation. In every cell the
same noise draws give three separate numbers:
  (a) linearisation error |dV_AC - dV_lin| per bus (absolute, and relative to the output radius)
  (b) miscoverage of the ellipse inside the linear model
  (c) total violation: AC samples outside the CD-HIM bounds

Writes results/c1_ac_sweep.csv and results/c1_first_exceedance.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import loadability_limit, monte_carlo, save_table, standard_bounds

EXP = CFG["experiments"]


def run(task):
    key, lf, tag, lam_max = task
    sysc = CFG["systems"][key]
    seed = seed_for("c1", key)                                   # common random numbers across all cells
    rows = []
    try:
        st = study(sysc["case"], load_factor=lf, threads=worker_threads())
    except RuntimeError as exc:                                  # no AC solution at this loading
        for delta in EXP["deltas"]:
            rows.append({"system": key, "load_factor": lf, "load_level": tag, "loadability_limit": lam_max, "delta": delta,
                         "n_samples": 0, "seed": seed, "status": f"blocked: {exc}"})
        return rows
    for delta in EXP["deltas"]:
        res = st.run(delta, method())
        res_c = st.run(delta, method(radius="coordinate"))
        t0 = perf_counter()
        cs = monte_carlo(st, delta, res, standard_bounds(res, res_c), EXP["n_sweep"], seed, batch=CFG["powerflow"]["batch"])
        ok = res.r_out > method().min_radius_pu
        row = {"system": key, "load_factor": lf, "load_level": tag, "loadability_limit": lam_max, "delta": delta,
               "n_samples": cs.n, "seed": seed, "status": "ok", "min_vm_base_pu": float(np.abs(st.V_star).min()),
               "nonconvergence_pct": pct(1 - cs.n_ac / cs.n), "newton_fallback": cs.solver.newton_fallback,
               "n_ellipse_branches": res.n_ellipse, "n_branches": res.n_branches,
               # (a) linearisation error
               "a_err_max_pu": cs.err_abs.max, "a_err_p99_pu": cs.err_abs.quantile(0.99),
               "a_err_over_radius_max": float(np.max(cs.err_bus_max[ok] / res.r_out[ok])) if cs.n_ac else np.nan,
               "a_err_over_radius_p99": cs.err_rel.quantile(0.99), "a_median_radius_pu": float(np.median(res.r_out[ok]))}
        b = cs.branch_summary("lin")                              # (b) ellipse miscoverage, linear model
        row.update({"b_ellipse_lin_mean_pct": pct(b["mean"]), "b_ellipse_lin_p99_pct": pct(b["p99"]),
                    "b_ellipse_lin_max_pct": pct(b["max"]), "b_ellipse_lin_n_over": b["n_over"]})
        if cs.n_ac:                                               # (c) total violation against AC
            for bound, short in (("cdhim_lmax", "lmax"), ("cdhim_coord", "coord"), ("affine_box", "affine"), ("cltm_3sigma", "cltm")):
                s = cs.bus_summary(bound, "ac")
                row.update({f"c_bus_{short}_mean_pct": pct(s["mean"]), f"c_bus_{short}_p99_pct": pct(s["p99"]),
                            f"c_bus_{short}_max_pct": pct(s["max"]), f"c_bus_{short}_n_over": s["n_over"],
                            f"c_bus_{short}_system_pct": pct(s["system"])})
            s = cs.branch_summary("ac")
            row.update({"c_ellipse_ac_mean_pct": pct(s["mean"]), "c_ellipse_ac_p99_pct": pct(s["p99"]),
                        "c_ellipse_ac_max_pct": pct(s["max"]), "c_ellipse_ac_n_over": s["n_over"],
                        "c_ellipse_ac_system_pct": pct(s["system"])})
        rows.append(row)
        print(f"[c1] {key:11s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s nonconv={row['nonconvergence_pct']:.2f}%", flush=True)
    return rows


def first_exceedance(rows):
    """First (load factor, delta), scanning load then delta, where each criterion passes 0.27 %."""
    crit = {"bus mean (lambda_max interval)": "c_bus_lmax_mean_pct", "bus max (lambda_max interval)": "c_bus_lmax_max_pct",
            "bus max (coordinate interval)": "c_bus_coord_max_pct", "ellipse mean (Mahalanobis)": "c_ellipse_ac_mean_pct",
            "ellipse max (Mahalanobis)": "c_ellipse_ac_max_pct", "system level (lambda_max interval)": "c_bus_lmax_system_pct"}
    out = []
    for key in CFG["systems"]:
        cells = sorted((r for r in rows if r["system"] == key and r["status"] == "ok" and r["n_samples"] and r["nonconvergence_pct"] < 100),
                       key=lambda r: (r["load_factor"], r["delta"]))
        for name, col in crit.items():
            over = [r for r in cells if r.get(col, 0) > 0.27]
            first = over[0] if over else None
            out.append({"system": key, "criterion": name, "n_cells": len(cells), "n_cells_over_0.27pct": len(over),
                        "first_load_factor": first["load_factor"] if first else "", "first_delta": first["delta"] if first else "",
                        "value_pct": first[col] if first else ""})
    return out


def main():
    tasks = []
    for key, sysc in CFG["systems"].items():
        lam_max, _ = loadability_limit(grid(sysc["case"]), tol=CFG["powerflow"]["tol_pu"])
        print(f"[c1] {key}: loadability limit = {lam_max:.3f}", flush=True)
        levels = [(lf, "fixed") for lf in EXP["load_factors"]] + [(round(EXP["near_limit_fraction"] * lam_max, 4), "near_limit")]
        tasks += [(key, lf, tag, lam_max) for lf, tag in levels]
    rows = [r for out in parallel(run, tasks, workers=5) for r in out]
    meta = {"n_samples_per_cell": EXP["n_sweep"], "load_scaling": "P_d, Q_d and scheduled P_g multiplied by the load factor",
            "near_limit_fraction": EXP["near_limit_fraction"], "radius_reference_for_relative_error": "CD-HIM lambda_max output radius"}
    print(save_table("c1_ac_sweep", rows, meta))
    print(save_table("c1_first_exceedance", first_exceedance(rows), meta))


if __name__ == "__main__":
    main()
