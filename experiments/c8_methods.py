"""C8 -- Comparison of the interval methods available in the repository, on the same linear model.

  affine_box   deterministic affine-arithmetic enclosure, half-width sum|c_k|
  cltm_3sigma  3-sigma box of the linear combination (the "CLTM Box" of the manuscript's Fig. 4; this is
               the CLT step of Ref. [10] applied to Jacobian sensitivities, NOT a re-implementation of
               Ref. [10], whose AC variant obtains the affine forms with an optimizing-scenarios step)
  cdhim_lmax   CD-HIM with the lambda_max radius (manuscript Eq. 11)
  cdhim_coord  CD-HIM with the exact coordinate projection
  ac_mc        nonlinear AC Monte Carlo (reference): central 99.73 % range and full range of the samples

For each: interval width, violation against AC, and run time. Refs. [11] and [12] are not implemented.
Run it alone (timings). Writes results/c8_methods.csv and results/c8_timing.csv
"""
import warnings
from time import perf_counter

import numpy as np

from common import CFG, grid, method, pct, seed_for, study
from cdhim import geometry as geo
from cdhim import sample_eps, save_table
from cdhim.core import run_cdhim
from cdhim.validation import VIOLATION_ATOL

N_PP = 200          # pandapower power flows timed to estimate the per-sample cost of the original MCS loop


def ac_samples(st, delta, n, seed, batch):
    rng = np.random.default_rng(seed)
    out = np.empty((st.C1.shape[0], n), dtype=np.float64)
    ok = np.zeros(n, bool)
    t = 0.0
    for s in range(0, n, batch):
        eps = sample_eps(rng, st.src.n, min(batch, n - s))
        t0 = perf_counter()
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        t += perf_counter() - t0
        out[:, s:s + eps.shape[1]] = np.abs(V[st.grid.pq]) - st.vm_star[:, None]
        ok[s:s + eps.shape[1]] = conv
    return out[:, ok], t, int((~ok).sum())


def pandapower_loop_time(case, delta, n=N_PP):
    """Per-sample wall time of the original protocol: pandapower.runpp(numba=False) in a Python loop."""
    import pandapower as pp
    from cdhim.network import _build_net, DATA_DIR

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net, _ = _build_net(case, DATA_DIR)
        p0, q0 = net.load.p_mw.values.copy(), net.load.q_mvar.values.copy()
        rng = np.random.default_rng(0)
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, numba=False)
        t = []
        for _ in range(n):
            net.load["p_mw"] = p0 * (1 + rng.uniform(-delta, delta, len(p0)))
            net.load["q_mvar"] = q0 * (1 + rng.uniform(-delta, delta, len(q0)))
            t0 = perf_counter()
            pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, numba=False)
            t.append(perf_counter() - t0)
    return float(np.mean(t)), float(np.std(t, ddof=1))


def timed(fn, repeats):
    t = []
    for _ in range(repeats):
        t0 = perf_counter()
        out = fn()
        t.append(perf_counter() - t0)
    return out, float(np.mean(t)), float(np.std(t, ddof=1))


def main():
    rows, trows = [], []
    rep = CFG["experiments"]["timing_repeats"]
    for key, sysc in CFG["systems"].items():
        delta, n = sysc["delta"], CFG["experiments"]["n_ac"]
        seed = seed_for("ac_design_point", key)
        grid(sysc["case"])                                        # load the case before any timing
        t_stage = []
        for _ in range(rep):
            st = study(sysc["case"], threads=1)
            t_stage.append(st.timings["jacobian"] + st.timings["lu"] + st.timings["solves"])
        t_lin = float(np.mean(t_stage))
        C = st.coefficients(delta)
        res, t_cd, sd_cd = timed(lambda: run_cdhim(C, st.fi, st.ti, method()), rep)
        res_c = run_cdhim(C, st.fi, st.ti, method(radius="coordinate"))
        _, t_aff, _ = timed(lambda: geo.affine_box_halfwidth(C), rep)
        _, t_cltm, _ = timed(lambda: 3.0 * np.sqrt(np.einsum("ij,ij->i", C, C) / 3.0), rep)
        dv, t_mc, n_fail = ac_samples(st, delta, n, seed, CFG["powerflow"]["batch"])
        t_pp, sd_pp = pandapower_loop_time(sysc["case"], delta)
        lo, hi = np.quantile(dv, [0.00135, 0.99865], axis=1)
        r_mc = 0.5 * (hi - lo)                                   # half-width of the central 99.73 % range
        r_full = np.max(np.abs(dv), axis=1)
        ok = res.r_affine > method().min_radius_pu
        radii = {"affine_box": res.r_affine, "cltm_3sigma": res.r_cltm, "cdhim_lmax": res.r_out, "cdhim_coord": res_c.r_out,
                 "ac_mc_central_99.73": r_mc, "ac_mc_full_range": r_full}
        guarantee = {"affine_box": "deterministic (linear model)", "cltm_3sigma": "99.73 % per bus (marginal)",
                     "cdhim_lmax": "hybrid: 99.73 % joint per branch / deterministic", "cdhim_coord": "hybrid: 99.73 % joint per branch / deterministic",
                     "ac_mc_central_99.73": "empirical", "ac_mc_full_range": "empirical"}
        t_method = {"affine_box": t_lin + t_aff, "cltm_3sigma": t_lin + t_cltm, "cdhim_lmax": t_lin + t_cd, "cdhim_coord": t_lin + t_cd,
                    "ac_mc_central_99.73": t_mc, "ac_mc_full_range": t_mc}
        for name, r in radii.items():
            out = np.abs(dv) > r[:, None] + VIOLATION_ATOL
            rate = out.mean(axis=1)
            rows.append({"system": key, "delta": delta, "method": name, "guarantee": guarantee[name],
                         "mean_width_pu": float(2 * r[ok].mean()), "median_width_pu": float(2 * np.median(r[ok])),
                         "median_ratio_to_ac_99.73_range": float(np.median(r[ok] / r_mc[ok])),
                         "mean_ratio_to_ac_99.73_range": float(np.mean(r[ok] / r_mc[ok])),
                         "bus_ac_violation_mean_pct": pct(rate.mean()), "bus_ac_violation_max_pct": pct(rate.max()),
                         "n_buses_over_0.27pct": int((rate > 0.0027).sum()), "system_level_violation_pct": pct(out.any(axis=0).mean()),
                         "time_s": t_method[name], "n_ac": dv.shape[1], "n_ac_not_converged": n_fail, "seed": seed})
        trows.append({"system": key, "n_ac": n, "linearisation_and_sensitivities_s": t_lin, "cdhim_branch_stage_s": t_cd,
                      "cdhim_branch_stage_std_s": sd_cd, "cdhim_total_s": t_lin + t_cd, "base_pf_pandapower_s": st.grid.base_pf_seconds,
                      "ac_mc_batched_s": t_mc, "ac_mc_batched_ms_per_sample": 1e3 * t_mc / n,
                      "pandapower_ms_per_pf": 1e3 * t_pp, "pandapower_ms_per_pf_std": 1e3 * sd_pp, "pandapower_n_timed": N_PP,
                      "pandapower_loop_extrapolated_s": t_pp * n,
                      "speedup_vs_batched_mc": t_mc / (t_lin + t_cd), "speedup_vs_pandapower_loop_extrapolated": t_pp * n / (t_lin + t_cd)})
        print(f"[c8] {key} done", flush=True)
    meta = {"n_ac": CFG["experiments"]["n_ac"], "threads": 1,
            "note": "all interval methods share the corrected linear model; time_s of an interval method = linearisation + "
                    "sensitivities + its own stage, excluding the base power flow; the pandapower loop time is an extrapolation "
                    f"from {N_PP} timed power flows"}
    print(save_table("c8_methods", rows, meta))
    print(save_table("c8_timing", trows, meta))


if __name__ == "__main__":
    main()
