"""X10 -- The voltage enclosure on a three-phase unbalanced network.

Network: the low-voltage side of the IEEE European LV test feeder as shipped with pandapower
(``pandapower.networks.ieee_european_lv_asymmetric``): 906 buses at 416 V, 905 line sections given by their
positive- and zero-sequence parameters, 55 single-phase loads connected phase to neutral (wye). The 11/0.416 kV
Dyn transformer is left out: an ideal three-phase source sits at the low-voltage busbar, with the voltage
magnitude of the original source. The model is in phase coordinates (cdhim.threephase): one node per bus and
phase, 3 x 3 line blocks with mutual coupling, neutral at ground potential, loads drawing conj(S / V).

Two versions of the same network are used. In the "merged" one the buses that only join two line sections are
removed (their impedances add; exact), which leaves 205 buses; all the checks are run on it. The "full" one keeps
the 906 buses and is used to confirm that the bounds are the same and to time the method at that size.

Uncertainty, as in X8: single-phase inverters in half disks (dp >= 0, |dp + j dq| <= C), here at every
SOURCE_EVERY-th load of the load table, on the phase of that load, all with the same rating, the total equal to
the total apparent load; the two scenarios of Sou et al. (F: forward, R: reverse flow) and several shares rho_U
of uncertain capacity.

Methods: the enclosure of this repository (cdhim.enclosure.voltage_enclosure) and the discs of Wang, Bernstein,
Le Boudec and Paolone (cdhim.discs). Their IEEE TSG 2018 paper is single-phase; the same condition for
multi-phase networks with wye-connected constant-power injections is Theorem 2.1 of C. Wang's thesis (EPFL 9456,
2019, Chapter 2; also their letter in IEEE TPWRS 32(4), 2017, not read here), which is what ``disc_bound``
evaluates when the grid is over node-phases. The radial method of Sou et al. (cdhim.distflow) is a
single-phase model and does not apply.

Reference: the best values ATTAINED for every node-phase (local search, nine boundary load flows, random samples),
as in X8; not exact extremes.

Writes results/x10_model_validation.csv  the phase-coordinate model against pandapower.runpp_3ph
       results/x10_three_phase.csv       enclosure and discs on the merged network
       results/x10_full_model.csv        the same enclosure on the 906-bus version: equal bounds, time
       results/x10_limits.csv            larger inverter ratings: where the fixed point and the discs stop closing

    python experiments/x10_three_phase.py                  # everything
    python experiments/x10_three_phase.py validation cases # only the named tables (validation, cases, full, limits)
"""
import sys
import warnings
from time import perf_counter

import numpy as np
from scipy.sparse.linalg import splu

from common import CFG, pct, seed_for
from cdhim import SecondOrder, save_table
from cdhim.discs import disc_bound
from cdhim.distflow import SCENARIOS, net_injections
from cdhim.enclosure import voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.powerflow import BatchACSolver, newton_pf, polar_jacobian
from cdhim.threephase import from_pandapower, lv_side, pandapower_phase_voltages
from cdhim.validation import VIOLATION_ATOL

SNAPSHOT = "on_peak_566"                                   # load snapshot used in the experiments
VALIDATION_SNAPSHOTS = ("on_peak_566", "off_peak_1", "off_peak_1440")
BASE_MVA_3PH = 0.1                                         # three-phase power base. The network ships with 100 MVA, which
                                                           # puts the loads at 1e-5 p.u., too close to the solver tolerance
SOURCE_EVERY = 5                                           # an inverter at loads 1, 6, 11, ... of the load table
RATING_FACTORS = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0)            # x10_limits: total rating over total load, scenario R
RHO_U_LIMITS = 0.8
RHO_U = CFG["experiments"]["x8_rho_u"]                     # same shares of uncertain capacity as X8
N_SAMPLES = CFG["experiments"]["x8_samples"]
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]
TOL = CFG["powerflow"]["tol_pu"]
NAME = "European LV, low-voltage side"


def stats(x):
    return pct(np.mean(x)), pct(np.max(x))


def pandapower_net(snapshot):
    import pandapower.networks as pn

    net = lv_side(pn.ieee_european_lv_asymmetric(snapshot))
    net.sn_mva = BASE_MVA_3PH
    return net


def unbalance(g, V):
    """Over the buses: largest spread of |V| between phases [p.u.], largest |V2| / |V1| [%] and largest |V0| / |V1| [%]."""
    v = V[g.node[(g.node >= 0).all(axis=1)]]
    a = np.exp(2j * np.pi / 3)
    v0 = v.sum(axis=1) / 3.0
    v1 = (v[:, 0] + a * v[:, 1] + a * a * v[:, 2]) / 3.0
    v2 = (v[:, 0] + a * a * v[:, 1] + a * v[:, 2]) / 3.0
    return float(np.ptp(np.abs(v), axis=1).max()), pct((np.abs(v2) / np.abs(v1)).max()), pct((np.abs(v0) / np.abs(v1)).max())


def node_of(g):
    return {(int(b), int(p)): i for i, (b, p) in enumerate(zip(g.bus_label, g.phase_of))}


def source_nodes(g, net):
    """Node-phase of every SOURCE_EVERY-th load of the load table. Every load of this network is single-phase."""
    key, nodes = node_of(g), []
    for ld in net.asymmetric_load.itertuples():
        s = np.abs([complex(ld.p_a_mw, ld.q_a_mvar), complex(ld.p_b_mw, ld.q_b_mvar), complex(ld.p_c_mw, ld.q_c_mvar)])
        if np.count_nonzero(s) != 1:
            raise ValueError("expected single-phase loads")
        nodes.append(key[(int(ld.bus), int(np.argmax(s)))])
    nodes = np.array(nodes[::SOURCE_EVERY])
    if len(np.unique(nodes)) != len(nodes):
        raise ValueError("two sources on the same node-phase")
    return nodes


def ratings(g, dg, factor=1.0):
    """Equal ratings whose total is ``factor`` times the total apparent load (the rule of Sou et al., Sec. VI-I)."""
    rating = np.zeros(g.n_bus)
    rating[dg] = factor * abs(g.s_load.sum()) / len(dg)
    return rating


# -------------------------------------------------------------------------------------------- model validation
def validation():
    """The phase-coordinate model against the three-phase power flow of pandapower, on three load snapshots."""
    import pandapower as pp
    import pandapower.networks as pn

    rows = []
    for snap in VALIDATION_SNAPSHOTS:
        original = pn.ieee_european_lv_asymmetric(snap)
        net = pandapower_net(snap)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pp.runpp_3ph(original)
            pp.runpp_3ph(net)
        if not net.converged:
            rows.append({"snapshot": snap, "status": "pandapower.runpp_3ph did not converge"})
            continue
        v_all = pandapower_phase_voltages(net)
        row_of = {int(b): i for i, b in enumerate(net.bus.index.values)}
        v_orig = pandapower_phase_voltages(original)
        o_of = {int(b): i for i, b in enumerate(original.bus.index.values)}
        cut = np.abs(np.abs(v_all) - np.abs(v_orig[[o_of[int(b)] for b in net.bus.index.values]])).max()
        full = None
        for merge in (False, True):
            t0 = perf_counter()
            g = from_pandapower(net, merge=merge, tol=TOL)
            t_build = perf_counter() - t0
            if full is None:
                full = g
            v_pp = v_all[[row_of[int(b)] for b in g.bus_label], g.phase_of]
            res = (v_pp * np.conj(g.Ybus @ v_pp) + g.s_load)[g.pq]                    # power mismatch of pandapower's solution here
            dx = splu(polar_jacobian(g.Ybus, g.V0, g.pvpq, g.pq)).solve(np.concatenate([res.real, res.imag]))
            first_order = np.abs(g.V0[g.pq]) + dx[len(g.pq):]                         # |V| moved by that mismatch, to first order
            spread, vuf, v0 = unbalance(g, g.V0)
            kf = node_of(full)
            same = np.array([kf[(int(b), int(p))] for b, p in zip(g.bus_label, g.phase_of)])
            load = [abs(g.s_load[g.phase_of == p].sum()) * g.base_mva * 1e3 for p in range(3)]
            rows.append({"snapshot": snap, "model": "merged" if merge else "full", "status": "ok", "buses": g.node.shape[0], "node_phases": g.n_bus,
                         "loads": len(net.asymmetric_load), "load_phase_a_kva": load[0], "load_phase_b_kva": load[1], "load_phase_c_kva": load[2],
                         "total_load_kva": abs(g.s_load.sum()) * g.base_mva * 1e3,
                         "max_abs_voltage_difference_pu": float(np.abs(g.V0 - v_pp).max()),
                         "max_magnitude_difference_pu": float(np.abs(np.abs(g.V0) - np.abs(v_pp)).max()),
                         "pandapower_residual_in_this_model_max_mva": float(np.abs(res).max() * g.base_mva),
                         "magnitude_difference_after_removing_that_residual_pu": float(np.abs(np.abs(v_pp[g.pq]) - first_order).max()),
                         "own_residual_max_mva": float(np.abs((g.V0 * np.conj(g.Ybus @ g.V0) + g.s_load)[g.pq]).max() * g.base_mva),
                         "merged_minus_full_max_pu": float(np.abs(g.V0 - full.V0[same]).max()) if merge else 0.0,
                         "phase_spread_max_pu": spread, "negative_sequence_over_positive_max_pct": vuf,
                         "zero_sequence_over_positive_max_pct": v0,
                         "pandapower_unbalance_percent_max": float(net.res_bus_3ph.unbalance_percent.max()),
                         "vmin_pu": float(np.abs(g.V0).min()), "vmax_pu": float(np.abs(g.V0).max()),
                         "leaving_the_transformer_out_changes_magnitudes_by_max_pu": float(cut),
                         "build_and_solve_ms": 1e3 * t_build})
            print(f"[x10] validation {snap} {'merged' if merge else 'full'}: {g.n_bus} node-phases, max |V - V_pandapower| "
                  f"{rows[-1]['max_abs_voltage_difference_pu']:.2e} p.u., after removing pandapower's residual "
                  f"{rows[-1]['magnitude_difference_after_removing_that_residual_pu']:.2e}, phase spread {spread:.4f} p.u.", flush=True)
    return rows


# ------------------------------------------------------------------------------------------------- experiments
def one_case(g, dg, s0, cap, seed, search=True):
    """Enclosure, attained range and discs for one set of injections; returns the row and the enclosure."""
    row = {}
    t0 = perf_counter()
    try:
        st = half_disk_study(g, s0, cap, tol=TOL)
    except RuntimeError as e:
        return {"status": str(e)}, None
    so = SecondOrder(st, 1.0)
    enc = voltage_enclosure(so, first_order_range(so.c_lin), node_budget=NODE_BUDGET)
    t_enc = perf_counter() - t0
    n_dg = len(dg)
    v_hat, ok, _, _ = newton_pf(g.Ybus, s0, g.V0, g.pv, g.pq, tol=TOL, max_iter=40)
    k = np.arange(9)
    S = np.repeat(s0[:, None], 9 + N_SAMPLES, axis=1)
    S[:, :9] += cap[:, None] * (np.sin(k * np.pi / 8) + 1j * np.cos(k * np.pi / 8))[None, :]
    S[dg, 9:] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(seed), n_dg, N_SAMPLES)
    V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, st.V_star, tol=TOL).solve(S)
    vs = np.abs(V[g.pq][:, conv])
    vmax, vmin = vs.max(axis=1), vs.min(axis=1)
    t_search = np.nan
    if search:
        t0 = perf_counter()
        vmax = np.maximum(attained_extreme(g, s0, cap, +1, tol=TOL), vmax)
        vmin = np.minimum(attained_extreme(g, s0, cap, -1, tol=TOL), vmin)
        t_search = perf_counter() - t0
    width = vmax - vmin
    corner = np.abs(V[g.pq, 8]) if conv[8] else np.full(len(g.pq), np.nan)             # (dp, dq) = (0, -C) at every source
    gap = corner - vmin
    spread, vuf, v0 = unbalance(g, v_hat) if ok else (np.nan, np.nan, np.nan)
    row.update({"status": "ok", "ac_failures_in_samples": int((~conv).sum()), "nominal_phase_spread_max_pu": spread,
                "nominal_negative_sequence_over_positive_max_pct": vuf, "nominal_zero_sequence_over_positive_max_pct": v0,
                "attained_width_max_pu": float(width.max()), "attained_width_min_pu": float(width.min()),
                "attained_vmin_pu": float(vmin.min()), "attained_vmax_pu": float(vmax.max()),
                "node_phases_with_minimum_not_at_zero_p_minus_C": int((gap > 1e-7).sum()),
                "largest_gap_to_that_minimum_pu": float(np.nanmax(np.maximum(gap, 0.0))),
                "largest_gap_over_attained_width_pct": pct(np.nanmax(np.maximum(gap, 0.0) / width)),
                "enclosure_fixed_point": bool(enc.converged), "enclosure_rel_radius": enc.rel_radius, "enclosure_ms": 1e3 * t_enc,
                "attained_search_s": t_search,
                "second_order_only_upper_violations": int((enc.hi2 < vmax - VIOLATION_ATOL).sum()),
                "second_order_only_lower_violations": int((enc.lo2 > vmin + VIOLATION_ATOL).sum())})
    if enc.converged:
        a, m = stats((enc.hi - vmax) / vmax)
        a2, m2 = stats((vmin - enc.lo) / vmin)
        ratio = (enc.hi - enc.lo) / width
        ph = g.phase_of[g.pq]
        row.update({"enclosure_upper_excess_ave_pct": a, "enclosure_upper_excess_max_pct": m, "enclosure_lower_slack_ave_pct": a2,
                    "enclosure_lower_slack_max_pct": m2, "enclosure_width_over_attained_median": float(np.median(ratio)),
                    "enclosure_width_over_attained_max": float(ratio.max()),
                    "enclosure_width_over_attained_median_phase_a": float(np.median(ratio[ph == 0])),
                    "enclosure_width_over_attained_median_phase_b": float(np.median(ratio[ph == 1])),
                    "enclosure_width_over_attained_median_phase_c": float(np.median(ratio[ph == 2])),
                    "enclosure_allowance_max_pu": float(enc.allowance.max()),
                    "enclosure_allowance_over_width_max_pct": pct((2.0 * enc.allowance / (enc.hi - enc.lo)).max()),
                    "enclosure_violations": int(((enc.hi < vmax - VIOLATION_ATOL) | (enc.lo > vmin + VIOLATION_ATOL)).sum())})
    t0 = perf_counter()
    db = disc_bound(g, s0, v_hat, cap) if ok else None
    row["discs_ms"] = 1e3 * (perf_counter() - t0)
    row["discs_ok"] = bool(db is not None and db.ok)
    if row["discs_ok"]:
        vh = np.abs(v_hat[g.pq])
        for tag, rad in (("discs", db.radius), ("discs_per_bus", db.radius_bus)):
            a, m = stats((vh + rad - vmax) / vmax)
            a2, m2 = stats((vmin - (vh - rad)) / vmin)
            r = 2 * rad / width
            row.update({f"{tag}_upper_excess_ave_pct": a, f"{tag}_upper_excess_max_pct": m, f"{tag}_lower_slack_ave_pct": a2,
                        f"{tag}_lower_slack_max_pct": m2, f"{tag}_width_over_attained_median": float(np.median(r)),
                        f"{tag}_width_over_attained_max": float(r.max()),
                        f"{tag}_violations": int(((vh + rad < vmax - VIOLATION_ATOL) | (vh - rad > vmin + VIOLATION_ATOL)).sum())})
        row["discs_rho"] = db.rho
    return row, enc


def header(g, dg, cap, tag, **extra):
    spread, vuf, v0 = unbalance(g, g.V0)
    return {"network": NAME, "model": tag, "buses": g.node.shape[0], "pq_node_phases": len(g.pq), "sources": len(dg), "symbols": 2 * len(dg),
            "sources_on_phases_a_b_c": "/".join(str(int(x)) for x in np.bincount(g.phase_of[dg], minlength=3)),
            "uncertain_capacity_over_load": float(cap.sum() / abs(g.s_load.sum())),
            "base_case_phase_spread_max_pu": spread, "base_case_negative_sequence_over_positive_max_pct": vuf,
            "base_case_zero_sequence_over_positive_max_pct": v0, **extra}


def cases():
    net = pandapower_net(SNAPSHOT)
    g = from_pandapower(net, merge=True, tol=TOL)
    dg = source_nodes(g, net)
    rating = ratings(g, dg)
    rows = []
    for scen in "FR":
        sc = SCENARIOS[scen]
        for rho_U in RHO_U:
            s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
            seed = seed_for("x10", scen, rho_U)
            row = header(g, dg, cap, "merged", scenario=scen, rho_U=rho_U, seed=seed)
            res, enc = one_case(g, dg, s0, cap, seed)
            row.update(res)
            rows.append(row)
            print(f"[x10] {scen} rho_U={rho_U:g} K={2 * len(dg)}: " + (f"fixed point {row['enclosure_fixed_point']}"
                  + (f" width/attained {row['enclosure_width_over_attained_median']:.4f} (max {row['enclosure_width_over_attained_max']:.4f})"
                     f" violations {row['enclosure_violations']}" if row["enclosure_fixed_point"] else "")
                  + (f" | discs per bus {row['discs_per_bus_width_over_attained_median']:.2f}" if row["discs_ok"] else " | discs fail")
                  + f" | minimum not at (0, -C): {row['node_phases_with_minimum_not_at_zero_p_minus_C']} of {len(g.pq)}"
                  if row["status"] == "ok" else row["status"]), flush=True)
    return rows


def full_model():
    """The enclosure on the 906-bus version against the merged one at the node-phases they share."""
    net = pandapower_net(SNAPSHOT)
    g, gf = from_pandapower(net, merge=True, tol=TOL), from_pandapower(net, merge=False, tol=TOL)
    dg = source_nodes(g, net)
    dgf = source_nodes(gf, net)
    kf = node_of(gf)
    same = np.array([kf[(int(b), int(p))] for b, p in zip(g.bus_label, g.phase_of)])
    pos = np.full(gf.n_bus, -1)
    pos[gf.pq] = np.arange(len(gf.pq))
    o = pos[same[g.pq]]
    rows = []
    for scen in "FR":
        sc = SCENARIOS[scen]
        for rho_U in RHO_U:
            seed = seed_for("x10", "full", scen, rho_U)
            s0, cap = net_injections(g, ratings(g, dg), sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
            s0f, capf = net_injections(gf, ratings(gf, dgf), sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
            enc = voltage_enclosure(SecondOrder(half_disk_study(g, s0, cap, tol=TOL), 1.0), None, node_budget=NODE_BUDGET)
            row = header(gf, dgf, capf, "full", scenario=scen, rho_U=rho_U, seed=seed)
            t0 = perf_counter()
            st = half_disk_study(gf, s0f, capf, tol=TOL)
            t1 = perf_counter()
            so = SecondOrder(st, 1.0)
            t2 = perf_counter()
            encf = voltage_enclosure(so, None, node_budget=NODE_BUDGET)
            t3 = perf_counter()
            fo = first_order_range(so.c_lin)
            r = np.abs(so.c_lin).sum(axis=1)
            lo, hi = encf.lo + r + fo[0], encf.hi - r + fo[1]                         # the same bounds with the half-disk first order
            S = np.repeat(s0f[:, None], N_SAMPLES, axis=1)
            S[dgf] += capf[dgf][:, None] * sample_half_disks(np.random.default_rng(seed), len(dgf), N_SAMPLES)
            V, conv, _ = BatchACSolver(gf.Ybus, gf.pv, gf.pq, st.V_star, tol=TOL).solve(S)
            vs = np.abs(V[gf.pq][:, conv])
            row.update({"enclosure_fixed_point": bool(encf.converged), "merged_fixed_point": bool(enc.converged),
                        "study_ms": 1e3 * (t1 - t0), "second_order_ms": 1e3 * (t2 - t1), "enclosure_ms": 1e3 * (t3 - t2),
                        "total_ms": 1e3 * (t3 - t0), "samples": int(conv.sum()), "ac_failures_in_samples": int((~conv).sum())})
            if encf.converged and enc.converged:
                row.update({"bounds_minus_merged_model_max_pu": float(max(np.abs(encf.lo[o] - enc.lo).max(), np.abs(encf.hi[o] - enc.hi).max())),
                            "allowance_minus_merged_model_max_pu": float(np.abs(encf.allowance[o] - enc.allowance).max()),
                            "violations_in_samples": int(((hi < vs.max(axis=1) - VIOLATION_ATOL) | (lo > vs.min(axis=1) + VIOLATION_ATOL)).sum())})
            rows.append(row)
            print(f"[x10] full model {scen} rho_U={rho_U:g}: {len(gf.pq)} PQ node-phases, fixed point {encf.converged}, {row['total_ms'] / 1e3:.1f} s,"
                  f" bounds minus merged {row.get('bounds_minus_merged_model_max_pu', float('nan')):.1e} p.u.", flush=True)
    return rows


def limits():
    """Scenario R with rho_U = RHO_U_LIMITS and the inverter ratings multiplied: where the methods stop closing."""
    net = pandapower_net(SNAPSHOT)
    g = from_pandapower(net, merge=True, tol=TOL)
    dg = source_nodes(g, net)
    sc = SCENARIOS["R"]
    rows = []
    for factor in RATING_FACTORS:
        s0, cap = net_injections(g, ratings(g, dg, factor), sc["rho_D"], sc["rho_G"], sc["pf"], RHO_U_LIMITS)
        seed = seed_for("x10", "limits", factor)
        row = header(g, dg, cap, "merged", scenario="R", rho_U=RHO_U_LIMITS, total_rating_over_load=factor, seed=seed)
        res, enc = one_case(g, dg, s0, cap, seed)
        row.update(res)
        rows.append(row)
        print(f"[x10] limits rating x{factor:g}: " + (f"fixed point {row['enclosure_fixed_point']} (radius {row['enclosure_rel_radius']:.3f}),"
              f" discs {row['discs_ok']}, attained |V| in [{row['attained_vmin_pu']:.3f}, {row['attained_vmax_pu']:.3f}]"
              + (f", width/attained {row['enclosure_width_over_attained_median']:.3f}" if row["enclosure_fixed_point"] else "")
              if row["status"] == "ok" else row["status"]), flush=True)
    return rows


PARTS = {"validation": ("x10_model_validation", validation), "cases": ("x10_three_phase", cases),
         "full": ("x10_full_model", full_model), "limits": ("x10_limits", limits)}


def main():
    meta = {"seed_base": CFG["experiments"]["seed"], "samples_per_case": N_SAMPLES, "rho_U": RHO_U, "tolerance_pu": TOL,
            "raak": "all symbols in one block search, node budget %d, forms normalised" % NODE_BUDGET,
            "network": "pandapower.networks.ieee_european_lv_asymmetric('%s') without its transformer: ideal source at the "
                       "low-voltage busbar; phase-coordinate model of cdhim.threephase" % SNAPSHOT,
            "base_mva_three_phase": BASE_MVA_3PH,
            "sources": "single-phase inverter at every %d-th load of the load table, on the phase of the load; equal ratings, "
                       "total = |sum of the loads| (times the factor in x10_limits)" % SOURCE_EVERY,
            "scenarios": SCENARIOS, "reference": "attained values (local search + 9 boundary load flows + random samples), not exact extremes",
            "limits": {"rating_factors": RATING_FACTORS, "rho_U": RHO_U_LIMITS, "scenario": "R"},
            "validation_snapshots": VALIDATION_SNAPSHOTS,
            "note": "times are single runs on a shared machine: indicative only"}
    wanted = sys.argv[1:] or list(PARTS)
    unknown = [w for w in wanted if w not in PARTS]
    if unknown:
        sys.exit(f"unknown table(s) {unknown}; choose from {list(PARTS)}")
    for part in wanted:
        table, fn = PARTS[part]
        print(save_table(table, fn(), meta))


if __name__ == "__main__":
    main()
