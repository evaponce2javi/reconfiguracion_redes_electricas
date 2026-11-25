"""X8 -- Verified voltage enclosures on distribution networks, against the published baselines.

Uncertainty: inverter-based generation in half disks (dp >= 0, |dp + j dq| <= C per source), added with the
benchmark rule of Sou, Khong and Sandberg (IEEE TPWRS 2026, Sec. VI-I) and their two scenarios, for several
shares rho_U of uncertain capacity.

Methods
  enclosure  second-order affine model + RAA-K bound of the quadratic term + fixed-point bound of the remainder
             in current-injection form (cdhim.enclosure.voltage_enclosure); any topology, PQ buses only.
  WCLF       Sou et al., reimplemented (cdhim.distflow): radial networks without shunt elements.
  discs      Wang, Bernstein, Le Boudec, Paolone (IEEE TSG 2018), Theorem 1, and its per-bus refinement
             (cdhim.discs): any topology, PQ buses only.

Reference: there is no exact maximum or minimum to compare with. For every bus the best values ATTAINED are
used (local maximisation and minimisation over the half disks, the nine boundary load flows of the WCLF method,
and random samples). A bound can be at most as loose as its distance to the attained value; a bound on the
wrong side of an attained value would be an error.

Writes results/x8_wclf_validation.csv  the reimplemented WCLF against Tables I, VI and VII of its paper
       results/x8_feeders.csv          the three methods on radial, meshed and shunt-compensated networks
       results/x8_load_boxes.csv       boxes of load uncertainty (as in X6): the two bounds of the remainder
       results/x8_raak_scale.csv       RAA-K with and without normalising the form, in the settings of X2 to X7
       results/x8_remainder_bounds.csv the two bounds of the remainder on every case of x8_feeders.csv
       results/x8_monotonicity.csv     is the minimum voltage of every bus at (dp, dq) = (0, -C)?

    python experiments/x8_feeders.py                 # everything
    python experiments/x8_feeders.py remainder scale # only the named tables (validation, feeders, boxes, scale, remainder, monotone)
"""
import sys
from time import perf_counter

import numpy as np

from common import CFG, grid, pct, seed_for, study
from cdhim import SecondOrder, load_grid, make_study, save_table, spectral_zones
from cdhim.discs import disc_bound
from cdhim.distflow import SCENARIOS, benchmark_ratings, feeder, net_injections, worst_case_load_flow
from cdhim.enclosure import Enclosure, InjectionEnclosure, voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.knomial import knomial_bound
from cdhim.powerflow import BatchACSolver, newton_pf
from cdhim.validation import VIOLATION_ATOL
from cdhim.zones import laplacian_embedding
from cdhim.second_order import bound_monomial
from x2_few_sources import R_FLOOR, vertex_batches

NODE_BUDGET = CFG["experiments"]["raak_node_budget"]
RHO_U = CFG["experiments"]["x8_rho_u"]
N_SAMPLES = CFG["experiments"]["x8_samples"]
TOL = CFG["powerflow"]["tol_pu"]
START_AMP = 5000.0                       # initial current bound of Algorithm 1, as suggested in the paper


def close_ties(net):
    net.line["in_service"] = True


def close_switches(net):
    net.switch["closed"] = True


def add_shunts(net):
    """Own modification of the 33-bus feeder: line charging on every line and two capacitor banks at feeder ends."""
    import pandapower as pp

    net.line["c_nf_per_km"] = 300.0
    for bus in (17, 32):
        pp.create_shunt(net, bus, q_mvar=-0.3)


# (name in the tables, case, modification, description)
NETWORKS = [("33", "matpower_case33bw", None, "MATPOWER 7.1 case33bw, radial"),
            ("69", "matpower_case69", None, "MATPOWER 7.1 case69, radial"),
            ("85", "matpower_case85", None, "MATPOWER 7.1 case85, radial"),
            ("141", "matpower_case141", None, "MATPOWER 7.1 case141, radial"),
            ("33 + ties", "matpower_case33bw", close_ties, "case33bw with its five tie lines closed"),
            ("33 + shunts", "matpower_case33bw", add_shunts, "case33bw, radial, with line charging (300 nF/km) and two 0.3 Mvar capacitor banks (own modification)"),
            ("CIGRE MV", "create_cigre_network_mv", close_switches, "pandapower CIGRE MV network, all switches closed: cables with charging, two 110/20 kV transformers"),
            ("Oberrhein", "mv_oberrhein", close_switches, "pandapower MV Oberrhein network, all switches closed: two 110/20 kV transformers, two slack buses")]

# Published values [%]: eps_ave and eps_max of (vU, vL). Table I: 33-bus feeder with the five sources of Fig. 10.
TABLE_1 = {0.0: ("7.5e-09", "1.3e-08"), 0.5: ("0.086", "0.17"), 1.0: ("0.18", "0.34")}
FIG_10 = ((9, 400), (17, 1000), (21, 500), (24, 800), (32, 1200))                 # (bus, kVA); bus 0 = substation
TABLES_6_7 = {"33": {"F": ("0.065", "0.14"), "R": ("0.074", "0.14")}, "69": {"F": ("0.012", "0.026"), "R": ("0.020", "0.049")},
              "85": {"F": ("0.14", "0.19"), "R": ("0.17", "0.20")}, "141": {"F": ("0.051", "0.084"), "R": ("0.058", "0.074")}}


def stats(x):
    return pct(np.mean(x)), pct(np.max(x))


def validation():
    """The reimplemented WCLF against the numbers printed in its paper."""
    rows = []
    g = load_grid("matpower_case33bw", tol=TOL)
    lab = np.asarray(g.bus_label)
    rating = np.zeros(g.n_bus)
    for bus, kva in FIG_10:
        rating[np.flatnonzero(lab == bus)[0]] = kva / 1e3 / g.base_mva
    for rho_U, (ave, top) in TABLE_1.items():
        s0, cap = net_injections(g, rating, 1.0, 0.5, 0.8, rho_U)
        wc = worst_case_load_flow(g, s0, cap, current_start_amp=START_AMP)
        a, m = stats((wc.v_max_U - wc.v_max_L) / wc.v_max_L)
        rows.append({"table": "I", "network": "33", "setup": "five sources of Fig. 10; rho_D = 1, rho_G = 0.5, pf = 0.8", "rho_U": rho_U,
                     "eps_ave_published_pct": ave, "eps_ave_pct": a, "eps_max_published_pct": top, "eps_max_pct": m,
                     "converged": wc.converged, "iterations": wc.iterations})
    for name, case, _, _ in NETWORKS[:4]:
        g = load_grid(case, tol=TOL)
        rating = benchmark_ratings(g)
        for scen, table in (("F", "VI"), ("R", "VII")):
            sc = SCENARIOS[scen]
            s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], sc["rho_U"])
            wc = worst_case_load_flow(g, s0, cap, current_start_amp=START_AMP)
            a, m = stats((wc.v_max_U - wc.v_max_L) / wc.v_max_L)
            rows.append({"table": table, "network": name, "setup": f"source at bus 1 + 10k, total rating = |sum of loads|; scenario {scen}", "rho_U": sc["rho_U"],
                         "eps_ave_published_pct": TABLES_6_7[name][scen][0], "eps_ave_pct": a, "eps_max_published_pct": TABLES_6_7[name][scen][1],
                         "eps_max_pct": m, "converged": wc.converged, "iterations": wc.iterations})
    return rows


def run_network(name, case, prepare, description):
    g = load_grid(case, tol=TOL, prepare=prepare)
    edges = {(min(a, b), max(a, b)) for a, b in zip(g.br_f, g.br_t) if a != b}
    try:
        feeder(g)
        wclf_reason = ""
    except ValueError as e:
        wclf_reason = str(e)
    rating = benchmark_ratings(g)
    dg = np.flatnonzero(rating > 0)
    n_dg = len(dg)
    rows = []
    for scen in "FR":
        sc = SCENARIOS[scen]
        for rho_U in RHO_U:
            s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
            seed = seed_for("x8", name, scen, rho_U)
            row = {"network": name, "scenario": scen, "rho_U": rho_U, "n_bus": g.n_bus, "loops": len(edges) - g.n_bus + 1, "sources": n_dg, "symbols": 2 * n_dg,
                   "uncertain_capacity_over_load": float(cap.sum() / abs(g.s_load.sum())), "seed": seed}
            # ---- enclosure of this repository
            t0 = perf_counter()
            try:
                st = half_disk_study(g, s0, cap, tol=TOL)
            except RuntimeError as e:
                row["status"] = str(e)
                rows.append(row)
                continue
            so = SecondOrder(st, 1.0)
            enc = voltage_enclosure(so, first_order_range(so.c_lin), node_budget=NODE_BUDGET)
            t_enc = perf_counter() - t0
            # ---- attained extremes: local search, boundary load flows, random samples
            v_hat, ok, _, _ = newton_pf(g.Ybus, s0, g.V0, g.pv, g.pq, tol=TOL, max_iter=40)
            k = np.arange(9)
            S = np.repeat(s0[:, None], 9 + N_SAMPLES, axis=1)
            S[:, :9] += cap[:, None] * (np.sin(k * np.pi / 8) + 1j * np.cos(k * np.pi / 8))[None, :]
            S[dg, 9:] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(seed), n_dg, N_SAMPLES)
            V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, st.V_star, tol=TOL).solve(S)
            vs = np.abs(V[g.pq][:, conv])
            corner = np.abs(V[g.pq, 8]) if conv[8] else np.full(len(g.pq), np.nan)       # (dp, dq) = (0, -C)
            vmax = np.maximum(attained_extreme(g, s0, cap, +1, tol=TOL), vs.max(axis=1))
            vmin = np.minimum(attained_extreme(g, s0, cap, -1, tol=TOL), vs.min(axis=1))
            width = vmax - vmin
            row.update({"status": "ok", "ac_failures_in_samples": int((~conv).sum()), "attained_width_max_pu": float(width.max()),
                        "attained_vmin_pu": float(vmin.min()), "attained_vmax_pu": float(vmax.max()),
                        "buses_with_minimum_not_at_zero_p_minus_C": int((corner > vmin + 1e-7).sum()),
                        "enclosure_fixed_point": bool(enc.converged), "enclosure_rel_radius": enc.rel_radius, "enclosure_ms": 1e3 * t_enc,
                        "second_order_only_upper_violations": int((enc.hi2 < vmax - VIOLATION_ATOL).sum()),
                        "second_order_only_lower_violations": int((enc.lo2 > vmin + VIOLATION_ATOL).sum()),
                        "second_order_only_lower_violation_max_pu": float(np.maximum(enc.lo2 - vmin, 0.0).max())})
            if enc.converged:
                a, m = stats((enc.hi - vmax) / vmax)
                a2, m2 = stats((vmin - enc.lo) / vmin)
                ratio = (enc.hi - enc.lo) / width
                row.update({"enclosure_upper_excess_ave_pct": a, "enclosure_upper_excess_max_pct": m, "enclosure_lower_slack_ave_pct": a2,
                            "enclosure_lower_slack_max_pct": m2, "enclosure_width_over_attained_median": float(np.median(ratio)),
                            "enclosure_width_over_attained_max": float(ratio.max()),
                            "enclosure_violations": int(((enc.hi < vmax - VIOLATION_ATOL) | (enc.lo > vmin + VIOLATION_ATOL)).sum())})
            # ---- WCLF (radial, no shunts)
            row["wclf_not_applicable"] = wclf_reason
            if not wclf_reason:
                t0 = perf_counter()
                wc = worst_case_load_flow(g, s0, cap, current_start_amp=START_AMP)
                t_wc = perf_counter() - t0
                pos = np.full(g.n_bus, -1)
                pos[g.pq] = np.arange(len(g.pq))
                o = pos[wc.nodes]
                a, m = stats((wc.v_max_U - wc.v_max_L) / wc.v_max_L)
                a2, m2 = stats((wc.v_max_U - vmax[o]) / vmax[o])
                row.update({"wclf_converged": wc.converged, "wclf_iterations": wc.iterations, "wclf_ms": 1e3 * t_wc,
                            "wclf_certified_gap_ave_pct": a, "wclf_certified_gap_max_pct": m,
                            "wclf_upper_excess_ave_pct": a2, "wclf_upper_excess_max_pct": m2,
                            "wclf_vmin_minus_attained_max_pu": float(np.abs(wc.v_min - vmin[o]).max()),
                            "wclf_violations": int(((wc.v_max_U < vmax[o] - VIOLATION_ATOL) | (wc.v_min > vmin[o] + VIOLATION_ATOL)).sum()),
                            "attained_minus_wclf_lower_bound_ave_pct": pct(np.mean((vmax[o] - wc.v_max_L) / wc.v_max_L))})
                if enc.converged:
                    row["enclosure_upper_bound_below_wclf_buses"] = int((enc.hi[o] < wc.v_max_U).sum())
            # ---- discs (centre: the nominal injections, |ds| <= C)
            t0 = perf_counter()
            db = disc_bound(g, s0, v_hat, cap) if ok else None
            row["discs_ms"] = 1e3 * (perf_counter() - t0)
            row["discs_ok"] = bool(db is not None and db.ok)
            if row["discs_ok"]:
                vh = np.abs(v_hat[g.pq])
                for tag, rad in (("discs", db.radius), ("discs_per_bus", db.radius_bus)):
                    a, m = stats((vh + rad - vmax) / vmax)
                    a2, m2 = stats((vmin - (vh - rad)) / vmin)
                    row.update({f"{tag}_upper_excess_ave_pct": a, f"{tag}_upper_excess_max_pct": m, f"{tag}_lower_slack_ave_pct": a2,
                                f"{tag}_lower_slack_max_pct": m2, f"{tag}_width_over_attained_median": float(np.median(2 * rad / width)),
                                f"{tag}_violations": int(((vh + rad < vmax - VIOLATION_ATOL) | (vh - rad > vmin + VIOLATION_ATOL)).sum())})
                row["discs_rho"] = db.rho
            rows.append(row)
            print(f"[x8] {name} {scen} rho_U={rho_U:g} K={2 * n_dg}: enclosure fixed point {enc.converged}"
                  + (f" width/attained {row['enclosure_width_over_attained_median']:.3f} upper excess {row['enclosure_upper_excess_ave_pct']:.4f} %"
                     f" violations {row['enclosure_violations']}" if enc.converged else "")
                  + (f" | WCLF upper excess {row['wclf_upper_excess_ave_pct']:.4f} % violations {row['wclf_violations']}" if not wclf_reason else " | WCLF n/a")
                  + (f" | discs width/attained {row['discs_per_bus_width_over_attained_median']:.2f}" if row["discs_ok"] else " | discs fail"), flush=True)
    return rows


def load_boxes():
    """Boxes of load uncertainty in 8 zones (the setting of X6) on the 33-bus feeder: the two bounds of the remainder."""
    rows = []
    for name, case, prepare, _ in (NETWORKS[0], NETWORKS[4]):
        g = load_grid(case, tol=TOL, prepare=prepare)
        st = make_study(g, zones=spectral_zones(g, 8, laplacian_embedding(g)), tol=TOL)
        for delta in CFG["experiments"]["x8_load_deltas"]:
            so = SecondOrder(st, delta)
            enc = voltage_enclosure(so, node_budget=NODE_BUDGET)
            old = Enclosure(so)
            fp_old = old.solve(1.0)
            vmin, vmax = np.full(len(g.pq), np.inf), np.full(len(g.pq), -np.inf)
            fails = 0
            for eps in vertex_batches(st.src.n):
                V, conv, _ = st.solver.solve(st.injections(eps, delta))
                fails += int((~conv).sum())
                vm = np.abs(V[g.pq][:, conv])
                vmin, vmax = np.minimum(vmin, vm.min(axis=1)), np.maximum(vmax, vm.max(axis=1))
            act = np.abs(so.c_lin).sum(axis=1) > 1e-9
            width = (vmax - vmin)[act]
            row = {"network": name, "K": st.src.n, "delta": delta, "vertices_without_ac_solution": fails,
                   "vertex_range_width_max_pu": float(width.max()), "rectangular_fixed_point": bool(fp_old.converged),
                   "injection_fixed_point": bool(enc.converged)}
            if enc.converged:
                ratio = ((enc.hi - enc.lo)[act]) / width
                row.update({"injection_width_over_vertex_range_median": float(np.median(ratio)), "injection_width_over_vertex_range_max": float(ratio.max()),
                            "injection_violations": int(((enc.hi < vmax - VIOLATION_ATOL) | (enc.lo > vmin + VIOLATION_ATOL))[act].sum()),
                            "injection_allowance_max_pu": float(enc.allowance[act].max())})
            if fp_old.converged:
                allow_old = old.modulus_allowance(fp_old, 1.0)[0]
                ratio = ((enc.hi2 - enc.lo2 + 2 * allow_old)[act]) / width
                row.update({"rectangular_width_over_vertex_range_median": float(np.median(ratio)), "rectangular_allowance_max_pu": float(allow_old[act].max())})
            rows.append(row)
            print(f"[x8] load boxes {name} delta={delta:g}: rectangular {fp_old.converged} injection {enc.converged}"
                  + (f" width/range {row['injection_width_over_vertex_range_median']:.4f}" if enc.converged else ""), flush=True)
    return rows


def raak_scale():
    """The RAA-K bound of the same quadratic forms with and without normalising them (zones of load uncertainty,
    as in X2 to X7). The search uses an absolute tolerance, so small forms lose their blocks unless normalised."""
    rows = []
    for key, ks, deltas in (("ieee33", (8, 16), (0.05, 0.15, 0.30)), ("pegase2869", (8, 16), (0.01, 0.02, 0.05))):
        case = CFG["systems"][key]["case"]
        g = grid(case)
        emb = laplacian_embedding(g)
        for k in ks:
            st = study(case, zones=spectral_zones(g, k, emb))
            n = st.src.n
            eps = np.hstack(list(vertex_batches(n))) if n <= 8 else None
            for delta in deltas:
                H = SecondOrder(st, delta).forms()
                H = H[np.abs(st.coefficients(delta)).sum(axis=1) > R_FLOOR]
                raw = np.array([knomial_bound(h, n, NODE_BUDGET)[:2] for h in H])
                nrm = np.array([knomial_bound(h, n, NODE_BUDGET, normalise=True)[:2] for h in H])
                af2 = np.array([bound_monomial(h) for h in H])
                w = lambda b: b[:, 1] - b[:, 0]
                ok = w(af2) > 0
                change = (w(raw)[ok] - w(nrm)[ok]) / w(raw)[ok]
                row = {"system": key, "K": n, "delta": delta, "design_delta": bool(np.isclose(delta, CFG["systems"][key]["delta"])),
                       "forms": int(ok.sum()), "largest_entry_median": float(np.median(np.abs(H).max(axis=(1, 2)))),
                       "raw_over_af2_median": float(np.median(w(raw)[ok] / w(af2)[ok])),
                       "normalised_over_af2_median": float(np.median(w(nrm)[ok] / w(af2)[ok])),
                       "forms_narrower_by_over_0p1pct": int((change > 1e-3).sum()), "largest_narrowing_pct": pct(change.max()),
                       "normalised_wider_than_raw": int((w(nrm) > w(raw) * (1 + 1e-12)).sum())}
                if eps is not None:
                    val = np.einsum("kv,ikl,lv->iv", eps, H, eps)
                    ref = val.max(axis=1) - val.min(axis=1)
                    good = ok & (ref > 0)
                    row.update({"raw_over_vertex_range_median": float(np.median(w(raw)[good] / ref[good])),
                                "normalised_over_vertex_range_median": float(np.median(w(nrm)[good] / ref[good])),
                                "normalised_short_of_vertex_range": int(((nrm[:, 0] > val.min(axis=1) + 1e-15) | (nrm[:, 1] < val.max(axis=1) - 1e-15)).sum())})
                rows.append(row)
                print(f"[x8] raak scale {key} K={n} delta={delta:g}: raw {row['raw_over_af2_median']:.3f} normalised {row['normalised_over_af2_median']:.3f} of AF2", flush=True)
    return rows


def remainder_bounds():
    """Section 22's bound of the remainder (rectangular coordinates) and the current-injection one, case by case."""
    rows = []
    for name, case, prepare, _ in NETWORKS:
        g = load_grid(case, tol=TOL, prepare=prepare)
        Y = abs(g.Ybus).tocoo()
        y_max = float(Y.data[Y.row != Y.col].max())
        rating = benchmark_ratings(g)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                row = {"network": name, "scenario": scen, "rho_U": rho_U, "largest_branch_admittance_pu": y_max}
                try:
                    so = SecondOrder(half_disk_study(g, s0, cap, tol=TOL), 1.0)
                except RuntimeError as e:
                    row["status"] = str(e)
                    rows.append(row)
                    continue
                old = Enclosure(so)
                fp_old = old.solve(1.0)
                fp_new = InjectionEnclosure(so, old.Hs).solve(1.0)
                m = so.m
                row.update({"status": "ok", "rectangular_spectral_radius_G": fp_old.spectral_radius, "rectangular_fixed_point": bool(fp_old.converged),
                            "injection_fixed_point": bool(fp_new.converged)})
                if fp_old.converged:
                    row["rectangular_third_order_bound_max_pu"] = float(np.hypot(fp_old.sigma[:m], fp_old.sigma[m:]).max())
                if fp_new.converged:
                    row["injection_third_order_bound_max_pu"] = float(fp_new.s_abs.max())
                rows.append(row)
                print(f"[x8] remainder {name} {scen} rho_U={rho_U:g}: rectangular {fp_old.converged} (radius {fp_old.spectral_radius:.3f}) injection {fp_new.converged}", flush=True)
    return rows


def monotonicity():
    """The WCLF method takes the minimum voltages from one load flow at (dp, dq) = (0, -C), which is exact if the
    voltages are monotone in the injections. Here: the voltage at that point against the lowest value found by
    the local search, bus by bus. A lower value is attained by a feasible injection, so it disproves monotonicity."""
    rows = []
    for name, case, prepare, _ in NETWORKS:
        g = load_grid(case, tol=TOL, prepare=prepare)
        rating = benchmark_ratings(g)
        dg = np.flatnonzero(rating > 0)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                row = {"network": name, "scenario": scen, "rho_U": rho_U, "pq_buses": len(g.pq)}
                V, ok, _, _ = newton_pf(g.Ybus, s0 - 1j * cap, g.V0, g.pv, g.pq, tol=TOL, max_iter=40)
                if not ok:
                    row["status"] = "no AC solution at (0, -C)"
                    rows.append(row)
                    continue
                corner = np.abs(V[g.pq])
                gap = corner - attained_extreme(g, s0, cap, -1, tol=TOL)
                lower = gap > 1e-7
                row.update({"status": "ok", "buses_with_lower_minimum": int(lower.sum()), "largest_gap_pu": float(max(gap.max(), 0.0)),
                            "largest_gap_over_attained_width_pct": pct((gap / np.maximum(attained_extreme(g, s0, cap, +1, tol=TOL) - (corner - gap), 1e-300)).max())
                            if lower.any() else 0.0})
                rows.append(row)
                print(f"[x8] monotone {name} {scen} rho_U={rho_U:g}: {row['buses_with_lower_minimum']} of {len(g.pq)} buses below the corner value, largest gap {row['largest_gap_pu']:.2e} p.u.", flush=True)
    return rows


PARTS = {"validation": ("x8_wclf_validation", lambda: validation()), "feeders": ("x8_feeders", lambda: [r for net in NETWORKS for r in run_network(*net)]),
         "boxes": ("x8_load_boxes", lambda: load_boxes()), "scale": ("x8_raak_scale", lambda: raak_scale()),
         "remainder": ("x8_remainder_bounds", lambda: remainder_bounds()), "monotone": ("x8_monotonicity", lambda: monotonicity())}


def main():
    meta = {"raak_node_budget": NODE_BUDGET, "samples_per_case": N_SAMPLES, "seed_rule": "common.seed_for('x8', network, scenario, rho_U), stored per row",
            "wclf_initial_current_A": START_AMP, "networks": {n: d for n, _, _, d in NETWORKS},
            "note": "bounds valid in exact arithmetic, no outward rounding; reference = best attained values, not exact extremes"}
    wanted = sys.argv[1:] or list(PARTS)
    unknown = [w for w in wanted if w not in PARTS]
    if unknown:
        sys.exit(f"unknown table(s) {unknown}; choose from {list(PARTS)}")
    for part in wanted:
        table, fn = PARTS[part]
        print(save_table(table, fn(), meta))


if __name__ == "__main__":
    main()
