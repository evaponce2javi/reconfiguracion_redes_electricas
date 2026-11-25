"""X12 -- The enclosure with directed rounding, and certified references for it.

    verified   the cases of X8 and the three-phase feeder of X10 (half disks): the verified enclosure
               (cdhim.verified) next to the one of X8 and to the Krawczyk baseline of X9 (cdhim.krawczyk), and
               certified inner bounds: the best attained point of every bus, with its power-flow solution verified.
               Then  (hi - lo) / (vmax_lo - vmin_hi)  is a proven upper bound of width / exact range.
    exact      boxes with few symbols: bracket of the exact extremes by adaptive bisection (cdhim.certified), for
               the load boxes of X8 (8 zones) and for the bounding boxes of the half disks on the 33-bus networks.

    subdivision  the same boxes, bisected until no leaf exceeds the attained values by more than 1 % of the width,
               with the verified enclosure in every leaf and with the per-bus disc of Wang et al. (cdhim.discs):
               subdividing the set is how that method is made tight in its own paper. Enclosures needed by each.

Writes results/x12_verified.csv, results/x12_exact_range.csv and results/x12_subdivision.csv.

    python experiments/x12_verified.py            # both
    python experiments/x12_verified.py exact      # only the named table
"""
import sys
from time import perf_counter

import numpy as np

from common import CFG, pct, seed_for
from cdhim import SecondOrder, load_grid, make_study, save_table, spectral_zones
from cdhim.certified import box_bracket, half_disk_inner, injection_matrix
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import InjectionEnclosure, state_forms, voltage_enclosure
from cdhim.halfdisk import attained_points, first_order_range, half_disk_study, sample_half_disks
from cdhim.krawczyk import Krawczyk
from cdhim.powerflow import BatchACSolver
from cdhim.verified import verified_enclosure
from cdhim.zones import laplacian_embedding
from x8_feeders import N_SAMPLES, NETWORKS, NODE_BUDGET, RHO_U, TOL

BRACKET_TOL = 1e-4              # stop when no leaf exceeds the attained values by more than this share of the width
BRACKET_BOXES = 3000            # enclosures per case at most
SUBDIVISION_TOL = 1e-2          # target of the subdivision comparison: excess of at most 1 % of the attained width
SUBDIVISION_BOXES = 3000        # enclosures per case and method at most
MIN_WIDTH = 1e-9                # buses whose attained range is narrower are left out of the ratios [p.u.]


def attained(g, st, s0, cap, seed):
    """Largest and smallest |V_i| found (local search, nine boundary flows, random samples) and where:
    (vmax, vmin, points_max, points_min), the points as u + j w of the unit half disk of every source."""
    dg = np.flatnonzero(cap > 0)
    n_dg = len(dg)
    k = np.arange(9)
    pts = np.hstack([np.repeat((np.sin(k * np.pi / 8) + 1j * np.cos(k * np.pi / 8))[None, :], n_dg, axis=0),
                     sample_half_disks(np.random.default_rng(seed), n_dg, N_SAMPLES)])
    pts.real = np.maximum(pts.real, 0.0)
    S = np.repeat(np.asarray(s0, dtype=complex)[:, None], pts.shape[1], axis=1)
    S[dg] += cap[dg][:, None] * pts
    V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, st.V_star, tol=TOL).solve(S)
    vs = np.where(conv[None, :], np.abs(V[g.pq]), np.nan)
    out = []
    for sign, pick in ((+1, np.nanargmax), (-1, np.nanargmin)):
        val, point = attained_points(g, s0, cap, sign, tol=TOL)
        col = pick(vs, axis=1)
        samp = vs[np.arange(len(g.pq)), col]
        better = sign * samp > sign * val
        val, point = np.where(better, samp, val), np.where(better[:, None], pts[:, col].T, point)
        out += [val, point]
    return out[0], out[2], out[1], out[3]


THREE_PHASE = "European LV, 3 phases"


def networks():
    """(name, grid, source ratings, tag of the seeds): the networks of X8 and the three-phase feeder of X10."""
    for name, case, prepare, _ in NETWORKS:
        g = load_grid(case, tol=TOL, prepare=prepare)
        yield name, g, benchmark_ratings(g), ("x8", name)
    import x10_three_phase as x10
    from cdhim.threephase import from_pandapower

    net = x10.pandapower_net(x10.SNAPSHOT)
    g = from_pandapower(net, merge=True, tol=TOL)
    yield THREE_PHASE, g, x10.ratings(g, x10.source_nodes(g, net)), ("x10",)


def verified_cases():
    rows = []
    for name, g, rating, tag in networks():
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                seed = seed_for(*tag, scen, rho_U)
                row = {"network": name, "scenario": scen, "rho_U": rho_U, "n_bus": g.n_bus, "pq_nodes": len(g.pq), "symbols": 2 * int((cap > 0).sum()), "seed": seed}
                try:
                    st = half_disk_study(g, s0, cap, tol=TOL)
                except RuntimeError as e:
                    row["status"] = str(e)
                    rows.append(row)
                    continue
                t0 = perf_counter()
                so = SecondOrder(st, 1.0)
                Hs = state_forms(so)
                t_model = perf_counter() - t0
                t0 = perf_counter()
                ref = voltage_enclosure(so, first_order_range(so.c_lin), node_budget=NODE_BUDGET)
                t_ref = perf_counter() - t0
                t0 = perf_counter()
                enc = verified_enclosure(so, half_disks=True, node_budget=NODE_BUDGET, Hs=Hs)
                t_ver = perf_counter() - t0
                sep = verified_enclosure(so, half_disks=True, node_budget=NODE_BUDGET, Hs=Hs, joint=False)
                t0 = perf_counter()
                kr = Krawczyk(so, "current", True, True)
                kra = kr.solve(first_order_range(kr.L))
                t_kr = perf_counter() - t0
                old = InjectionEnclosure(so, Hs).solve(1.0)
                row.update({"status": "ok", "x8_fixed_point": bool(ref.converged), "verified": bool(enc.verified), "unique": bool(enc.unique),
                            "reason": enc.reason, "model_ms": 1e3 * t_model, "x8_enclosure_ms": 1e3 * t_ref, "verified_ms": 1e3 * t_ver,
                            "krawczyk_closes": bool(kra.closes), "krawczyk_ms": 1e3 * t_kr})
                vmax, vmin, pmax, pmin = attained(g, st, s0, cap, seed)
                width = vmax - vmin
                act = width > MIN_WIDTH
                t0 = perf_counter()
                vmin_hi, vmax_lo = half_disk_inner(st, pmax, pmin, TOL)
                row.update({"inner_ms": 1e3 * (perf_counter() - t0), "buses": int(act.sum()), "attained_width_max_pu": float(width.max()),
                            "points_not_certified": int((~np.isfinite(vmin_hi) | ~np.isfinite(vmax_lo)).sum()),
                            "certified_minus_attained_max_pu": float(np.nanmax(np.maximum(vmin_hi - vmin, vmax - vmax_lo)))})
                if enc.verified:
                    cert = vmax_lo - vmin_hi
                    ratio_c, ratio_a = ((enc.hi - enc.lo) / cert)[act], ((enc.hi - enc.lo) / width)[act]
                    row.update({"tube_radius_max_pu": float(enc.s_abs.max()), "rounding_share_of_tube": enc.rounding_share,
                                "width_over_certified_inner_median": float(np.nanmedian(ratio_c)), "width_over_certified_inner_max": float(np.nanmax(ratio_c)),
                                "width_over_attained_median": float(np.median(ratio_a)), "width_over_attained_max": float(ratio_a.max()),
                                "contradictions_with_certified_points": int(((enc.lo > vmin_hi) | (enc.hi < vmax_lo)).sum())})
                    if old.converged:
                        row["tube_radius_over_x8_remainder_bound_median"] = float(np.median(enc.s_abs / np.maximum(old.s_abs[g.pq], 1e-300)))
                    if ref.converged:
                        rr = ((enc.hi - enc.lo) / (ref.hi - ref.lo))[act]
                        row.update({"width_over_x8_width_median": float(np.median(rr)), "width_over_x8_width_min": float(rr.min()),
                                    "width_over_x8_width_max": float(rr.max())})
                    # where the excess comes from, and what the joint bound of first and second order buys
                    ex = lambda b: float(np.median(((b.hi - b.lo) / cert)[act])) - 1.0
                    row.update({"separate_bound_width_over_certified_inner_median": 1.0 + ex(sep), "joint_bound_share_upper": enc.joint_share[1],
                                "joint_bound_share_lower": enc.joint_share[0],
                                "upper_excess_over_width_median_pct": pct(np.nanmedian(((enc.hi - vmax_lo) / cert)[act])),
                                "lower_excess_over_width_median_pct": pct(np.nanmedian(((vmin_hi - enc.lo) / cert)[act])),
                                "tube_over_width_median_pct": pct(np.nanmedian((enc.s_abs / cert)[act]))})
                if kra.closes:
                    cert = vmax_lo - vmin_hi
                    kw = kra.vm_hi - kra.vm_lo
                    row.update({"krawczyk_width_over_certified_inner_median": float(np.nanmedian((kw / cert)[act])),
                                "krawczyk_width_over_certified_inner_max": float(np.nanmax((kw / cert)[act])),
                                "krawczyk_contradictions_with_certified_points": int(((kra.vm_lo > vmin_hi) | (kra.vm_hi < vmax_lo)).sum())})
                    if enc.verified:
                        rk = (kw / (enc.hi - enc.lo))[act]
                        row.update({"krawczyk_width_over_verified_median": float(np.median(rk)), "krawczyk_width_over_verified_min": float(rk.min()),
                                    "krawczyk_width_over_verified_max": float(rk.max()),
                                    "krawczyk_lower_bound_tighter_share": float(np.mean((kra.vm_lo > enc.lo)[act])),
                                    "krawczyk_upper_bound_tighter_share": float(np.mean((kra.vm_hi < enc.hi)[act]))})
                rows.append(row)
                print(f"[x12] {name} {scen} rho_U={rho_U:g}: verified {enc.verified} unique {enc.unique}"
                      + (f" width/certified inner {row['width_over_certified_inner_median']:.4f} (max {row['width_over_certified_inner_max']:.4f};"
                         f" separate bound {row['separate_bound_width_over_certified_inner_median']:.4f})"
                         f" contradictions {row['contradictions_with_certified_points']}" if enc.verified else f" ({enc.reason})")
                      + (f" | Krawczyk {row['krawczyk_width_over_certified_inner_median']:.4f}" if kra.closes else " | Krawczyk does not close")
                      + f" | {1e3 * t_ver:.0f} ms", flush=True)
    return rows


def bracket_row(g, sbus, M, base):
    t0 = perf_counter()
    br = box_bracket(g, sbus, M, BRACKET_TOL, BRACKET_BOXES, NODE_BUDGET, TOL, MIN_WIDTH)
    row = dict(base, K=M.shape[1], leaves=br.leaves, enclosures=br.evaluated, unverified_leaves=br.unverified, seconds=perf_counter() - t0,
               largest_leaf_excess_over_attained_width=br.gap)
    inner, outer = br.vmax_lo - br.vmin_hi, br.vmax_hi - br.vmin_lo
    act = inner > MIN_WIDTH
    row["buses"] = int(act.sum())
    row["points_not_certified"] = int((~np.isfinite(br.vmin_hi) | ~np.isfinite(br.vmax_lo)).sum())
    row["root_verified"] = bool(np.all(np.isfinite(br.root_lo)))
    if br.unverified == 0 and act.any():
        w = (br.root_hi - br.root_lo)
        row.update({"exact_width_max_pu": float(inner[act].max()),
                    "bracket_gap_median_pct": pct(np.median((outer / inner - 1.0)[act])), "bracket_gap_max_pct": pct((outer / inner - 1.0)[act].max()),
                    "contradictions": int(((br.vmin_lo > br.vmin_hi) | (br.vmax_hi < br.vmax_lo)).sum())})
        if row["root_verified"]:
            row.update({"enclosure_over_exact_at_least_median": float(np.median((w / outer)[act])), "enclosure_over_exact_at_most_median": float(np.median((w / inner)[act])),
                        "enclosure_over_exact_at_least_max": float((w / outer)[act].max()), "enclosure_over_exact_at_most_max": float((w / inner)[act].max())})
    return row


def boxes():
    """The boxes with few symbols: (grid, centre injections, injection matrix, labels)."""
    for name, case, prepare, _ in (NETWORKS[0], NETWORKS[4]):
        g = load_grid(case, tol=TOL, prepare=prepare)
        st = make_study(g, zones=spectral_zones(g, 8, laplacian_embedding(g)), tol=TOL)
        for delta in CFG["experiments"]["x8_load_deltas"]:
            yield g, st.sbus, injection_matrix(st, delta), {"set": "load box, 8 zones", "network": name, "scenario": "", "size": delta}
    for name, case, prepare, _ in (NETWORKS[0], NETWORKS[4], NETWORKS[5]):
        g = load_grid(case, tol=TOL, prepare=prepare)
        rating = benchmark_ratings(g)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                st = half_disk_study(g, s0, cap, tol=TOL)
                yield g, st.sbus, injection_matrix(st), {"set": "bounding boxes of the half disks", "network": name, "scenario": scen, "size": rho_U}


def exact_ranges():
    rows = []
    for g, sbus, M, base in boxes():
        rows.append(bracket_row(g, sbus, M, base))
        print(f"[x12] exact {base}: {rows[-1]['enclosures']} enclosures, bracket gap {rows[-1].get('bracket_gap_max_pct', float('nan')):.4f} %", flush=True)
    return rows


def subdivision():
    rows = []
    for g, sbus, M, base in boxes():
        row = dict(base, K=M.shape[1], target_excess=SUBDIVISION_TOL)
        for leaf in ("verified", "krawczyk", "discs"):
            t0 = perf_counter()
            br = box_bracket(g, sbus, M, SUBDIVISION_TOL, SUBDIVISION_BOXES, NODE_BUDGET, TOL, MIN_WIDTH, leaf=leaf)
            inner, outer = br.vmax_lo - br.vmin_hi, br.vmax_hi - br.vmin_lo
            act = inner > MIN_WIDTH
            row.update({f"{leaf}_enclosures": br.evaluated, f"{leaf}_reached_target": bool(br.gap <= SUBDIVISION_TOL), f"{leaf}_largest_excess": br.gap,
                        f"{leaf}_root_ok": bool(np.all(np.isfinite(br.root_lo))), f"{leaf}_unverified_leaves": br.unverified, f"{leaf}_seconds": perf_counter() - t0})
            if np.all(np.isfinite(br.root_lo)):
                row[f"{leaf}_root_over_exact_median"] = float(np.median(((br.root_hi - br.root_lo) / inner)[act]))
            if br.unverified == 0:
                row[f"{leaf}_final_over_exact_median"] = float(np.median((outer / inner)[act]))
        rows.append(row)
        print(f"[x12] subdivision {base}: verified {row['verified_enclosures']} enclosures (target reached {row['verified_reached_target']}),"
              f" Krawczyk {row['krawczyk_enclosures']} ({row['krawczyk_reached_target']}),"
              f" discs {row['discs_enclosures']} (target reached {row['discs_reached_target']}, excess {row['discs_largest_excess']:.3g})", flush=True)
    return rows


PARTS = {"verified": ("x12_verified", verified_cases), "exact": ("x12_exact_range", exact_ranges), "subdivision": ("x12_subdivision", subdivision)}


def main():
    meta = {"script": "experiments/x12_verified.py", "pf_tol": TOL, "rho_U": RHO_U, "samples": N_SAMPLES, "raak_node_budget": NODE_BUDGET,
            "bracket_tol": BRACKET_TOL, "bracket_boxes": BRACKET_BOXES, "min_width_pu": MIN_WIDTH,
            "subdivision_tol": SUBDIVISION_TOL, "subdivision_boxes": SUBDIVISION_BOXES,
            "networks": {**{n: d for n, _, _, d in NETWORKS}, THREE_PHASE: "low-voltage side of the IEEE European LV feeder, merged (X10)"},
            "note": "directed rounding by stepping to the adjacent floating-point number and a priori bounds of sums (cdhim.verified); "
                    "the reference points are verified power-flow solutions (cdhim.certified)"}
    wanted = sys.argv[1:] or list(PARTS)
    unknown = [w for w in wanted if w not in PARTS]
    if unknown:
        sys.exit(f"unknown table(s) {unknown}; choose from {list(PARTS)}")
    for part in wanted:
        table, fn = PARTS[part]
        print(save_table(table, fn(), meta))


if __name__ == "__main__":
    main()
