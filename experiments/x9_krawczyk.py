"""X9 -- The Krawczyk operator as a baseline of the voltage enclosure, on the cases of X8.

Is the enclosure of X8 (second-order affine model + bound of the remainder) better than a well-built interval
method? The baseline is the parametric Krawczyk operator of cdhim.krawczyk: first-order term exact, remainder
by the operator. Eight variants: formulation (power mismatch, current injection) x matrix (derivative, slope)
x coordinates (global axes, axes turned bus by bus to the voltage of the centre).

Both methods get the same information, in two modes:
  box        the bounding boxes of the half disks (what an interval method takes);
  half disk  the first-order term maximised over the half disks themselves, higher orders over the boxes
             (what X8 does; for Krawczyk, ``first_order`` of ``Krawczyk.solve``).

Reference: as in X8, the best values ATTAINED over the half disks (local search, nine boundary load flows,
random samples with the seeds of X8). They are inside the boxes too, so no bound may be on the wrong side of
them in either mode; the ratio "width over attained range" is only meaningful in the half-disk mode.

Writes results/x9_krawczyk.csv              one row per case and variant, with the enclosure of X8 on the same row
       results/x9_krawczyk_load_boxes.csv   boxes of load uncertainty in 8 zones (the cases of x8_load_boxes.csv),
                                            against the exact range over the vertices
       results/x9_krawczyk_subdivision.csv  Krawczyk on the box split into 2^K pieces (networks with few symbols)

    python experiments/x9_krawczyk.py                    # everything
    python experiments/x9_krawczyk.py feeders            # only the named tables (feeders, boxes, subdivision)
"""
import itertools
import sys
from time import perf_counter
from types import SimpleNamespace

import numpy as np

from common import CFG, pct, seed_for
from cdhim import SecondOrder, load_grid, make_study, save_table, spectral_zones
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.krawczyk import Krawczyk
from cdhim.powerflow import BatchACSolver
from cdhim.study import custom_study
from cdhim.validation import VIOLATION_ATOL
from cdhim.zones import laplacian_embedding
from x2_few_sources import vertex_batches
from x8_feeders import NETWORKS, NODE_BUDGET, N_SAMPLES, RHO_U, TOL

VARIANTS = list(itertools.product(("power", "current"), (False, True), (False, True)))      # formulation, slope, rotate
LOAD_DELTAS = CFG["experiments"]["x8_load_deltas"]
LOAD_ZONES = 8
TIGHTER_ATOL = 1e-12                    # p.u.; a bound counts as tighter than another only beyond this
SUBDIVISION_MAX_SYMBOLS = 8             # networks with at most this many symbols get the split box
SUBDIVISION_PIECES = 256                # and as many splits per symbol as keep the number of pieces within this
SUBDIVISION_VARIANT = ("current", True, True)


def variant_fields(formulation, slope, rotate):
    return {"formulation": "power mismatch" if formulation == "power" else "current injection",
            "matrix": "slope" if slope else "derivative", "rotated": bool(rotate)}


def compare(tag, lo, hi, closes, ref_lo, ref_hi, enc_lo, enc_hi, enc_ok, act=None):
    """Statistics of the bounds [lo, hi] of |V| against the reference range and against the enclosure of X8."""
    out = {f"{tag}_closes": bool(closes)}
    if not closes:
        return out
    act = np.ones(len(lo), bool) if act is None else act
    width = (ref_hi - ref_lo)[act]
    ratio = (hi - lo)[act] / width
    out.update({f"{tag}_width_over_reference_median": float(np.median(ratio)), f"{tag}_width_over_reference_max": float(ratio.max()),
                f"{tag}_upper_excess_ave_pct": pct(np.mean(((hi - ref_hi) / ref_hi)[act])),
                f"{tag}_lower_slack_ave_pct": pct(np.mean(((ref_lo - lo) / ref_lo)[act])),
                f"{tag}_violations": int(((hi < ref_hi - VIOLATION_ATOL) | (lo > ref_lo + VIOLATION_ATOL))[act].sum())})
    if enc_ok:
        rel = ((hi - lo) / (enc_hi - enc_lo))[act]
        out.update({f"{tag}_width_over_enclosure_min": float(rel.min()), f"{tag}_width_over_enclosure_median": float(np.median(rel)),
                    f"{tag}_width_over_enclosure_max": float(rel.max()),
                    f"{tag}_lower_bound_tighter_share": float(np.mean((lo > enc_lo + TIGHTER_ATOL)[act])),
                    f"{tag}_upper_bound_tighter_share": float(np.mean((hi < enc_hi - TIGHTER_ATOL)[act]))})
    return out


def attained(g, st, s0, cap, seed):
    """Best attained |V| of every PQ bus over the half disks, exactly as X8 computes them. (vmin, vmax, AC failures)."""
    dg = np.flatnonzero(cap > 0)
    k = np.arange(9)
    S = np.repeat(s0[:, None], 9 + N_SAMPLES, axis=1)
    S[:, :9] += cap[:, None] * (np.sin(k * np.pi / 8) + 1j * np.cos(k * np.pi / 8))[None, :]
    S[dg, 9:] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(seed), len(dg), N_SAMPLES)
    V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, st.V_star, tol=TOL).solve(S)
    vs = np.abs(V[g.pq][:, conv])
    vmax = np.maximum(attained_extreme(g, s0, cap, +1, tol=TOL), vs.max(axis=1))
    vmin = np.minimum(attained_extreme(g, s0, cap, -1, tol=TOL), vs.min(axis=1))
    return vmin, vmax, int((~conv).sum())


def cases():
    """(network, grid, scenario, rho_U, s0, cap, seed) for every case of x8_feeders.csv."""
    for name, case, prepare, _ in NETWORKS:
        g = load_grid(case, tol=TOL, prepare=prepare)
        rating = benchmark_ratings(g)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                yield name, g, scen, rho_U, s0, cap, seed_for("x8", name, scen, rho_U)


def box_mode(so, enc):
    """The enclosure with the first-order term bounded over the box, from the half-disk one (same higher orders)."""
    lo1, hi1 = first_order_range(so.c_lin)
    r = np.abs(so.c_lin).sum(axis=1)
    return enc.lo - lo1 - r, enc.hi - hi1 + r


def feeders():
    rows = []
    for name, g, scen, rho_U, s0, cap, seed in cases():
        base = {"network": name, "scenario": scen, "rho_U": rho_U, "n_bus": g.n_bus, "symbols": 2 * int((cap > 0).sum()), "seed": seed}
        try:
            st = half_disk_study(g, s0, cap, tol=TOL)
        except RuntimeError as e:
            rows.append({**base, "status": str(e)})
            continue
        so = SecondOrder(st, 1.0)
        t0 = perf_counter()
        enc = voltage_enclosure(so, first_order_range(so.c_lin), node_budget=NODE_BUDGET)        # the call of X8
        t_enc = perf_counter() - t0
        enc_box = box_mode(so, enc)
        vmin, vmax, fails = attained(g, st, s0, cap, seed)
        base.update({"status": "ok", "ac_failures_in_samples": fails, "attained_width_max_pu": float((vmax - vmin).max()), "enclosure_ms": 1e3 * t_enc})
        base.update(compare("enclosure_half_disk", enc.lo, enc.hi, enc.converged, vmin, vmax, None, None, False))
        base.update(compare("enclosure_box", *enc_box, enc.converged, vmin, vmax, None, None, False))
        for variant in VARIANTS:
            t0 = perf_counter()
            kr = Krawczyk(so, *variant)
            box = kr.solve()
            hd = kr.solve(first_order_range(kr.L))
            row = {**base, **variant_fields(*variant), "krawczyk_ms": 1e3 * (perf_counter() - t0),
                   "box_operator_evaluations": box.iterations, "box_rel_radius": box.rel_radius, "box_inflation": box.inflation,
                   "half_disk_operator_evaluations": hd.iterations, "half_disk_rel_radius": hd.rel_radius, "half_disk_inflation": hd.inflation}
            row.update(compare("box", box.vm_lo, box.vm_hi, box.closes, vmin, vmax, *enc_box, enc.converged))
            row.update(compare("half_disk", hd.vm_lo, hd.vm_hi, hd.closes, vmin, vmax, enc.lo, enc.hi, enc.converged))
            rows.append(row)
        best = [r for r in rows[-len(VARIANTS):] if r["half_disk_closes"]]
        print(f"[x9] {name} {scen} rho_U={rho_U:g}: enclosure {enc.converged}"
              + (f" {base['enclosure_half_disk_width_over_reference_median']:.4f}" if enc.converged else "")
              + f" | Krawczyk closes in {len(best)} of {len(VARIANTS)} variants"
              + (f", best width/attained {min(r['half_disk_width_over_reference_median'] for r in best):.4f}"
                 f", violations {sum(r['half_disk_violations'] + (r['box_violations'] if r['box_closes'] else 0) for r in best)}" if best else ""), flush=True)
    return rows


def load_boxes():
    """Boxes of load uncertainty in 8 zones (the cases of x8_load_boxes.csv) against the exact range over the vertices."""
    rows = []
    for name, case, prepare, _ in (NETWORKS[0], NETWORKS[4]):
        g = load_grid(case, tol=TOL, prepare=prepare)
        st = make_study(g, zones=spectral_zones(g, LOAD_ZONES, laplacian_embedding(g)), tol=TOL)
        for delta in LOAD_DELTAS:
            so = SecondOrder(st, delta)
            enc = voltage_enclosure(so, node_budget=NODE_BUDGET)
            vmin, vmax = np.full(len(g.pq), np.inf), np.full(len(g.pq), -np.inf)
            fails = 0
            for eps in vertex_batches(st.src.n):
                V, conv, _ = st.solver.solve(st.injections(eps, delta))
                fails += int((~conv).sum())
                vm = np.abs(V[g.pq][:, conv])
                vmin, vmax = np.minimum(vmin, vm.min(axis=1)), np.maximum(vmax, vm.max(axis=1))
            act = np.abs(so.c_lin).sum(axis=1) > 1e-9
            base = {"network": name, "K": st.src.n, "delta": delta, "vertices_without_ac_solution": fails,
                    "vertex_range_width_max_pu": float((vmax - vmin)[act].max())}
            base.update(compare("enclosure", enc.lo, enc.hi, enc.converged, vmin, vmax, None, None, False, act))
            for variant in VARIANTS:
                t0 = perf_counter()
                res = Krawczyk(so, *variant).solve()
                row = {**base, **variant_fields(*variant), "krawczyk_ms": 1e3 * (perf_counter() - t0),
                       "krawczyk_operator_evaluations": res.iterations, "krawczyk_rel_radius": res.rel_radius}
                row.update(compare("krawczyk", res.vm_lo, res.vm_hi, res.closes, vmin, vmax, enc.lo, enc.hi, enc.converged, act))
                rows.append(row)
            best = [r for r in rows[-len(VARIANTS):] if r["krawczyk_closes"]]
            print(f"[x9] load boxes {name} delta={delta:g}: enclosure {enc.converged}"
                  + (f" {base['enclosure_width_over_reference_median']:.4f}" if enc.converged else "")
                  + (f" | best Krawczyk {min(r['krawczyk_width_over_reference_median'] for r in best):.4f}" if best else " | Krawczyk does not close"), flush=True)
    return rows


def subdivision():
    """Krawczyk on the box of the symbols split into equal pieces: the union of the enclosures of the pieces.

    Its excess is of second order in the size of a piece, so splitting tightens it, at a cost that grows as
    (splits per symbol)^K. Box mode only, against the enclosure of X8 in box mode and the attained values.
    """
    rows = []
    for name, g, scen, rho_U, s0, cap, seed in cases():
        K = 2 * int((cap > 0).sum())
        if K > SUBDIVISION_MAX_SYMBOLS:
            continue
        try:
            st = half_disk_study(g, s0, cap, tol=TOL)
        except RuntimeError as e:
            rows.append({"network": name, "scenario": scen, "rho_U": rho_U, "symbols": K, "status": str(e)})
            continue
        so = SecondOrder(st, 1.0)
        t0 = perf_counter()
        enc = voltage_enclosure(so, first_order_range(so.c_lin), node_budget=NODE_BUDGET)
        t_enc = perf_counter() - t0
        enc_box = box_mode(so, enc)
        vmin, vmax, _ = attained(g, st, s0, cap, seed)
        src = st.src
        splits = 1
        while (2 * splits) ** K <= SUBDIVISION_PIECES:
            splits *= 2
            t0 = perf_counter()
            lo, hi, ok, why = np.full(len(g.pq), np.inf), np.full(len(g.pq), -np.inf), True, ""
            centres = (2.0 * np.arange(splits) + 1.0) / splits - 1.0
            for c in itertools.product(centres, repeat=K):
                c = np.array(c)
                centre = st.injections(c[:, None], 1.0)[:, 0]
                try:
                    sub = custom_study(g, centre, src.bus, src.is_q, src.unit / splits, tol=TOL)
                except RuntimeError as e:
                    ok, why = False, str(e)
                    break
                res = Krawczyk(SimpleNamespace(st=sub, delta=1.0), *SUBDIVISION_VARIANT).solve()
                if not res.closes:
                    ok, why = False, "the operator does not close on a piece"
                    break
                lo, hi = np.minimum(lo, res.vm_lo), np.maximum(hi, res.vm_hi)
            row = {"network": name, "scenario": scen, "rho_U": rho_U, "symbols": K, "status": "ok", "splits_per_symbol": splits,
                   "pieces": splits ** K, **variant_fields(*SUBDIVISION_VARIANT), "krawczyk_ms": 1e3 * (perf_counter() - t0), "piece_failure": why,
                   "enclosure_ms": 1e3 * t_enc}
            row.update(compare("enclosure_box", *enc_box, enc.converged, vmin, vmax, None, None, False))
            row.update(compare("box", lo, hi, ok, vmin, vmax, *enc_box, enc.converged))
            rows.append(row)
            print(f"[x9] subdivision {name} {scen} rho_U={rho_U:g} K={K}: {splits ** K} pieces, closes {ok}"
                  + (f", width over enclosure (box mode) {row['box_width_over_enclosure_median']:.4f}" if ok and enc.converged else ""), flush=True)
    return rows


PARTS = {"feeders": ("x9_krawczyk", feeders), "boxes": ("x9_krawczyk_load_boxes", load_boxes), "subdivision": ("x9_krawczyk_subdivision", subdivision)}


def main():
    meta = {"raak_node_budget": NODE_BUDGET, "samples_per_case": N_SAMPLES, "rho_U": RHO_U, "load_deltas": LOAD_DELTAS, "load_zones": LOAD_ZONES,
            "seed_rule": "common.seed_for('x8', network, scenario, rho_U): the reference values are those of X8; stored per row",
            "variants": [variant_fields(*v) for v in VARIANTS], "tighter_atol_pu": TIGHTER_ATOL,
            "subdivision": {"max_symbols": SUBDIVISION_MAX_SYMBOLS, "max_pieces": SUBDIVISION_PIECES, "variant": variant_fields(*SUBDIVISION_VARIANT)},
            "networks": {n: d for n, _, _, d in NETWORKS},
            "note": "bounds valid in exact arithmetic, no outward rounding; reference = best attained values (exact vertex range in the load boxes); "
                    "times are indicative, other jobs were running"}
    wanted = sys.argv[1:] or list(PARTS)
    unknown = [w for w in wanted if w not in PARTS]
    if unknown:
        sys.exit(f"unknown table(s) {unknown}; choose from {list(PARTS)}")
    for part in wanted:
        table, fn = PARTS[part]
        print(save_table(table, fn(), meta))


if __name__ == "__main__":
    main()
