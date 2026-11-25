"""X5 -- The closest published method (Liao et al., IET GTD 2017) against the direct |V| form.

Liao et al. expand the state x = (e, f) to second order in the perturbation variables,
    x = x_c + (dx/dxi) dxi + 1/2 (d2x/dxi2) [dxi_i dxi_j],                       (their Eq. 26)
with -J dx/dxi_j = dF/dxi_j (Eq. 23) and -J d2x/dxi_i dxi_j = (dx/dxi_i)' H(F) (dx/dxi_j) (Eq. 25), and
evaluate it with interval arithmetic. Reimplemented here as: the quadratic form of every state component
bounded monomial by monomial, with [dxi]^2 = [0, 1] (their Eq. 14). The paper does not say how the
intervals of e and f become an interval of |V|; interval arithmetic on sqrt([e]^2 + [f]^2) is used.

Compared with the direct form of this package, d|V_i| = c_i' eps + eps' H_i eps (the same second-order
information, projected on |V_i| before bounding), bounded monomial by monomial (AF2) and with RAA-K.

1. IEEE 9, +-10 % on the active power of the PV generators and on P and Q of the loads: compared with
   Table 1 of the paper, with all 2^K vertices and with Monte Carlo.
2. IEEE 57, +-20 %, same variables.
3. PEGASE-2869 with 8 zone symbols (both balancing models), reference from results/x2_bus.csv.

The paper counts 12 perturbation variables on IEEE 9 without listing them; five readings of the uncertain
set are tried against its Table 1 (results/x5_liao_variants.csv).

Writes results/x5_liao_table1.csv, results/x5_liao_variants.csv and results/x5_liao_summary.csv
"""
import numpy as np
import pandas as pd

from common import CFG, ROOT, grid, method, pct, seed_for, study
from cdhim import SecondOrder, make_study, participation_factors, sample_eps, save_table, spectral_zones
from cdhim.core import run_cdhim
from cdhim.enclosure import state_forms
from cdhim.knomial import knomial_bound
from cdhim.powerflow import solve_blocked
from cdhim.study import Sources, voltage_coefficients
from cdhim.validation import VIOLATION_ATOL
from cdhim.zones import laplacian_embedding
from x2_few_sources import BATCH, R_FLOOR, model_minimiser_vertices, vertex_batches

N_MC = CFG["experiments"]["n_ac"]
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]

# Table 1 of Liao et al. (2017), IEEE 9-bus, +-10 %: bus -> (nominal, proposed lower, upper, Monte Carlo lower, upper)
LIAO_TABLE1 = {4: (0.9870, 0.9717, 0.9923, 0.9793, 0.9910), 5: (0.9755, 0.9672, 0.9839, 0.9711, 0.9805),
               6: (1.0034, 1.0001, 1.0066, 1.0010, 1.0062), 7: (0.9856, 0.9793, 0.9926, 0.9823, 0.9896),
               8: (0.9962, 0.9917, 1.0088, 0.9931, 1.0058), 9: (0.9576, 0.9470, 0.9682, 0.9485, 0.9649)}


def injection_study(case, gen_p=True, load_p=True, load_q=True):
    """Study whose symbols are P of the PV generators and P, Q of the loads (the variables of Liao et al.)."""
    g = grid(case)
    st = make_study(g, tol=CFG["powerflow"]["tol_pu"])
    n_th = len(g.pvpq)
    row_p, row_q = np.full(g.n_bus, -1), np.full(g.n_bus, -1)
    row_p[g.pvpq] = np.arange(n_th)
    row_q[g.pq] = n_th + np.arange(len(g.pq))
    pg = g.s_gen.real
    gen = g.pv[np.abs(pg[g.pv]) > 1e-6] if gen_p else g.pv[:0]
    pd_, qd = g.s_load.real, g.s_load.imag
    bp = np.flatnonzero((np.abs(pd_) > 1e-6) & (row_p >= 0) & load_p)
    bq = np.flatnonzero((np.abs(qd) > 1e-6) & (row_q >= 0) & load_q)
    bus = np.concatenate([gen, bp, bq])
    is_q = np.concatenate([np.zeros(len(gen) + len(bp), bool), np.ones(len(bq), bool)])
    unit = np.concatenate([pg[gen], -pd_[bp], -qd[bq]])
    st.src = Sources(bus, is_q, unit, np.where(is_q, row_q[bus], row_p[bus]))
    st.C1 = voltage_coefficients(st.lin, g, st.src)
    return g, st


def monomial_all(H):
    d = np.einsum("ikk->ik", H)
    off = np.abs(H).sum(axis=(1, 2)) - np.abs(d).sum(axis=1)
    return -np.maximum(-d, 0.0).sum(axis=1) - off, np.maximum(d, 0.0).sum(axis=1) + off


def sq(lo, hi):
    """Interval square (Eq. 14 of the paper)."""
    a, b = lo * lo, hi * hi
    return np.where((lo <= 0) & (hi >= 0), 0.0, np.minimum(a, b)), np.maximum(a, b)


def liao_state_intervals(so):
    """Eq. (26) of Liao et al.: interval of every state component (e, f of the non-slack buses)."""
    g = so.st.grid
    Hs = state_forms(so)
    lo2, hi2 = monomial_all(Hs)
    rad = np.abs(so.C).sum(axis=1)
    x0 = np.concatenate([so.st.V_star[g.pvpq].real, so.st.V_star[g.pvpq].imag])
    return x0 - rad + lo2, x0 + rad + hi2, Hs


def liao_intervals(so):
    """|V| interval of every PQ bus: intervals of e and f, then interval arithmetic on sqrt(e^2 + f^2)."""
    g, m = so.st.grid, so.m
    lo, hi, Hs = liao_state_intervals(so)
    ie, jf = so.col[g.pq], m + so.col[g.pq]
    # the same second-order information projected on |V| must be the adjoint form of this package
    er, ei = (so.st.V_star[g.pq].real / so.vm)[:, None, None], (so.st.V_star[g.pq].imag / so.vm)[:, None, None]
    proj = er * Hs[ie] + ei * Hs[jf] + so.c_perp[:, :, None] * so.c_perp[:, None, :] / (2.0 * so.vm)[:, None, None]
    H = so.forms()
    assert np.allclose(proj, H, rtol=1e-6, atol=1e-9 * np.abs(H).max()), "forward (Liao) and adjoint second-order forms differ"
    e2, f2 = sq(lo[ie], hi[ie]), sq(lo[jf], hi[jf])
    return np.sqrt(e2[0] + f2[0]), np.sqrt(e2[1] + f2[1])


def direct_intervals(so, C):
    H = so.forms()
    r = np.abs(C).sum(axis=1)
    out = {"af2": monomial_all(H)}
    vals = np.array([knomial_bound(h, H.shape[1], NODE_BUDGET)[:2] for h in H])
    out["raak"] = (vals[:, 0], vals[:, 1])
    v = so.vm
    return {name: (v - r + lo, v + r + hi) for name, (lo, hi) in out.items()}, H, r


def reference(st, C, H, delta, key, enumerate_all):
    """(lower, upper, description, MC 0.135 % and 99.865 % quantiles) of the exact AC |V| of every PQ bus."""
    g = st.grid
    K = st.src.n
    lo, hi = np.full(len(g.pq), np.inf), np.full(len(g.pq), -np.inf)
    x_lo, x_hi = np.full(2 * len(g.pvpq), np.inf), np.full(2 * len(g.pvpq), -np.inf)
    if enumerate_all:
        batches, what = vertex_batches(K), f"all 2^{K} vertices"
    else:
        first = -np.sign(C)
        first[first == 0] = 1.0
        cand = np.unique(np.vstack([first, -first, model_minimiser_vertices(C, H), model_minimiser_vertices(-C, -H)]), axis=0).T
        batches, what = (cand[:, s:s + BATCH] for s in range(0, cand.shape[1], BATCH)), "extreme and model-optimal vertices of every bus"
    n_fail = 0
    for eps in batches:
        V, conv, _ = st.solver.solve(st.injections(np.ascontiguousarray(eps), delta))
        n_fail += int((~conv).sum())
        vm = np.abs(V[g.pq][:, conv])
        lo, hi = np.minimum(lo, vm.min(axis=1)), np.maximum(hi, vm.max(axis=1))
        x = np.vstack([V[g.pvpq][:, conv].real, V[g.pvpq][:, conv].imag])
        x_lo, x_hi = np.minimum(x_lo, x.min(axis=1)), np.maximum(x_hi, x.max(axis=1))
    seed = seed_for("x5", key)
    rng = np.random.default_rng(seed)
    samples = []
    for start in range(0, N_MC, BATCH):
        V, conv, _ = st.solver.solve(st.injections(sample_eps(rng, K, min(BATCH, N_MC - start)), delta))
        samples.append(np.abs(V[g.pq][:, conv]).astype(np.float32))
    vm = np.concatenate(samples, axis=1)
    return lo, hi, what, vm.min(axis=1), vm.max(axis=1), seed, n_fail, x_lo, x_hi


def summarise(key, label, st, so, C, delta, ref, extra):
    g = st.grid
    lo_r, hi_r, what, mc_lo, mc_hi, seed, n_fail, x_lo, x_hi = ref
    direct, H, r = direct_intervals(so, C)
    intervals = {"liao": liao_intervals(so), **direct}
    act = r > R_FLOOR
    K = st.src.n
    row = {"case": label, "system": key, "K": K, "delta": delta, "n_pq": int(act.sum()), "reference": what, "seed": seed,
           "reference_solves_without_solution": n_fail,
           "angle_radius_median_rad": float(np.median((np.abs(so.c_perp).sum(axis=1) / so.vm)[act])),
           "linear_solves_liao": K + K * (K + 1) // 2, "linear_solves_direct": K + int(act.sum()), **extra}
    for name, (lo, hi) in intervals.items():
        w = (hi - lo)[act]
        row[f"{name}_width_over_reference_median"] = float(np.median(w / (hi_r - lo_r)[act]))
        row[f"{name}_width_over_reference_max"] = float((w / (hi_r - lo_r)[act]).max())
        row[f"{name}_width_over_mc_median"] = float(np.median(w / (mc_hi - mc_lo)[act]))
        short = np.maximum(lo - lo_r, hi_r - hi)[act]
        row[f"{name}_buses_not_containing_reference"] = int((short > VIOLATION_ATOL).sum())
        row[f"{name}_shortfall_max_pu"] = float(max(short.max(), 0.0))
    if x_lo is not None:                                    # Eq. (26) itself, on the state components it is written for
        s_lo, s_hi, _ = liao_state_intervals(so)
        span = x_hi - x_lo
        ok = span > R_FLOOR
        row["liao_state_width_over_reference_median"] = float(np.median(((s_hi - s_lo) / span)[ok]))
        row["liao_state_width_over_reference_max"] = float(((s_hi - s_lo) / span)[ok].max())
        short = np.maximum(s_lo - x_lo, x_hi - s_hi)[ok]
        row["liao_state_components_not_containing_reference"] = int((short > VIOLATION_ATOL).sum())
        row["liao_state_components"] = int(ok.sum())
        row["liao_state_shortfall_max_pu"] = float(max(short.max(), 0.0))
    row["liao_over_af2_width_median"] = float(np.median(((intervals["liao"][1] - intervals["liao"][0]) / (intervals["af2"][1] - intervals["af2"][0]))[act]))
    row["liao_over_af2_width_max"] = float(((intervals["liao"][1] - intervals["liao"][0]) / (intervals["af2"][1] - intervals["af2"][0]))[act].max())
    print(f"[x5] {label}: Liao / reference width median {row['liao_width_over_reference_median']:.3f}, direct AF2 {row['af2_width_over_reference_median']:.3f}, "
          f"RAA-K {row['raak_width_over_reference_median']:.3f}", flush=True)
    return row, intervals, (lo_r, hi_r, mc_lo, mc_hi)


def main():
    rows, table1, variants = [], [], []
    for case, delta, enum in (("case9", 0.10, True), ("case57", 0.20, False)):
        g, st = injection_study(case)
        so = SecondOrder(st, delta)
        C = st.coefficients(delta)
        assert np.allclose(so.c_lin, C, rtol=1e-6, atol=1e-10)
        ref = reference(st, C, so.forms(), delta, case, enum)
        row, iv, (lo_r, hi_r, mc_lo, mc_hi) = summarise(case, f"{case}, +-{100 * delta:.0f} % on generator P and load P, Q", st, so, C, delta, ref,
                                                         {"balancing": "slack"})
        rows.append(row)
        if case == "case9":
            paper = np.array([LIAO_TABLE1[int(g.bus_label[b]) + 1] for b in g.pq])
            for name, flags in (("generator P, load P and Q", (True, True, True)), ("load P and Q", (False, True, True)),
                                ("generator P, load P", (True, True, False)), ("load P", (False, True, False)), ("generator P", (True, False, False))):
                gv, sv = injection_study(case, *flags)
                sov = SecondOrder(sv, delta)
                Cv = sv.coefficients(delta)
                rv = reference(sv, Cv, sov.forms(), delta, case, True)
                lv = liao_intervals(sov)
                dv, _, _ = direct_intervals(sov, Cv)
                variants.append({"uncertain_variables": name, "K": sv.src.n,
                                 "max_abs_diff_exact_range_vs_paper_mc": float(max(np.abs(rv[0] - paper[:, 3]).max(), np.abs(rv[1] - paper[:, 4]).max())),
                                 "max_abs_diff_mc_here_vs_paper_mc": float(max(np.abs(rv[3] - paper[:, 3]).max(), np.abs(rv[4] - paper[:, 4]).max())),
                                 "max_abs_diff_liao_here_vs_paper": float(max(np.abs(lv[0] - paper[:, 1]).max(), np.abs(lv[1] - paper[:, 2]).max())),
                                 "max_abs_diff_direct_af2_vs_paper": float(max(np.abs(dv["af2"][0] - paper[:, 1]).max(), np.abs(dv["af2"][1] - paper[:, 2]).max())),
                                 "paper_interval_contains_exact_range_buses": int(((paper[:, 1] <= rv[0] + 1e-4) & (paper[:, 2] >= rv[1] - 1e-4)).sum())})
            for k, b in enumerate(g.pq):
                label = int(g.bus_label[b]) + 1                                  # pandapower index 0..8 -> bus 1..9 of the paper
                nom, plo, phi, pmlo, pmhi = LIAO_TABLE1[label]
                table1.append({"bus": label, "nominal_paper": nom, "nominal_here": float(np.abs(st.V_star[b])),
                               "liao_paper_lower": plo, "liao_paper_upper": phi, "liao_here_lower": iv["liao"][0][k], "liao_here_upper": iv["liao"][1][k],
                               "direct_af2_lower": iv["af2"][0][k], "direct_af2_upper": iv["af2"][1][k],
                               "direct_raak_lower": iv["raak"][0][k], "direct_raak_upper": iv["raak"][1][k],
                               "exact_all_vertices_lower": lo_r[k], "exact_all_vertices_upper": hi_r[k],
                               "mc_paper_lower": pmlo, "mc_paper_upper": pmhi, "mc_here_lower": float(mc_lo[k]), "mc_here_upper": float(mc_hi[k])})
    bus = pd.read_csv(ROOT / "results" / "x2_bus.csv")
    key = "pegase2869"
    sysc = CFG["systems"][key]
    g = grid(sysc["case"])
    emb = laplacian_embedding(g)
    for balancing in CFG["experiments"]["balancing"][key]:
        st = study(sysc["case"], zones=spectral_zones(g, 8, emb), participation=participation_factors(g) if balancing == "participation" else None)
        delta = sysc["delta"]
        so = SecondOrder(st, delta)
        C = st.coefficients(delta)
        b = bus[(bus.system == key) & (bus.balancing == balancing) & (bus.K == 8) & np.isclose(bus.delta, delta)].set_index("bus")
        res = run_cdhim(C, st.fi, st.ti, method())
        lab = [int(g.bus_label[i]) for i in g.pq]
        lo_r, hi_r = np.full(len(g.pq), np.nan), np.full(len(g.pq), np.nan)
        mc_lo, mc_hi = lo_r.copy(), hi_r.copy()
        act = np.flatnonzero(res.r_affine > R_FLOOR)
        sel = b.loc[[lab[i] for i in act]]
        lo_r[act], hi_r[act] = so.vm[act] + sel.vertex_dv_min.to_numpy(), so.vm[act] + sel.vertex_dv_max.to_numpy()
        mc_lo[act], mc_hi[act] = so.vm[act] + sel.mc_q00135.to_numpy(), so.vm[act] + sel.mc_q99865.to_numpy()
        ref = (lo_r, hi_r, "all 2^8 vertices (results/x2_bus.csv); Monte Carlo column is the 99.73 % range", mc_lo, mc_hi, "", 0, None, None)
        row, _, _ = summarise(key, f"PEGASE-2869, 8 zones, {balancing}", st, so, C, delta, ref, {"balancing": balancing})
        rows.append(row)
    meta = {"n_mc": N_MC, "raak_node_budget": NODE_BUDGET,
            "note": "Liao et al. reimplemented from the paper; the step from [e], [f] to [|V|] is not specified there (interval arithmetic used)"}
    print(save_table("x5_liao_table1", table1, meta))
    print(save_table("x5_liao_variants", variants, meta))
    print(save_table("x5_liao_summary", rows, meta))


if __name__ == "__main__":
    main()
