"""X4 -- Critical checks of the reference used in X2 and X3.

A. Independent solver. At the vertex where the second-order bound is most violated, the same scenario is
   built in pandapower by scaling the load and static-generator elements zone by zone (and moving the
   generator set-points in the participation model) and solved with its Newton-Raphson. Compared with
   the batched solver of this package and with the bounds.

B. Is the extreme of the box at a vertex? X2 took the range of the exact AC solution over the 2^K
   vertices as the reference. Here the second-order model is minimised / maximised over the whole box
   (coordinate descent, interior points allowed) for every bus, and the exact AC power flow is solved
   at that point. Any bus where it beats the best vertex shows that the vertex range is too narrow.

C. Polishing with the exact model. For the buses where B gains most, a bounded quasi-Newton search
   (L-BFGS-B, finite-difference gradient) on the exact AC voltage, started from the best point known.

Needs results/x2_bus.csv. Deterministic. Writes results/x4_independent.csv and results/x4_interior.csv
"""
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from common import CFG, ROOT, grid, method, pct, study
from cdhim import SecondOrder, participation_factors, save_table, spectral_zones
from cdhim.knomial import knomial_bound
from cdhim.network import DATA_DIR, _build_net
from cdhim.powerflow import mismatch, stack_mismatch
from cdhim.second_order import bound_monomial
from cdhim.zones import laplacian_embedding
from x2_few_sources import BATCH, K_ENUM, R_FLOOR, model_minimiser_vertices, vertex_batches
from x3_worst_vertex import ac_at_own_vertex

N_POLISH = 10          # buses per case and side sent to the exact-model search
GAIN_TOL = 1e-9        # p.u.; an interior point "beats" the vertices only above the accuracy of the AC solver


def build(key, balancing, k, embedding):
    sysc = CFG["systems"][key]
    g = grid(sysc["case"])
    zones = spectral_zones(g, k, embedding)
    alpha = participation_factors(g) if balancing == "participation" else None
    st = study(sysc["case"], zones=zones, participation=alpha)
    delta = sysc["delta"]
    so = SecondOrder(st, delta)
    return g, st, so, delta, st.run(delta, method())


# ----------------------------------------------------------------------------- A. pandapower
def pandapower_solution(g, st, eps, delta):
    """|V| of every internal bus from pandapower's own Newton-Raphson for the scenario ``eps``."""
    import pandapower as pp

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net, _ = _build_net(g.name, DATA_DIR)
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, max_iteration=50)
    lookup = np.asarray(net._pd2ppc_lookups["bus"])
    src = st.src
    rel = delta * src.expand(eps[:, None])[:, 0]                       # relative change of every uncertain injection
    fp, fq = np.ones(g.n_bus), np.ones(g.n_bus)
    fp[src.bus[~src.is_q]] += rel[~src.is_q]
    fq[src.bus[src.is_q]] += rel[src.is_q]
    for table in ("load", "sgen"):
        el = getattr(net, table)
        if len(el):
            ib = lookup[el.bus.values]
            el["p_mw"] = el.p_mw.values * fp[ib]
            el["q_mvar"] = el.q_mvar.values * fq[ib]
    if src.G is not None:                                              # generator response of the participation model
        dp = np.asarray(delta * (src.G @ eps)).ravel() * g.base_mva
        ib = lookup[net.gen.bus.values]
        first = np.unique(ib, return_index=True)[1]                    # one generator per bus takes the change
        net.gen.loc[net.gen.index[first], "p_mw"] = net.gen.p_mw.values[first] + dp[ib[first]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, max_iteration=50, tolerance_mva=1e-9)
    idx = np.searchsorted(net.bus.index.values, g.bus_label)
    return net.res_bus.vm_pu.values[idx]


def independent_check(key, balancing, k, embedding):
    g, st, so, delta, res = build(key, balancing, k, embedding)
    C, r, H = st.coefficients(delta), res.r_affine, so.forms()
    active = np.flatnonzero(r > R_FLOOR)
    lam, vec = np.linalg.eigh(H)
    s2 = np.abs(vec).sum(axis=1) ** 2
    lo = {"af2": np.array([bound_monomial(h)[0] for h in H]), "raak": np.array([knomial_bound(h, k, CFG["experiments"]["raak_node_budget"])[0] for h in H]),
          "spectral": -(np.maximum(-lam, 0.0) * s2).sum(axis=1)}
    worst = (-np.inf, None, None)
    for eps in vertex_batches(st.src.n):                                # vertex and bus with the largest shortfall of the tightest bound
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        dv = np.abs(V[g.pq]) - st.vm_star[:, None]
        short = (-r + lo["af2"])[:, None] - dv
        short[r <= R_FLOOR] = -np.inf
        i, j = np.unravel_index(np.argmax(short), short.shape)
        if short[i, j] > worst[0]:
            worst = (short[i, j], i, eps[:, j].copy())
    _, i, eps = worst
    S = st.injections(eps[:, None], delta)
    V, conv, _ = st.solver.solve(S)
    resid = np.abs(stack_mismatch(mismatch(g.Ybus, V, S), g.pvpq, g.pq)).max()
    vm_batch = np.abs(V[:, 0])
    vm_pp = pandapower_solution(g, st, eps, delta)
    dv_b, dv_p = vm_batch[g.pq[i]] - st.vm_star[i], vm_pp[g.pq[i]] - st.vm_star[i]
    row = {"system": key, "balancing": balancing, "K": st.src.n, "delta": delta, "bus": int(g.bus_label[g.pq[i]]),
           "batch_solver_max_mismatch_pu": float(resid), "max_abs_vm_difference_pandapower_pu": float(np.abs(vm_pp - vm_batch).max()),
           "r_affine": r[i], "dv_ac_batch": dv_b, "dv_ac_pandapower": dv_p, "affine_lower_edge": -r[i],
           "affine_box_violated_by_pu_pandapower": -r[i] - dv_p}
    for name, v in lo.items():
        row[f"{name}_lower_bound"] = -r[i] + v[i]
        row[f"{name}_violated_by_pu_batch"] = -r[i] + v[i] - dv_b
        row[f"{name}_violated_by_pu_pandapower"] = -r[i] + v[i] - dv_p
    print(f"[x4] A {key} {balancing}: bus {row['bus']}, |vm_pp - vm_batch| max {row['max_abs_vm_difference_pandapower_pu']:.2e}", flush=True)
    return row


# ----------------------------------------------------------------------------- B, C. interior extremes
def box_minimiser(C, H, X, sweeps=300):
    """Coordinate descent of c'x + x'Hx on [-1, 1]^K for every bus at once (exact 1-D minimiser, interior allowed)."""
    X = X.copy()
    d = np.einsum("ikk->ik", H)
    for _ in range(sweeps):
        moved = 0.0
        for k in range(C.shape[1]):
            a = d[:, k]
            b = C[:, k] + 2.0 * (np.einsum("il,il->i", H[:, k, :], X) - a * X[:, k])
            best = np.where(a - b <= a + b, -1.0, 1.0)                 # value at -1 is a - b, at +1 is a + b
            inner = np.clip(-b / (2.0 * np.where(a > 0, a, 1.0)), -1.0, 1.0)
            take = (a > 0) & (a * inner * inner + b * inner < a + b * best)
            best = np.where(take, inner, best)
            moved = max(moved, float(np.abs(best - X[:, k]).max()))
            X[:, k] = best
        if moved < 1e-13:
            break
    return X


def polish(st, i, x0, delta, sign):
    """Bounded quasi-Newton search of sign * d|V_i| on the exact AC model. Returns (x, value of d|V_i|)."""
    k = len(x0)
    h = 1e-4

    def fun(x):
        step = np.where(x + h <= 1.0, h, -h)
        pts = np.column_stack([x] + [x + step[j] * np.eye(k)[j] for j in range(k)])
        V, conv, _ = st.solver.solve(st.injections(pts, delta))
        assert conv.all()
        dv = np.abs(V[st.grid.pq[i]]) - st.vm_star[i]
        return sign * dv[0], sign * (dv[1:] - dv[0]) / step

    out = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=[(-1.0, 1.0)] * k, options={"maxiter": 60, "ftol": 1e-15, "gtol": 1e-12})
    V, conv, _ = st.solver.solve(st.injections(out.x[:, None], delta))
    return out.x, float(np.abs(V[st.grid.pq[i], 0]) - st.vm_star[i])


def interior_check(key, balancing, k, embedding, bus_ref):
    g, st, so, delta, res = build(key, balancing, k, embedding)
    C, r, H = st.coefficients(delta), res.r_affine, so.forms()
    active = np.flatnonzero(r > R_FLOOR)
    ref = bus_ref.set_index("bus").loc[[int(g.bus_label[g.pq[i]]) for i in active]]
    rows = []
    for side, sgn, true in (("lowest", 1.0, ref.vertex_dv_min.to_numpy()), ("highest", -1.0, ref.vertex_dv_max.to_numpy())):
        Xv = model_minimiser_vertices(sgn * C, sgn * H)                # best vertex of the second-order model
        Xb = box_minimiser(sgn * C, sgn * H, Xv)                       # its optimum over the whole box
        interior = (np.abs(Xb) < 1.0 - 1e-9).any(axis=1)[active]
        dv, n_solves = ac_at_own_vertex(st, Xb, active, delta)
        gain = sgn * (true - dv)                                       # > 0: the interior candidate beats every vertex
        order = np.argsort(-gain)[:N_POLISH]
        pol = np.zeros(len(order))
        for t, j in enumerate(order):                                  # exact-model search from the best point known
            _, v = polish(st, active[j], Xb[active[j]], delta, sgn)
            pol[t] = sgn * (true[j] - v)
        rr = r[active]
        rows.append({"system": key, "balancing": balancing, "K": st.src.n, "delta": delta, "side": side, "n_buses": len(active),
                     "buses_model_optimum_interior": int(interior.sum()),
                     "buses_beating_all_vertices": int((gain > GAIN_TOL).sum()),
                     "buses_beating_all_vertices_pct": pct((gain > GAIN_TOL).mean()),
                     "gain_over_vertices_max_pu": float(gain.max()), "gain_over_vertices_over_r_max_pct": pct((gain / rr).max()),
                     "buses_gain_over_1e-5pu": int((gain > 1e-5).sum()), "buses_gain_over_1e-4pu": int((gain > 1e-4).sum()),
                     "vertex_range_too_narrow_median_pct": pct(np.median(np.maximum(gain, 0.0) / (ref.vertex_dv_max - ref.vertex_dv_min).to_numpy())),
                     "ac_solves": n_solves, "polished_buses": len(order),
                     "polish_gain_max_pu": float(pol.max()), "polish_extra_over_model_point_max_pu": float((pol - gain[order]).max())})
        print(f"[x4] B {key} {balancing} K={st.src.n} {side}: {rows[-1]['buses_beating_all_vertices']} buses beat the vertices, "
              f"max {rows[-1]['gain_over_vertices_max_pu']:.2e} pu; after polishing {rows[-1]['polish_gain_max_pu']:.2e}", flush=True)
    return rows


def main():
    bus = pd.read_csv(ROOT / "results" / "x2_bus.csv")
    ind, inter = [], []
    for key, counts in CFG["experiments"]["zone_counts"].items():
        emb = laplacian_embedding(grid(CFG["systems"][key]["case"]))
        delta = CFG["systems"][key]["delta"]
        for balancing in CFG["experiments"]["balancing"][key]:
            ind.append(independent_check(key, balancing, 8, emb))
            for k in [c for c in counts if c <= K_ENUM]:
                ref = bus[(bus.system == key) & (bus.balancing == balancing) & (bus.K == k) & np.isclose(bus.delta, delta)]
                inter += interior_check(key, balancing, k, emb, ref)
    meta = {"note": "deterministic; reference = exact AC at all 2^K vertices (results/x2_bus.csv)", "gain_tolerance_pu": GAIN_TOL}
    print(save_table("x4_independent", ind, meta))
    print(save_table("x4_interior", inter, meta))


if __name__ == "__main__":
    main()
