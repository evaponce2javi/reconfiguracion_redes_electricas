"""X2 -- Gate 2 of CHECKLIST_AA.md: does the second-order affine power flow close with few sources?

Uncertainty model: K zone symbols. Every load of a zone scales by (1 + delta * eps_k), P and Q together.
Two balancing models on PEGASE: "slack" (the slack bus absorbs the whole imbalance, as in the manuscript)
and "participation" (generation follows the demand change in fixed shares, see generator_response).

For every case (system, balancing, K, delta):

1. first- and second-order model of d|V| for all PQ buses; H_i is K x K and is assembled for every bus;
2. enclosures of eps' H_i eps for every bus: standard affine product term by term, AF2 monomial-wise,
   the group's RAA-K rewriting (cdhim.knomial, ported from ibex-lib) with cliques of up to 2, 3, 6 and K
   variables, and spectral blocks (1 / 8 / all) bounded with the same k-nomial theorem;
3. reference: exact AC power flow at ALL 2^K vertices of the noise box when K <= enumerate_vertices_up_to;
   otherwise at the extreme vertices of every bus and at the minimiser of its second-order model;
4. Monte Carlo with exact AC samples (containment, 99.73 % range, 3-sigma box).

The interval tested for bus i is  [V* - r_i + lo_i,  V* + r_i + hi_i]  with r_i the affine radius and
[lo_i, hi_i] an enclosure of the quadratic form. It is a second-order bound: the O(delta^3) remainder is
measured, not bounded, and nothing is rounded outward. "greedy_k2" is the simplified binomial pass used
in X1; "raak_*" is the group's operator.

Writes results/x2_summary.csv and results/x2_bus.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, grid, method, pct, seed_for, study
from cdhim import participation_factors, sample_eps, save_table, spectral_zones
from cdhim.knomial import knomial_bound
from cdhim.validation import VIOLATION_ATOL
from cdhim.zones import laplacian_embedding
from x1_aa_second_order import SecondOrder, bound_greedy_binomial, bound_knomial, l1, neg, pos

N_MC = CFG["experiments"]["n_ac"]
BATCH = CFG["powerflow"]["batch"]
ZONES = CFG["experiments"]["zone_counts"]
K_ENUM = CFG["experiments"]["enumerate_vertices_up_to"]
BALANCING = CFG["experiments"]["balancing"]
SWEEP = CFG["experiments"]["x2_delta_sweep"]          # extra half-widths, run with the smallest K only
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]  # cap on the clique search of RAA-K for max_k >= 6
R_FLOOR = 1e-6
BOUNDS = ("natural", "binomial_terms", "monomial", "greedy_k2", "raak_k2", "raak_k3", "raak_k6", "raak_kn",
          "knomial_1", "knomial_8", "spectral")
RAAK = ("raak_k2", "raak_k3", "raak_k6", "raak_kn")
MAIN = ("monomial", "knomial_8", "spectral")          # written per bus (AF2 and the two k-nomial variants)
LEVELS = (1e-4, 1e-3)                                 # p.u.; shortfalls above these are counted separately


# ----------------------------------------------------------------------------- quadratic forms
def quadratic_forms(so):
    """H_i for every PQ bus, shape (n_pq, K, K)."""
    return so.forms()


def all_bounds(so, H, lam, vec):
    """name -> (lo, hi) arrays over the PQ buses."""
    tb = so.term_bounds()
    out = {"natural": tb["natural"], "binomial_terms": tb["binomial"]}
    d = np.einsum("ikk->ik", H)
    off = np.abs(H).sum(axis=(1, 2)) - np.abs(d).sum(axis=1)
    out["monomial"] = (-neg(d).sum(axis=1) - off, pos(d).sum(axis=1) + off)
    s2 = np.abs(vec).sum(axis=1) ** 2                           # |v_k|_1^2 of every eigenvector
    out["spectral"] = (-(neg(lam) * s2).sum(axis=1), (pos(lam) * s2).sum(axis=1))
    k = H.shape[1]
    loops = {"greedy_k2": lambda i: bound_greedy_binomial(H[i]), "knomial_1": lambda i: bound_knomial(H[i], lam[i], vec[i], 1),
             "knomial_8": lambda i: bound_knomial(H[i], lam[i], vec[i], min(8, k)),
             "raak_k2": lambda i: knomial_bound(H[i], 2)[:2], "raak_k3": lambda i: knomial_bound(H[i], 3)[:2],
             "raak_k6": lambda i: knomial_bound(H[i], min(6, k), NODE_BUDGET)[:2],
             "raak_kn": lambda i: knomial_bound(H[i], k, NODE_BUDGET)[:2]}
    for name, fn in loops.items():
        vals = np.array([fn(i) for i in range(len(H))])
        out[name] = (vals[:, 0], vals[:, 1])
    return out


def model_minimiser_vertices(C, H):
    """Vertex where the second-order model c'x + x'Hx of each bus is lowest (1-flip descent from -sign(c))."""
    X = -np.sign(C)
    X[X == 0] = 1.0
    d = np.einsum("ikk->ik", H)
    for _ in range(50 * C.shape[1]):
        HX = np.einsum("ikl,il->ik", H, X)
        gain = -2.0 * X * C - 4.0 * X * HX + 4.0 * d            # change of the model when x_k flips
        k = gain.argmin(axis=1)
        go = gain[np.arange(len(X)), k] < -1e-15
        if not go.any():
            break
        X[go, k[go]] *= -1.0
    return X


def vertex_batches(k):
    """All 2^k vertices of [-1, 1]^k, in batches of columns."""
    for start in range(0, 2 ** k, BATCH):
        idx = np.arange(start, min(start + BATCH, 2 ** k), dtype=np.int64)
        yield 2.0 * ((idx[None, :] >> np.arange(k)[:, None]) & 1) - 1.0


class Scenarios:
    """Streaming statistics of exact AC solutions against the affine box and the second-order intervals."""

    def __init__(self, st, C, Hflat, r, bounds, keep=False):
        self.st, self.C, self.Hflat, self.r = st, C, Hflat, r
        n = len(r)
        self.lo = {"affine": -r, **{b: -r + bounds[b][0] for b in bounds}}
        self.hi = {"affine": r.copy(), **{b: r + bounds[b][1] for b in bounds}}
        self.n = self.n_fail = 0
        self.dv_min, self.dv_max = np.full(n, np.inf), np.full(n, -np.inf)
        self.q_min, self.q_max = np.zeros(n), np.zeros(n)
        self.e2_max, self.e1_max = np.zeros(n), np.zeros(n)
        self.out = {b: np.zeros(n, np.int64) for b in self.lo}
        self.short = {b: np.zeros(n) for b in self.lo}
        self.any = {b: 0 for b in self.lo}
        self.samples = [] if keep else None
        self.seconds = 0.0

    def add(self, eps, delta):
        t0 = perf_counter()
        V, conv, _ = self.st.solver.solve(self.st.injections(eps, delta))
        self.seconds += perf_counter() - t0
        self.n_fail += int((~conv).sum())
        eps = eps[:, conv]
        if not eps.shape[1]:
            return
        self.n += eps.shape[1]
        k = eps.shape[0]
        dv = np.abs(V[self.st.grid.pq][:, conv]) - self.st.vm_star[:, None]
        q = self.Hflat @ (eps[:, None, :] * eps[None, :, :]).reshape(k * k, -1)
        e1 = dv - self.C @ eps
        np.minimum(self.dv_min, dv.min(axis=1), out=self.dv_min)
        np.maximum(self.dv_max, dv.max(axis=1), out=self.dv_max)
        np.minimum(self.q_min, q.min(axis=1), out=self.q_min)
        np.maximum(self.q_max, q.max(axis=1), out=self.q_max)
        np.maximum(self.e1_max, np.abs(e1).max(axis=1), out=self.e1_max)
        np.maximum(self.e2_max, np.abs(e1 - q).max(axis=1), out=self.e2_max)
        active = (self.r > R_FLOOR)[:, None]
        for b in self.lo:
            gap = np.maximum(self.lo[b][:, None] - dv, dv - self.hi[b][:, None])
            bad = (gap > VIOLATION_ATOL) & active
            self.out[b] += bad.sum(axis=1)
            self.any[b] += int(bad.any(axis=0).sum())
            np.maximum(self.short[b], np.where(active, gap, 0.0).max(axis=1), out=self.short[b])
        if self.samples is not None:
            self.samples.append(dv.astype(np.float32))


# ----------------------------------------------------------------------------- one case
def run_case(key, balancing, k, embedding, delta=None):
    sysc = CFG["systems"][key]
    design = delta is None
    delta = sysc["delta"] if design else delta
    g = grid(sysc["case"])
    zones = spectral_zones(g, k, embedding)
    alpha = participation_factors(g) if balancing == "participation" else None
    st = study(sysc["case"], zones=zones, participation=alpha)
    t_first = sum(st.timings.values())
    res = st.run(delta, method())
    C, r = st.coefficients(delta), res.r_affine
    K = st.src.n
    t0 = perf_counter()
    so = SecondOrder(st, delta)
    assert np.allclose(so.c_lin, C, rtol=1e-6, atol=1e-9 * np.abs(C).max()), "rectangular and polar sensitivities differ"
    H = quadratic_forms(so)
    t_second = perf_counter() - t0
    check = sample_eps(np.random.default_rng(seed_for("x2", key, balancing, k, "check")), K, 4)
    Hflat = H.reshape(len(H), K * K)
    q_direct = Hflat @ (check[:, None, :] * check[None, :, :]).reshape(K * K, -1)
    assert np.allclose(q_direct, so.quadratic(check), rtol=1e-7, atol=1e-13), "H_i is inconsistent with the term-by-term form"
    t0 = perf_counter()
    lam, vec = np.linalg.eigh(H)
    bounds = all_bounds(so, H, lam, vec)
    t_bounds = perf_counter() - t0
    nsd = lam.max(axis=1) <= 1e-9 * np.abs(lam).max(axis=1)
    active = r > R_FLOOR
    tag = f"[x2] {key} {balancing} K={K} delta={delta:g}"

    # ---- vertices: all of them, or the adversarial ones
    vx = Scenarios(st, C, Hflat, r, bounds)
    if K <= K_ENUM:
        mode = "all 2^K"
        for i, eps in enumerate(vertex_batches(K)):
            vx.add(eps, delta)
            if i % 16 == 15:
                print(f"{tag}: {vx.n + vx.n_fail}/{2 ** K} vertices", flush=True)
    else:
        mode = "extreme and model-minimiser vertices of every bus"
        lo_v = -np.sign(C)
        lo_v[lo_v == 0] = 1.0
        cand = np.unique(np.vstack([lo_v[active], -lo_v[active], model_minimiser_vertices(C, H)[active]]), axis=0).T
        for s in range(0, cand.shape[1], BATCH):
            vx.add(np.ascontiguousarray(cand[:, s:s + BATCH]), delta)
    print(f"{tag}: {vx.n} vertices ({mode}), {vx.n_fail} without AC solution, {vx.seconds:.0f} s", flush=True)

    # ---- Monte Carlo
    seed = seed_for("x2", key, balancing, k, delta)
    rng = np.random.default_rng(seed)
    mc = Scenarios(st, C, Hflat, r, bounds, keep=True)
    for start in range(0, N_MC, BATCH):
        mc.add(sample_eps(rng, K, min(BATCH, N_MC - start)), delta)
    dv_mc = np.concatenate(mc.samples, axis=1)
    tail = (1.0 - CFG["method"]["confidence"]) / 2.0
    q_lo, q_hi = np.quantile(dv_mc, [tail, 1.0 - tail], axis=1)
    viol_3s = np.abs(dv_mc) > res.r_cltm[:, None] + VIOLATION_ATOL
    print(f"{tag}: {mc.n} AC samples, {mc.n_fail} without solution", flush=True)

    # ---- per-bus table
    rows = []
    for i in np.flatnonzero(active):
        row = {"system": key, "balancing": balancing, "K": K, "delta": delta, "bus": int(g.bus_label[g.pq[i]]), "r_affine": r[i],
               "r_3sigma": res.r_cltm[i], "h_negative_semidefinite": bool(nsd[i]),
               "vertex_dv_min": vx.dv_min[i], "vertex_dv_max": vx.dv_max[i], "vertex_q_min": vx.q_min[i], "vertex_q_max": vx.q_max[i],
               "vertex_lin_err_max": vx.e1_max[i], "vertex_third_order_max": vx.e2_max[i],
               "mc_q00135": q_lo[i], "mc_q99865": q_hi[i], "mc_third_order_max": mc.e2_max[i]}
        for b in BOUNDS:
            row[f"{b}_lo"], row[f"{b}_hi"] = bounds[b][0][i], bounds[b][1][i]
        for b in ("affine",) + MAIN:
            row[f"{b}_vertices_outside"], row[f"{b}_vertex_shortfall"] = int(vx.out[b][i]), vx.short[b][i]
            row[f"{b}_mc_outside"] = int(mc.out[b][i])
        rows.append(row)

    # ---- summary
    a = active
    ref_w = (vx.dv_max - vx.dv_min)[a]
    mc_w = (q_hi - q_lo)[a]
    exc_lo = ((-r - vx.dv_min) / r)[a]
    depth = -vx.q_min                                           # deepest second-order term reached at the vertices
    has = a & (depth > 0)
    exact = a & nsd & (depth > 0) if K <= K_ENUM else np.zeros_like(a)   # exact range of the form: [min over vertices, 0]
    s = {"system": key, "balancing": balancing, "K": K, "delta": delta, "design_delta": design, "n_buses": int(a.sum()),
         "vertex_set": mode, "n_vertices": vx.n, "vertices_without_ac_solution": vx.n_fail, "n_mc": mc.n,
         "mc_without_ac_solution": mc.n_fail, "seed": seed, "first_order_seconds": t_first, "second_order_seconds": t_second,
         "bounds_seconds": t_bounds, "vertex_solve_seconds": vx.seconds,
         "buses_h_negative_semidefinite_pct": pct(nsd[a].mean()),
         "angle_radius_median_rad": float(np.median((l1(so.c_perp) / so.vm)[a])),
         "max_abs_dv_at_vertices_pu": float(np.maximum(-vx.dv_min, vx.dv_max)[a].max()),
         "r_affine_median_pu": float(np.median(r[a])), "r_affine_max_pu": float(r[a].max()),
         "affine_over_3sigma_radius_median": float(np.median((r / res.r_cltm)[a])),
         "affine_width_over_mc9973_median": float(np.median(2.0 * r[a] / mc_w)),
         "sigma3_width_over_mc9973_median": float(np.median(2.0 * res.r_cltm[a] / mc_w)),
         "sigma3_mc_violation_max_pct": pct(viol_3s[a].mean(axis=1).max()),
         "sigma3_mc_system_pct": pct(viol_3s[a].any(axis=0).mean()),
         "sigma3_buses_wider_than_affine": int((res.r_cltm > r)[a].sum()),
         "affine_low_excess_over_r_median_pct": pct(np.median(exc_lo)), "affine_low_excess_over_r_p90_pct": pct(np.quantile(exc_lo, 0.9)),
         "affine_low_excess_over_r_max_pct": pct(exc_lo.max()),
         "second_order_depth_over_r_median_pct": pct(np.median((depth / r)[a])),
         "second_order_depth_over_r_p90_pct": pct(np.quantile((depth / r)[a], 0.9)),
         "second_order_depth_over_r_max_pct": pct((depth / r)[a].max()), "second_order_depth_max_pu": float(depth[a].max()),
         "third_order_over_r_median_pct": pct(np.median((vx.e2_max / r)[a])),
         "third_order_over_r_p90_pct": pct(np.quantile((vx.e2_max / r)[a], 0.9)),
         "third_order_over_r_max_pct": pct((vx.e2_max / r)[a].max()), "third_order_max_pu": float(vx.e2_max[a].max()),
         "third_over_second_order_median_pct": pct(np.median((vx.e2_max / depth)[has])),
         "third_over_second_order_p90_pct": pct(np.quantile((vx.e2_max / depth)[has], 0.9)),
         "n_buses_exact_form_range": int(exact.sum())}
    for b in ("affine",) + BOUNDS:
        s[f"{b}_buses_with_vertex_outside"] = int((vx.out[b][a] > 0).sum())
        s[f"{b}_vertex_shortfall_max_pu"] = float(vx.short[b][a].max())
        s[f"{b}_vertex_shortfall_over_r_max_pct"] = pct((vx.short[b] / r)[a].max())
        for lev in LEVELS:
            s[f"{b}_buses_shortfall_over_{lev:g}pu"] = int((vx.short[b][a] > lev).sum())
        s[f"{b}_mc_buses_outside"] = int((mc.out[b][a] > 0).sum())
        s[f"{b}_mc_samples_outside_pct"] = pct(mc.any[b] / mc.n)
    for b in BOUNDS:
        lo, hi = bounds[b]
        w = (2.0 * r + hi - lo)[a]
        s[f"{b}_width_over_reference_median"] = float(np.median(w / ref_w))
        s[f"{b}_width_over_reference_max"] = float((w / ref_w).max())
        s[f"{b}_width_over_mc9973_median"] = float(np.median(w / mc_w))
        s[f"{b}_remainder_over_r_median_pct"] = pct(np.median(((hi - lo) / r)[a]))
        if exact.any():
            s[f"{b}_form_width_over_exact_median"] = float(np.median((hi - lo)[exact] / depth[exact]))
            s[f"{b}_form_width_over_exact_max"] = float(((hi - lo)[exact] / depth[exact]).max())
            s[f"{b}_positive_side_over_depth_median_pct"] = pct(np.median(hi[exact] / depth[exact]))
    for b in ("greedy_k2",) + RAAK + ("knomial_1", "knomial_8", "spectral"):   # same bus: AF2 width over this bound's width
        s[f"af2_over_{b}_form_width_median"] = float(np.median(((bounds["monomial"][1] - bounds["monomial"][0])
                                                                / (bounds[b][1] - bounds[b][0]))[a]))
    s["reference_width_over_affine_width_median"] = float(np.median(ref_w / (2.0 * r[a])))
    for b in ("spectral", "monomial", "raak_kn"):                         # the gate criteria of the checklist
        s[f"gate1_contains_all_{b}"] = bool(s[f"{b}_buses_with_vertex_outside"] == 0 and s[f"{b}_mc_buses_outside"] == 0)
        s[f"gate2_width_ok_{b}"] = bool(s[f"{b}_width_over_reference_median"] <= 1.5)
    return rows, s


def main():
    bus, summ = [], []
    for key, counts in ZONES.items():
        emb = laplacian_embedding(grid(CFG["systems"][key]["case"]))
        for balancing in BALANCING[key]:
            cases = [(k, None) for k in counts] + [(min(counts), d) for d in SWEEP[key]]
            for k, delta in cases:
                rows, s = run_case(key, balancing, k, emb, delta)
                bus += rows
                summ.append(s)
                print("  " + "  ".join(f"{n}={s[n]}" for n in ("gate1_contains_all_monomial", "gate1_contains_all_spectral",
                                                                "spectral_vertex_shortfall_max_pu", "spectral_width_over_reference_median",
                                                                "af2_over_spectral_form_width_median")), flush=True)
    meta = {"n_mc": N_MC, "enumerate_vertices_up_to": K_ENUM,
            "seeds": {f"{s['system']}/{s['balancing']}/K={s['K']}/delta={s['delta']:g}": s["seed"] for s in summ},
            "raak_node_budget": NODE_BUDGET,
            "note": "second-order bound; third order measured, not bounded; no outward rounding; raak_* = group's RAA-K (cdhim.knomial)"}
    print(save_table("x2_bus", bus, meta))
    print(save_table("x2_summary", summ, meta))


if __name__ == "__main__":
    main()
