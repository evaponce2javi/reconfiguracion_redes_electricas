"""X1 (exploratory, not requested by the reviewers) -- what affine arithmetic can add to the method.

In rectangular coordinates x = (e, f) the power-flow equations are exactly quadratic,
    g(x* + D) = g(x*) + J D + q(D),     q_m(D) = dV_m conj((Y dV)_m)   (|dV_m|^2 at PV buses),
so the deviation of the state solves D = C eps - J^-1 q(D) with no truncated higher-order terms and
    d|V_i| = c_i' eps + eps' H_i eps + O(delta^3).
The first term is the affine form the method already computes. The second is a quadratic form in the
noise symbols: the object that AF2 / RAA bilinear operators, diagonal budgeting and the k-nomial
theorem bound. The script measures

1. how much of the AC linearisation error the quadratic form explains (all PQ buses);
2. enclosures of its range that need no H_i (all buses): term-by-term products of affine forms
   ("natural" affine arithmetic) and the same terms as scaled binomial squares (k-nomial blocks);
3. enclosures from the assembled H_i (a subset of buses): monomial-wise (k = 1), greedy binomial
   blocks on the diagonal budget (k = 2), spectral blocks (k = n), against an inner estimate;
4. the exact AC solution at the two vertices of the noise box where the affine form of each bus is
   extreme: is the affine box an enclosure of the AC solution set?
5. probabilistic intervals with and without the second-order centre shift tr(H_i)/3, and a
   distribution-free (Hoeffding) radius;
6. the dependency that the shared noise symbols carry on each branch: zonotope of (V_i, V_j) against
   its box, and the radius of the voltage difference V_i - V_j.

Only the leading (second-order) part of the remainder is bounded: this is a feasibility study, not a
verified enclosure. A verified one needs the fixed point  r = -J^-1 q(C eps + r)  and outward rounding.

Writes results/x1_bus_remainder.csv, x1_quadratic_bounds.csv, x1_branch_dependency.csv, x1_summary.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, method, pct, seed_for, study
from cdhim import geometry as geo
from cdhim import sample_eps, save_table
from cdhim.validation import VIOLATION_ATOL

N_MC = CFG["experiments"]["n_ac"]
BATCH = CFG["powerflow"]["batch"]
ALPHA = 1.0 - CFG["method"]["confidence"]
N_FORMS = {"ieee33": None, "pegase2869": (8, 16)}   # (largest error, random)       # buses whose H_i is assembled (None = all)
R_FLOOR = 1e-6                                     # p.u.; buses below it are left out of the statistics
CLTM_TERM, CLTM_M = 1e-6, 6                        # Ref. [10]: a term is effective above 1e-6; CLT needs M = 6 of them


from cdhim.second_order import (SecondOrder, bound_greedy_binomial, bound_knomial, bound_monomial,  # noqa: F401
                                bound_spectral, l1, neg, pos)


def _ascend(H, x, max_steps=200000):
    """Greedy coordinate ascent of x'Hx on [-1, 1]^n (exact 1-D maximiser, interior points allowed)."""
    d = np.diag(H).copy()
    u = H @ x
    concave = d < 0
    dc = np.where(concave, d, -1.0)
    for _ in range(max_steps):
        up, dn = 1.0 - x, -1.0 - x
        g_up, g_dn = (2.0 * u + d * up) * up, (2.0 * u + d * dn) * dn
        step = np.where(g_up >= g_dn, up, dn)
        gain = np.maximum(g_up, g_dn)
        t = np.clip(x - u / dc, -1.0, 1.0) - x
        g_in = (2.0 * u + d * t) * t
        use = concave & (g_in > gain)
        step, gain = np.where(use, t, step), np.where(use, g_in, gain)
        k = int(np.argmax(gain))
        if gain[k] <= 1e-14 * max(1.0, abs(float(x @ u))):
            break
        x[k] += step[k]
        u += step[k] * H[:, k]
    return float(x @ H @ x)


def inner_range(H, lam, vec, rng, n_start=2):
    """Attained values of eps' H eps: a lower estimate of the width of the exact range."""
    n = len(H)
    lo = hi = 0.0                                            # eps = 0
    starts_hi = [np.sign(vec[:, -1]) + (vec[:, -1] == 0)] + [np.sign(rng.standard_normal(n)) for _ in range(n_start)]
    starts_lo = [np.sign(vec[:, 0]) + (vec[:, 0] == 0)] + [np.sign(rng.standard_normal(n)) for _ in range(n_start)]
    for x in starts_hi:
        hi = max(hi, _ascend(H, x.astype(float)))
    for x in starts_lo:
        lo = min(lo, -_ascend(-H, x.astype(float)))
    return lo, hi


# ----------------------------------------------------------------------------- branch dependency
def zonotope_area(a, b):
    """Area of {(a'eps, b'eps): eps in [-1, 1]^n} for every row of (a, b), by sorting the generators."""
    flip = (b < 0) | ((b == 0) & (a < 0))
    a, b = np.where(flip, -a, a), np.where(flip, -b, b)
    order = np.argsort(np.arctan2(b, a), axis=1, kind="stable")
    a, b = np.take_along_axis(a, order, 1), np.take_along_axis(b, order, 1)
    px = np.concatenate([np.zeros((len(a), 1)), 2.0 * np.cumsum(a, axis=1)], axis=1) - a.sum(axis=1, keepdims=True)
    py = np.concatenate([np.zeros((len(b), 1)), 2.0 * np.cumsum(b, axis=1)], axis=1) - b.sum(axis=1, keepdims=True)
    return np.abs((px[:, :-1] * py[:, 1:] - px[:, 1:] * py[:, :-1]).sum(axis=1))


def branch_dependency(key, C, res):
    rows = []
    ok = np.flatnonzero(res.valid)
    r = res.r_affine
    for s in range(0, len(ok), 256):
        b = ok[s:s + 256]
        i, j = res.fi[b], res.ti[b]
        a_zon = zonotope_area(C[i], C[j])
        a_box = geo.box_area(r[i], r[j])
        a_ell = geo.ellipse_area(res.lam_min[b], res.lam_max[b], res.chi2)
        diff = C[i] - C[j]
        r_diff, s_diff = l1(diff), np.sqrt((diff ** 2).sum(axis=1) / 3.0)
        for t in range(len(b)):
            rows.append({"system": key, "branch": int(b[t]), "rho": res.rho[b[t]], "cdhim_gain_pct": pct(res.gain[b[t]]),
                         "cdhim_mode": "ellipse" if res.ellipse[b[t]] else "box",
                         "area_affine_box": a_box[t], "area_zonotope": a_zon[t], "area_ellipse_9973": a_ell[t],
                         "zonotope_gain_pct": pct(1.0 - a_zon[t] / a_box[t]),
                         "drop_radius_interval": r[i[t]] + r[j[t]], "drop_radius_affine": r_diff[t],
                         "drop_radius_3sigma": 3.0 * s_diff[t],
                         "drop_reduction_pct": pct(1.0 - r_diff[t] / (r[i[t]] + r[j[t]]))})
    return rows


# ----------------------------------------------------------------------------- extreme vertices
def vertex_check(st, so, C, delta):
    """Exact AC |V_i| at eps = -sign(c_i) and +sign(c_i), where the affine form of bus i is extreme."""
    n_pq = C.shape[0]
    out = {}
    for side, sgn in (("lo", -1.0), ("hi", 1.0)):
        dv, q2, ok = np.full(n_pq, np.nan), np.empty(n_pq), np.zeros(n_pq, bool)
        for s in range(0, n_pq, BATCH):
            idx = np.arange(s, min(s + BATCH, n_pq))
            cols = np.arange(len(idx))
            eps = sgn * np.sign(C[idx]).T
            V, conv, _ = st.solver.solve(st.injections(eps, delta))
            dv[idx] = np.abs(V[st.grid.pq[idx], cols]) - st.vm_star[idx]
            q2[idx] = so.quadratic(eps)[idx, cols]
            ok[idx] = conv
        out[side] = (dv, q2, ok)
    return out


# ----------------------------------------------------------------------------- driver
def run_system(key):
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    st = study(sysc["case"])
    g = st.grid
    res = st.run(delta, method())
    t0 = perf_counter()
    so = SecondOrder(st, delta)
    C = st.coefficients(delta)
    scale = np.abs(C).max()
    assert np.allclose(so.c_lin, C, rtol=1e-6, atol=1e-9 * scale), "rectangular and polar sensitivities differ"
    tb = so.term_bounds()
    shift = so.trace() / 3.0                                     # E[eps' H eps] for independent U[-1, 1]
    t_model = perf_counter() - t0
    n_pq = len(g.pq)
    radii = {"affine": res.r_affine, "sigma3": res.r_cltm, "coord": res.r_coord,
             "hoeffding": np.sqrt(2.0 * np.log(2.0 / ALPHA) * (C ** 2).sum(axis=1))}
    lo_b, hi_b = tb["binomial"]
    vx = vertex_check(st, so, C, delta)
    n_terms = (np.abs(C) > CLTM_TERM).sum(axis=1)

    seed = seed_for("x1", key)
    rng = np.random.default_rng(seed)
    e1_max, e2_max = np.zeros(n_pq), np.zeros(n_pq)
    e1_min_s, e1_max_s, e1_sum = np.full(n_pq, np.inf), np.full(n_pq, -np.inf), np.zeros(n_pq)
    viol = {f"{name}_{c}": np.zeros(n_pq, np.int64) for name in radii for c in ("plain", "shift")}
    viol["affine_plus_binomial"] = np.zeros(n_pq, np.int64)
    any_viol = {k: 0 for k in viol}
    n_ac = 0
    for start in range(0, N_MC, BATCH):
        eps = sample_eps(rng, st.src.n, min(BATCH, N_MC - start))
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        eps = eps[:, conv]
        dv = np.abs(V[g.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        e1 = dv - C @ eps
        e2 = e1 - so.quadratic(eps)
        np.maximum(e1_max, np.abs(e1).max(axis=1), out=e1_max)
        np.maximum(e2_max, np.abs(e2).max(axis=1), out=e2_max)
        np.minimum(e1_min_s, e1.min(axis=1), out=e1_min_s)
        np.maximum(e1_max_s, e1.max(axis=1), out=e1_max_s)
        e1_sum += e1.sum(axis=1)
        for name, r in radii.items():
            for c, centre in (("plain", 0.0), ("shift", shift[:, None])):
                out = np.abs(dv - centre) > r[:, None] + VIOLATION_ATOL
                viol[f"{name}_{c}"] += out.sum(axis=1)
                any_viol[f"{name}_{c}"] += int(out.any(axis=0).sum())
        out = (dv > (res.r_affine + hi_b)[:, None] + VIOLATION_ATOL) | (dv < (-res.r_affine + lo_b)[:, None] - VIOLATION_ATOL)
        viol["affine_plus_binomial"] += out.sum(axis=1)
        any_viol["affine_plus_binomial"] += int(out.any(axis=0).sum())
        print(f"[x1] {key}: {start + eps.shape[1]}/{N_MC} AC samples", flush=True)

    bus_rows = []
    for k in range(n_pq):
        row = {"system": key, "bus": int(g.bus_label[g.pq[k]]), "r_affine": res.r_affine[k], "r_3sigma": res.r_cltm[k],
               "r_coord": res.r_coord[k], "r_hoeffding": radii["hoeffding"][k],
               "lin_err_max": e1_max[k], "second_order_err_max": e2_max[k],
               "remainder_mc_min": e1_min_s[k], "remainder_mc_max": e1_max_s[k], "remainder_mc_mean": e1_sum[k] / n_ac,
               "centre_shift_predicted": shift[k], "angle_radius_rad": float(np.abs(so.c_perp[k]).sum() / so.vm[k]),
               "n_terms_above_1e-6": int(n_terms[k])}
        for side in ("lo", "hi"):
            dv_v, q_v, ok_v = vx[side]
            sgn = -1.0 if side == "lo" else 1.0
            row[f"vertex_{side}_converged"] = bool(ok_v[k])
            row[f"vertex_{side}_dv_ac"] = dv_v[k]
            row[f"vertex_{side}_dv_affine"] = sgn * res.r_affine[k]
            row[f"vertex_{side}_dv_second_order"] = sgn * res.r_affine[k] + q_v[k]
        for name, part in so.natural_parts.items():
            row[f"natural_width_{name}"] = part[k]
        for name, (lo, hi) in tb.items():
            row[f"quad_{name}_lo"], row[f"quad_{name}_hi"] = lo[k], hi[k]
        row.update({f"viol_{name}_pct": pct(v[k] / n_ac) for name, v in viol.items()})
        bus_rows.append(row)

    # ---- assembled quadratic forms on a subset of buses
    active = np.flatnonzero(res.r_affine > R_FLOOR)
    if N_FORMS[key] is None:
        subset, why = active, ["all"] * len(active)
    else:
        n_worst, n_rnd = N_FORMS[key]
        worst = active[np.argsort(-(e1_max[active] / res.r_affine[active]))[:n_worst]]
        rest = np.setdiff1d(active, worst)
        rnd = np.random.default_rng(seed_for("x1", key, "subset")).choice(rest, n_rnd, replace=False)
        subset, why = np.concatenate([worst, rnd]), ["largest relative linearisation error"] * n_worst + ["random"] * n_rnd
    form_rows = []
    rng_in = np.random.default_rng(seed_for("x1", key, "inner"))
    check = sample_eps(np.random.default_rng(seed_for("x1", key, "check")), st.src.n, 3)
    q_check = so.quadratic(check)
    for k, reason in zip(subset, why):
        H = so.assemble(k)
        assert np.allclose(np.einsum("is,is->s", check, H @ check), q_check[k], rtol=1e-8, atol=1e-14), "H_i is inconsistent"
        lam, vec = np.linalg.eigh(H)
        in_lo, in_hi = inner_range(H, lam, vec, rng_in)
        row = {"system": key, "bus": int(g.bus_label[g.pq[k]]), "selected_as": reason, "n_symbols": len(H),
               "r_affine": res.r_affine[k], "lin_err_max": e1_max[k], "second_order_err_max": e2_max[k],
               "inner_lo": in_lo, "inner_hi": in_hi, "trace": float(np.trace(H)),
               "eig_for_90pct_of_abs_spectrum": int(np.searchsorted(np.cumsum(np.sort(np.abs(lam))[::-1]) / np.abs(lam).sum(), 0.9) + 1)}
        bounds = {"natural": (tb["natural"][0][k], tb["natural"][1][k]), "binomial_terms": (lo_b[k], hi_b[k]),
                  "monomial": bound_monomial(H), "greedy_k2": bound_greedy_binomial(H),
                  "knomial_1": bound_knomial(H, lam, vec, 1), "knomial_8": bound_knomial(H, lam, vec, 8),
                  "spectral": bound_spectral(lam, vec)}
        for name, (lo, hi) in bounds.items():
            tol = 1e-9 * (in_hi - in_lo)
            assert lo <= in_lo + tol and hi >= in_hi - tol, f"{name} does not enclose an attained value"
            row[f"{name}_lo"], row[f"{name}_hi"] = lo, hi
            row[f"{name}_over_inner"] = (hi - lo) / (in_hi - in_lo)
        row["depth_over_r_affine"] = -in_lo / res.r_affine[k]            # attained second-order sag / affine radius
        row["negative_share_of_abs_spectrum"] = float(neg(lam).sum() / np.abs(lam).sum())
        form_rows.append(row)
        print(f"[x1] {key}: form of bus {row['bus']} done", flush=True)

    dep_rows = branch_dependency(key, C, res)

    # ---- summary
    import pandas as pd
    b, f, d = pd.DataFrame(bus_rows), pd.DataFrame(form_rows), pd.DataFrame(dep_rows)
    a = b[b.r_affine > R_FLOOR]
    w_nat, w_bin = a.quad_natural_hi - a.quad_natural_lo, a.quad_binomial_hi - a.quad_binomial_lo
    summ = {"system": key, "delta": delta, "n_pq": n_pq, "n_buses_above_floor": len(a), "n_symbols": st.src.n,
            "n_ac": n_ac, "n_not_converged": N_MC - n_ac, "seed": seed, "model_seconds": t_model,
            "lin_err_over_r_affine_median_pct": pct((a.lin_err_max / a.r_affine).median()),
            "lin_err_over_r_affine_max_pct": pct((a.lin_err_max / a.r_affine).max()),
            "lin_err_over_r_3sigma_median_pct": pct((a.lin_err_max / a.r_3sigma).median()),
            "lin_err_over_r_3sigma_max_pct": pct((a.lin_err_max / a.r_3sigma).max()),
            "err_reduction_by_quadratic_median": float((a.lin_err_max / a.second_order_err_max).median()),
            "err_reduction_by_quadratic_min": float((a.lin_err_max / a.second_order_err_max).min()),
            "natural_width_over_r_affine_median_pct": pct((w_nat / a.r_affine).median()),
            "natural_width_over_r_affine_max_pct": pct((w_nat / a.r_affine).max()),
            "binomial_width_over_r_affine_median_pct": pct((w_bin / a.r_affine).median()),
            "binomial_width_over_r_affine_max_pct": pct((w_bin / a.r_affine).max()),
            "binomial_over_natural_width_median": float((w_bin / w_nat).median()),
            "buses_remainder_mc_outside_binomial": int(((a.remainder_mc_min < a.quad_binomial_lo - VIOLATION_ATOL)
                                                         | (a.remainder_mc_max > a.quad_binomial_hi + VIOLATION_ATOL)).sum()),
            "shift_vs_mc_mean_max_abs_diff_over_r_3sigma_pct": pct(((a.centre_shift_predicted - a.remainder_mc_mean).abs() / a.r_3sigma).max()),
            "shift_over_r_3sigma_median_pct": pct((a.centre_shift_predicted.abs() / a.r_3sigma).median()),
            "shift_over_r_3sigma_max_pct": pct((a.centre_shift_predicted.abs() / a.r_3sigma).max()),
            "hoeffding_over_affine_median": float((a.r_hoeffding / a.r_affine).median()),
            "hoeffding_over_3sigma": float(np.sqrt(2.0 * np.log(2.0 / ALPHA) * 3.0) / 3.0)}
    # exact AC at the extreme vertex of each bus: positive excess = AC solution outside the affine box
    vl, vh = a[a.vertex_lo_converged], a[a.vertex_hi_converged]
    exc_lo = (vl.vertex_lo_dv_affine - vl.vertex_lo_dv_ac) / vl.r_affine
    exc_hi = (vh.vertex_hi_dv_ac - vh.vertex_hi_dv_affine) / vh.r_affine
    pred_lo = (vl.vertex_lo_dv_ac - vl.vertex_lo_dv_second_order).abs() / (vl.vertex_lo_dv_ac - vl.vertex_lo_dv_affine).abs()
    parts = a[["natural_width_products", "natural_width_pv_squares", "natural_width_curvature"]]
    share = parts.div(parts.sum(axis=1), axis=0)
    summ.update({
        "vertex_lo_not_converged": int((~a.vertex_lo_converged).sum()), "vertex_hi_not_converged": int((~a.vertex_hi_converged).sum()),
        "vertex_lo_buses_outside_affine_box": int((vl.vertex_lo_dv_ac < vl.vertex_lo_dv_affine - VIOLATION_ATOL).sum()),
        "vertex_lo_excess_over_r_median_pct": pct(exc_lo.median()), "vertex_lo_excess_over_r_p90_pct": pct(exc_lo.quantile(0.9)),
        "vertex_lo_excess_over_r_max_pct": pct(exc_lo.max()),
        "vertex_lo_buses_excess_over_10pct": int((exc_lo > 0.10).sum()), "vertex_lo_buses_excess_over_100pct": int((exc_lo > 1.0).sum()),
        "vertex_lo_second_order_misses_median_pct_of_excess": pct(pred_lo.median()),
        "vertex_lo_second_order_misses_p90_pct_of_excess": pct(pred_lo.quantile(0.9)),
        "vertex_hi_buses_outside_affine_box": int((vh.vertex_hi_dv_ac > vh.vertex_hi_dv_affine + VIOLATION_ATOL).sum()),
        "vertex_hi_excess_over_r_max_pct": pct(exc_hi.max()), "vertex_hi_slack_over_r_median_pct": pct((-exc_hi).median()),
        "natural_share_products_median_pct": pct(share.natural_width_products.median()),
        "natural_share_pv_squares_median_pct": pct(share.natural_width_pv_squares.median()),
        "natural_share_curvature_median_pct": pct(share.natural_width_curvature.median()),
        "angle_radius_median_rad": float(a.angle_radius_rad.median()), "angle_radius_max_rad": float(a.angle_radius_rad.max()),
        "magnitude_radius_median_pu": float(a.r_affine.median()),
        "cltm_buses_below_6_effective_terms": int((b["n_terms_above_1e-6"] < CLTM_M).sum()),
        "cltm_buses_below_6_effective_terms_above_floor": int((a["n_terms_above_1e-6"] < CLTM_M).sum())})
    for name in viol:
        col = a[f"viol_{name}_pct"]
        summ[f"viol_{name}_mean_pct"], summ[f"viol_{name}_max_pct"] = col.mean(), col.max()
        summ[f"viol_{name}_buses_over_nominal"] = int((col > pct(ALPHA)).sum())
        summ[f"viol_{name}_system_pct"] = pct(any_viol[name] / n_ac)
    for name in ("natural", "binomial_terms", "monomial", "greedy_k2", "knomial_1", "knomial_8", "spectral"):
        summ[f"forms_{name}_over_inner_median"] = float(f[f"{name}_over_inner"].median())
        summ[f"forms_{name}_over_inner_max"] = float(f[f"{name}_over_inner"].max())
        summ[f"forms_{name}_spurious_rise_over_depth_median_pct"] = pct((f[f"{name}_hi"] / -f.inner_lo).median())
    for group, sub in f.groupby("selected_as"):
        tag = {"all": "all", "random": "random", "largest relative linearisation error": "worst"}[group]
        summ[f"forms_{tag}_n"] = len(sub)
        summ[f"forms_{tag}_depth_over_r_affine_median_pct"] = pct(sub.depth_over_r_affine.median())
        summ[f"forms_{tag}_depth_over_r_affine_max_pct"] = pct(sub.depth_over_r_affine.max())
    summ["forms_true_rise_over_depth_max_pct"] = pct((f.inner_hi / -f.inner_lo).max())
    summ["forms_negative_share_of_spectrum_median_pct"] = pct(f.negative_share_of_abs_spectrum.median())
    summ.update({"forms_n": len(f), "forms_eig90_median": float(f.eig_for_90pct_of_abs_spectrum.median()),
                 "zonotope_gain_median_pct": float(d.zonotope_gain_pct.median()),
                 "zonotope_gain_min_pct": float(d.zonotope_gain_pct.min()),
                 "zonotope_over_ellipse_area_median": float((d.area_zonotope / d.area_ellipse_9973).median()),
                 "branches_zonotope_smaller_than_ellipse": int((d.area_zonotope < d.area_ellipse_9973).sum()),
                 "n_branches": len(d),
                 "drop_reduction_median_pct": float(d.drop_reduction_pct.median()),
                 "drop_reduction_min_pct": float(d.drop_reduction_pct.min()),
                 "drop_reduction_max_pct": float(d.drop_reduction_pct.max())})
    return bus_rows, form_rows, dep_rows, summ


def main():
    bus, forms, dep, summ = [], [], [], []
    for key in CFG["systems"]:
        b, f, d, s = run_system(key)
        bus += b
        forms += f
        dep += d
        summ.append(s)
    meta = {"n_mc": N_MC, "seeds": {s["system"]: s["seed"] for s in summ}, "confidence": CFG["method"]["confidence"],
            "note": "exploratory; second-order (leading) remainder only, floating point without outward rounding"}
    for name, rows in (("x1_bus_remainder", bus), ("x1_quadratic_bounds", forms), ("x1_branch_dependency", dep), ("x1_summary", summ)):
        print(save_table(name, rows, meta))
    for s in summ:
        print("\n".join(f"  {k:55s} {v}" for k, v in s.items()))


if __name__ == "__main__":
    main()
