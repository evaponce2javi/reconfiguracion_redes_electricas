"""X6 -- Contribution 1: can a fixed-point argument turn the second-order bound into an enclosure?

See cdhim/enclosure.py for the mathematics. For each case and half-width delta:

1. spectral radius of G = sum_k |J^-1 dJ(C_k)|: the linear part of the remainder map contracts only below 1;
2. the smallest rho with rho >= b2 + G rho + N(rho), if the iteration converges (then, by Brouwer, a
   solution with |D - C eps| <= rho exists for every eps of the box);
3. the enclosure of |V_i| that follows: second-order form + allowance for the higher orders;
4. a check against the exact AC solution at all 2^K vertices: the state remainder must be inside rho, the
   part above second order inside sigma, and |V_i| inside the enclosure. A violation would mean an error in
   the derivation or in the code, because the test is a theorem in exact arithmetic.

The enclosure is not verified: nothing is rounded outward.

Writes results/x6_fixed_point.csv
"""
import numpy as np

from common import CFG, grid, method, pct, study
from cdhim import SecondOrder, participation_factors, save_table, spectral_zones
from cdhim.enclosure import Enclosure, monomial_bounds, polar_contraction_radius, polar_true_radius
from cdhim.knomial import knomial_bound
from cdhim.validation import VIOLATION_ATOL
from cdhim.zones import laplacian_embedding
from x2_few_sources import R_FLOOR, vertex_batches

K = 8
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]
# (case, balancing models, half-widths)
CASES = [("case33bw", ("slack",), (0.05, 0.15, 0.30)),
         ("case118", ("slack", "participation"), (0.02, 0.05, 0.10)),
         ("case300", ("slack", "participation"), (0.02, 0.05, 0.10)),
         ("pglib_opf_case2869_pegase", ("slack", "participation"), (0.01, 0.02, 0.05))]


def run_case(case, balancing, deltas):
    g = grid(case)
    alpha = participation_factors(g) if balancing == "participation" else None
    st = study(case, zones=spectral_zones(g, K, laplacian_embedding(g)), participation=alpha)
    d0 = deltas[-1]
    so = SecondOrder(st, d0)
    enc = Enclosure(so)
    true0 = enc.true_radius()                                    # at d0; linear in delta, like G
    polar0 = polar_contraction_radius(st, d0)
    polar_true = {d: polar_true_radius(st, d) for d in deltas}   # not linear in delta: computed at each half-width
    H0 = so.forms()
    C0 = st.coefficients(d0)
    m = so.m
    x_star = np.concatenate([st.V_star[g.pvpq].real, st.V_star[g.pvpq].imag])
    raak0 = np.array([knomial_bound(h, H0.shape[1], NODE_BUDGET)[:2] for h in H0])
    rows = []
    for delta in deltas:
        sc = delta / d0
        fp = enc.solve(sc)
        C, H = C0 * sc, H0 * sc * sc
        r = np.abs(C).sum(axis=1)
        act = r > R_FLOOR
        allow, e3, s_mag = enc.modulus_allowance(fp, sc)
        lo2 = {"af2": monomial_bounds(H), "raak": (raak0[:, 0] * sc * sc, raak0[:, 1] * sc * sc)}
        # exact AC at every vertex
        dv_min, dv_max = np.full(len(r), np.inf), np.full(len(r), -np.inf)
        out_state = out_sigma = 0
        worst_state = worst_sigma = 0.0
        r_true = np.zeros(2 * m)
        s_true = np.zeros(2 * m)
        n_fail = 0
        Hs = enc.Hs.reshape(2 * m, -1) * sc * sc
        for eps in vertex_batches(st.src.n):
            V, conv, _ = st.solver.solve(st.injections(eps, delta))
            n_fail += int((~conv).sum())
            eps = eps[:, conv]
            dv = np.abs(V[g.pq][:, conv]) - st.vm_star[:, None]
            dv_min, dv_max = np.minimum(dv_min, dv.min(axis=1)), np.maximum(dv_max, dv.max(axis=1))
            D = np.vstack([V[g.pvpq][:, conv].real, V[g.pvpq][:, conv].imag]) - x_star[:, None]
            rem = D - (so.C * sc) @ eps
            third = rem - Hs @ (eps[:, None, :] * eps[None, :, :]).reshape(st.src.n ** 2, -1)
            r_true, s_true = np.maximum(r_true, np.abs(rem).max(axis=1)), np.maximum(s_true, np.abs(third).max(axis=1))
        row = {"case": case, "balancing": balancing, "K": st.src.n, "delta": delta, "n_states": 2 * m, "n_buses": int(act.sum()),
               "vertices_without_ac_solution": n_fail,
               "spectral_radius_G": fp.spectral_radius, "spectral_radius_G_polar_first_order": polar0 * sc,
               "spectral_radius_true_max_over_64_vertices": true0 * sc,
               "spectral_radius_true_polar_max_over_64_vertices": polar_true[delta],
               "quadratic_part": "assembled" if enc.quad is not None else "interval",
               "fixed_point_found": bool(fp.converged), "iterations": fp.iterations,
               "angle_radius_median_rad": float(np.median((np.abs(so.c_perp).sum(axis=1) * sc / so.vm)[act])),
               "state_remainder_true_max_pu": float(r_true.max()), "state_third_order_true_max_pu": float(s_true.max()),
               "b2_max_pu": float(fp.b2.max())}
        if fp.converged:
            row.update({"rho_max_pu": float(fp.rho.max()), "sigma_max_pu": float(fp.sigma.max()),
                        "rho_over_true_remainder_median": float(np.median((fp.rho / np.maximum(r_true, 1e-300))[r_true > 1e-12])),
                        "sigma_over_true_third_order_median": float(np.median((fp.sigma / np.maximum(s_true, 1e-300))[s_true > 1e-12])),
                        "states_outside_rho": int((r_true > fp.rho + VIOLATION_ATOL).sum()),
                        "states_outside_sigma": int((s_true > fp.sigma + VIOLATION_ATOL).sum()),
                        "buses_without_allowance": int((~np.isfinite(allow[act])).sum()),
                        "allowance_over_r_median_pct": pct(np.median((allow / r)[act])), "allowance_over_r_max_pct": pct((allow / r)[act].max()),
                        "allowance_modulus_share_median_pct": pct(np.median((e3 / allow)[act])),
                        "third_order_true_over_r_max_pct": pct(((np.maximum(-(dv_min + r - lo2["af2"][0]), 0.0)) / r)[act].max())})
            ref_w = (dv_max - dv_min)[act]
            for name, (lo, hi) in lo2.items():
                low, up = -r + lo - allow, r + hi + allow
                short = np.maximum(low - dv_min, dv_max - up)[act]
                row[f"{name}_enclosure_buses_violated"] = int((short > VIOLATION_ATOL).sum())
                row[f"{name}_enclosure_width_over_reference_median"] = float(np.median(((up - low)[act]) / ref_w))
                row[f"{name}_enclosure_width_over_reference_max"] = float((((up - low)[act]) / ref_w).max())
                low2, up2 = -r + lo, r + hi
                row[f"{name}_second_order_only_buses_violated"] = int((np.maximum(low2 - dv_min, dv_max - up2)[act] > VIOLATION_ATOL).sum())
                row[f"{name}_second_order_only_width_over_reference_median"] = float(np.median(((up2 - low2)[act]) / ref_w))
        rows.append(row)
        print(f"[x6] {case} {balancing} delta={delta:g}: rho(G)={fp.spectral_radius:.3f} polar={polar0 * sc:.3f} true={true0 * sc:.3f} true polar={polar_true[delta]:.3f} "
              f"fixed point={fp.converged}"
              + (f" allowance/r median {row['allowance_over_r_median_pct']:.3g} % width/ref {row['raak_enclosure_width_over_reference_median']:.3f}"
                 f" violated {row['raak_enclosure_buses_violated']}" if fp.converged else ""), flush=True)
    return rows


def main():
    rows = []
    for case, models, deltas in CASES:
        for balancing in models:
            rows += run_case(case, balancing, deltas)
    print(save_table("x6_fixed_point", rows, {"K": K, "raak_node_budget": NODE_BUDGET,
                                              "note": "deterministic; valid in exact arithmetic, no outward rounding; reference = exact AC at all 2^K vertices"}))


if __name__ == "__main__":
    main()
