"""X14 -- Voltage-controlled (PV) buses: the verified enclosure in voltage-current form (cdhim.verified_vi).

    distribution   the 33-bus networks of X8 (radial, with ties, with shunts) with the sources of the X8 rule turned
                   into voltage-controlled generators. Set-up: every source regulates the voltage magnitude of its
                   bus to the value it has when all sources inject the centre of their X8 range with no reactive
                   power (s0 + C/2); its active power is uncertain over the same interval as in X8,
                   [certain part, certain part + C], one symbol per source; its reactive power is whatever the
                   voltage control needs (no limits). Scenarios F and R and the shares rho_U of X8.
    transmission   case30, case118 and case300 of pandapower with uncertain demand in 8 spectral zones (as in X2 to
                   X7) and half-widths of 2, 5, 10 and 20 %.
    forms          the 80 half-disk cases of X8 (PQ networks) with the current form (cdhim.verified) and with the
                   voltage-current form: where each one closes and the ratio of their radii.

Per case: verified / unique / reason; width over the exact range at the vertices of the box when K <= 10 (the
extremes need not be at vertices: it is a range attained, not the exact one); violations against those vertices;
and two diagnostics that are NOT rigorous (plain floating point):
  linear_radius_G   spectral radius of G = sum_k |J^-1 dJ(E_k)| of cdhim.enclosure.Enclosure: below 1 the linear
                    part of the remainder contracts in the entrywise majorant;
  majorant radii    spectral radius of the linear part of the map of the radii (z, zeta) that verify_tube_vi
                    iterates, [[R^V + T^V D1, T^V D2], [R^I + T^I D1, T^I D2]] with D1 = diag(Ibar (+ 2 dl at PV
                    buses)) and D2 = diag(dl). At or above 1 that iteration has no fixed point however small the
                    defect is. Reported with and without the rounding terms R^V, R^I, and for its two diagonal
                    blocks.

Writes results/x14_pv_distribution.csv, results/x14_pv_transmission.csv and results/x14_forms.csv.

    python experiments/x14_pv.py                  # the three tables
    python experiments/x14_pv.py transmission     # only the named ones
"""
import sys
from dataclasses import replace
from time import perf_counter

import numpy as np
import scipy.sparse as sp

from common import CFG, grid, pct, study
from cdhim import SecondOrder, load_grid, save_table, spectral_zones
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import injection_changes, rect_jacobian, state_forms
from cdhim.halfdisk import half_disk_study
from cdhim.powerflow import newton_pf, solve_blocked
from cdhim.study import custom_study
from cdhim.validation import VIOLATION_ATOL
from cdhim.verified import _cform_sup, _sparse_abs, abs_up, add_up, err_const, matmul_up, mul_up, sum_up, up, verified_enclosure
from cdhim.verified_vi import block_frobenius_up, block_spectral_up, verified_enclosure_vi
from cdhim.zones import laplacian_embedding
from x8_feeders import NETWORKS, NODE_BUDGET, RHO_U, TOL

TRANSMISSION = ("case30", "case118", "case300")
ZONES = 8
DELTAS = (0.02, 0.05, 0.10, 0.20)
MAX_VERTEX_SYMBOLS = 10               # the range at the vertices is computed up to this many symbols
ACTIVE = 1e-9                         # buses whose first-order range is narrower are left out of the ratios [p.u.]


def pv_variant(g, s0, cap):
    """The grid with the sources (cap > 0) as PV buses, the centre injections and the sources: (grid, s_c, dg).

    The set-points are the voltage magnitudes of the AC solution of the original (PQ) grid at s_c = s0 + cap/2."""
    dg = np.flatnonzero(cap > 0)
    s_c = np.asarray(s0, dtype=complex) + cap / 2.0
    V, ok, _, _ = newton_pf(g.Ybus, s_c, g.V0, g.pv, g.pq, tol=TOL, max_iter=40)
    if not ok:
        raise RuntimeError("no AC solution at the centre of the uncertainty box")
    v_set = np.asarray(g.v_set, dtype=float).copy()
    v_set[dg] = np.abs(V[dg])
    pv = np.sort(np.concatenate([g.pv, dg])).astype(int)
    return replace(g, pv=pv, pq=np.setdiff1d(g.pq, dg), v_set=v_set, V0=V), s_c, dg


def linear_radius(so) -> float:
    """Spectral radius of G = sum_k |J^-1 dJ(E_k)|, the linear part of cdhim.enclosure.Enclosure (same formula)."""
    g, n2 = so.st.grid, 2 * so.m
    G = np.zeros((n2, n2))
    for k in range(so.E.shape[1]):
        G += np.abs(solve_blocked(so.lu, rect_jacobian(g, so.E[:, k]).toarray()))
    return perron(G)


def perron(A, iterations: int = 300) -> float:
    """Perron root of a non-negative matrix by power iteration (as in cdhim.enclosure.Enclosure)."""
    x = np.ones(A.shape[0])
    lam = 0.0
    for _ in range(iterations):
        y = A @ x
        lam = float(np.linalg.norm(y))
        if lam == 0.0:
            return 0.0
        x = y / lam
    return lam


def majorant_matrices(so, Hs=None):
    """The matrices of the radii map of cdhim.verified_vi.verify_tube_vi and its coefficients, recomputed here for
    diagnosis with the same formulas: dict(RV, TV, RI, TI, delta, Ibar, phi_low, is_pv). The directed rounding of
    the module is kept, but this is a copy made to look inside a verification that failed, not a proof of anything."""
    g, st = so.st.grid, so.st
    ns, m, npv = g.pvpq, so.m, len(g.pv)
    is_pv = np.r_[np.ones(npv, bool), np.zeros(m - npv, bool)]
    pqf = (~is_pv).astype(float)
    Hs = state_forms(so) if Hs is None else Hs
    E, B = so.E[ns], Hs[:m] + 1j * Hs[m:]
    K = E.shape[1]
    V = st.V_star
    Y = g.Ybus.tocsr()
    Yn = Y[ns].tocsr()
    YL = Yn[:, ns].tocsr()
    p = max(int(np.diff(Yn.indptr).max()), 1)
    aYn, aYL = _sparse_abs(Yn), _sparse_abs(YL)
    Vn = V[ns]
    vr, vi = Vn.real, Vn.imag
    cur = Yn @ V
    e_cur = up(err_const(p) * matmul_up(aYn, abs_up(V), p))
    aE, aB = abs_up(E), abs_up(B)
    YE = YL @ E
    e_YE = up(err_const(p) * matmul_up(aYL, aE, p))
    YB = (YL @ B.reshape(m, K * K)).reshape(m, K, K)
    e_YB = up(err_const(p) * matmul_up(aYL, aB.reshape(m, K * K), p)).reshape(m, K, K)
    alpha, beta = sum_up(aE, axis=1), _cform_sup(B)
    YEn = sum_up(add_up(abs_up(YE), e_YE), axis=1)
    YBn = add_up(_cform_sup(YB), sum_up(e_YB, axis=(1, 2)))
    delta, Ibar = add_up(alpha, beta), add_up(YEn, YBn, e_cur)
    high = add_up(mul_up(alpha, YBn), mul_up(beta, YEn), mul_up(beta, YBn))
    high = np.where(is_pv, add_up(high, 2.0 * mul_up(alpha, beta), mul_up(beta, beta)), high)

    Yr = sp.bmat([[YL.real, -YL.imag], [YL.imag, YL.real]], format="csr")
    d11, d12, d21, d22 = vr, vi, -vi * pqf, vr * pqf
    ir, ii = cur.real, cur.imag
    k11, k12, k21, k22 = ir, ii, np.where(is_pv, 2.0 * vr, ii), np.where(is_pv, 2.0 * vi, -ir)
    Dv = sp.bmat([[sp.diags(d11), sp.diags(d12)], [sp.diags(d21), sp.diags(d22)]], format="csr")
    Kc = sp.bmat([[sp.diags(k11), sp.diags(k12)], [sp.diags(k21), sp.diags(k22)]], format="csr")
    C = np.linalg.inv((Dv @ Yr + Kc).toarray())
    aC, aYr = np.abs(C), abs(Yr)
    pc = max(int(np.diff(Yr.tocsc().indptr).max()), 1)
    pr = max(int(np.diff(Yr.indptr).max()), 1)

    def two_terms(A, c1, c2, c3, c4):
        a1, a2 = A[:, :m], A[:, m:]
        b1, b2 = np.abs(a1), np.abs(a2)
        return (np.hstack([a1 * c1 + a2 * c3, a1 * c2 + a2 * c4]),
                up(err_const(2) * np.hstack([add_up(mul_up(b1, np.abs(c1)), mul_up(b2, np.abs(c3))), add_up(mul_up(b1, np.abs(c2)), mul_up(b2, np.abs(c4)))])))

    CD, e_CD = two_terms(C, d11, d12, d21, d22)
    P1 = (Yr.T.tocsr() @ CD.T).T
    e_P1 = add_up(up(err_const(pc) * matmul_up(np.abs(CD), aYr, pc)), matmul_up(e_CD, aYr, pc))
    P2, e_P2 = two_terms(C, k11, k12, k21, k22)
    eye = np.eye(2 * m)
    res = add_up(np.abs(eye - P1 - P2), e_P1, e_P2, up(err_const(3) * add_up(eye, np.abs(P1), np.abs(P2))))
    return dict(RV=block_frobenius_up(res, m), TV=block_spectral_up(C, m),
                TI=add_up(block_spectral_up(Yr @ C, m), block_frobenius_up(up(err_const(pr) * matmul_up(aYr, aC, pr)), m)),
                RI=block_frobenius_up(matmul_up(aYr, res, pr), m), delta=delta, Ibar=Ibar, high=high, is_pv=is_pv)


def majorant_radii(so, Hs=None) -> dict:
    """Spectral radii of the linear part of the radii map (see the module docstring); plain floating point."""
    q = majorant_matrices(so, Hs)
    d1 = q["Ibar"] + np.where(q["is_pv"], 2.0 * q["delta"], 0.0)
    d2 = q["delta"]
    m = len(d2)
    zero = np.zeros((m, m))
    L0 = np.block([[q["TV"] * d1, q["TV"] * d2], [q["TI"] * d1, q["TI"] * d2]])
    L = L0 + np.block([[q["RV"], zero], [q["RI"], zero]])
    return {"majorant_radius": perron(L), "majorant_radius_without_rounding_terms": perron(L0),
            "majorant_radius_voltage_block": perron(q["TV"] * d1), "majorant_radius_current_block": perron(q["TI"] * d2),
            "rounding_block_RI_row_sum_max": float(q["RI"].sum(axis=1).max()), "rounding_block_RV_row_sum_max": float(q["RV"].sum(axis=1).max()),
            "TV_row_sum_max": float(q["TV"].sum(axis=1).max()), "TI_row_sum_max": float(q["TI"].sum(axis=1).max()),
            "Ibar_max_pu": float(q["Ibar"].max()), "model_deviation_max_pu": float(d2.max())}


def vertex_range(st, delta):
    """(vmin, vmax, vertices without AC solution) of |V| at the PQ buses over the vertices of the box."""
    K = st.src.n
    g = st.grid
    eps = np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** K)] for j in range(K)])
    vmin, vmax, fails = np.full(len(g.pq), np.inf), np.full(len(g.pq), -np.inf), 0
    for s in range(0, eps.shape[1], 256):
        V, conv, _ = st.solver.solve(st.injections(eps[:, s:s + 256], delta))
        fails += int((~conv).sum())
        if conv.any():
            vm = np.abs(V[g.pq][:, conv])
            vmin, vmax = np.minimum(vmin, vm.min(axis=1)), np.maximum(vmax, vm.max(axis=1))
    return vmin, vmax, fails


def measure(so, st, delta, base):
    """The verified enclosure in voltage-current form of one case, with its checks and diagnostics."""
    g = st.grid
    Hs = state_forms(so)
    t0 = perf_counter()
    enc = verified_enclosure_vi(so, node_budget=NODE_BUDGET, Hs=Hs)
    row = dict(base, pv_buses=len(g.pv), pq_buses=len(g.pq), symbols=st.src.n, verified=bool(enc.verified), unique=bool(enc.unique),
               reason=enc.reason, verified_ms=1e3 * (perf_counter() - t0), linear_radius_G=linear_radius(so))
    row.update(majorant_radii(so, Hs))
    if enc.verified:
        row["tube_radius_max_pu"] = float(enc.s_abs.max())
    if st.src.n <= MAX_VERTEX_SYMBOLS:
        vmin, vmax, fails = vertex_range(st, delta)
        act = np.abs(so.c_lin).sum(axis=1) > ACTIVE
        row.update({"vertices_without_ac_solution": fails, "buses": int(act.sum()), "vertex_range_width_max_pu": float((vmax - vmin)[act].max())})
        if enc.verified:
            ratio = ((enc.hi - enc.lo) / (vmax - vmin))[act]
            row.update({"width_over_vertex_range_median": float(np.median(ratio)), "width_over_vertex_range_max": float(ratio.max()),
                        "violations": int(((enc.hi < vmax - VIOLATION_ATOL) | (enc.lo > vmin + VIOLATION_ATOL)).sum())})
    return row


def distribution():
    rows = []
    for name, case, prepare, _ in (NETWORKS[0], NETWORKS[4], NETWORKS[5]):
        g = load_grid(case, tol=TOL, prepare=prepare)
        rating = benchmark_ratings(g)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                base = {"network": name, "scenario": scen, "rho_U": rho_U}
                try:
                    gpv, s_c, dg = pv_variant(g, s0, cap)
                    st = custom_study(gpv, s_c, dg, np.zeros(len(dg), bool), cap[dg] / 2.0, tol=TOL)
                except RuntimeError as e:
                    rows.append(dict(base, status=str(e)))
                    continue
                row = measure(SecondOrder(st, 1.0), st, 1.0, dict(base, status="ok"))
                row["set_point_min_pu"], row["set_point_max_pu"] = float(gpv.v_set[dg].min()), float(gpv.v_set[dg].max())
                rows.append(row)
                report("distribution", row)
    return rows


def transmission():
    rows = []
    for case in TRANSMISSION:
        g = grid(case)
        st = study(case, zones=spectral_zones(g, ZONES, laplacian_embedding(g)))
        for delta in DELTAS:
            row = measure(SecondOrder(st, delta), st, delta, {"network": case, "delta": delta, "status": "ok",
                                                              "largest_admittance_pu": float(abs(g.Ybus).max())})
            rows.append(row)
            report("transmission", row)
    return rows


def forms():
    rows = []
    for name, case, prepare, _ in NETWORKS:
        g = load_grid(case, tol=TOL, prepare=prepare)
        rating = benchmark_ratings(g)
        for scen in "FR":
            sc = SCENARIOS[scen]
            for rho_U in RHO_U:
                s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
                row = {"network": name, "scenario": scen, "rho_U": rho_U, "n_bus": g.n_bus, "symbols": 2 * int((cap > 0).sum()),
                       "largest_admittance_pu": float(abs(g.Ybus).max())}
                try:
                    st = half_disk_study(g, s0, cap, tol=TOL)
                except RuntimeError as e:
                    rows.append(dict(row, status=str(e)))
                    continue
                so = SecondOrder(st, 1.0)
                Hs = state_forms(so)
                t0 = perf_counter()
                a = verified_enclosure(so, half_disks=True, node_budget=NODE_BUDGET, Hs=Hs)
                t1 = perf_counter()
                b = verified_enclosure_vi(so, half_disks=True, node_budget=NODE_BUDGET, Hs=Hs)
                t2 = perf_counter()
                row.update({"status": "ok", "current_form_verified": bool(a.verified), "current_form_unique": bool(a.unique), "current_form_reason": a.reason,
                            "vi_form_verified": bool(b.verified), "vi_form_unique": bool(b.unique), "vi_form_reason": b.reason,
                            "current_form_ms": 1e3 * (t1 - t0), "vi_form_ms": 1e3 * (t2 - t1), "linear_radius_G": linear_radius(so)})
                row.update(majorant_radii(so, Hs))
                if a.verified and b.verified:
                    r = b.s_abs / a.s_abs
                    w = (b.hi - b.lo) / (a.hi - a.lo)
                    row.update({"vi_over_current_radius_median": float(np.median(r)), "vi_over_current_radius_max": float(r.max()),
                                "vi_over_current_width_min": float(w.min()), "vi_over_current_width_max": float(w.max())})
                rows.append(row)
                print(f"[x14] forms {name} {scen} rho_U={rho_U:g}: current {a.verified}/{a.unique} VI {b.verified}/{b.unique} {b.reason}"
                      f" | majorant radius {row['majorant_radius']:.3f} (without rounding terms {row['majorant_radius_without_rounding_terms']:.3f})"
                      + (f" | radius VI/current {row['vi_over_current_radius_median']:.4f}" if a.verified and b.verified else ""), flush=True)
    return rows


def report(part, row):
    print(f"[x14] {part} {row['network']} {row.get('scenario', '')} {row.get('rho_U', row.get('delta'))}: PV {row['pv_buses']} PQ {row['pq_buses']} K {row['symbols']}"
          f" verified {row['verified']} unique {row['unique']} {row['reason']} | G radius {row['linear_radius_G']:.3f} majorant {row['majorant_radius']:.3f}"
          f" (voltage block {row['majorant_radius_voltage_block']:.3f}, current block {row['majorant_radius_current_block']:.3f})"
          + (f" | width/vertex range {row['width_over_vertex_range_median']:.4f} (max {row['width_over_vertex_range_max']:.4f}) violations {row['violations']}"
             if "width_over_vertex_range_median" in row else ""), flush=True)


PARTS = {"distribution": ("x14_pv_distribution", distribution), "transmission": ("x14_pv_transmission", transmission), "forms": ("x14_forms", forms)}


def main():
    meta = {"script": "experiments/x14_pv.py", "pf_tol": TOL, "rho_U": RHO_U, "zones": ZONES, "deltas": list(DELTAS), "raak_node_budget": NODE_BUDGET,
            "transmission_cases": list(TRANSMISSION), "networks": {n: d for n, _, _, d in NETWORKS},
            "note": "linear_radius_G and the majorant radii are diagnostics in plain floating point; verified / unique come from cdhim.verified_vi "
                    "(and cdhim.verified in x14_forms) with directed rounding"}
    wanted = sys.argv[1:] or list(PARTS)
    unknown = [w for w in wanted if w not in PARTS]
    if unknown:
        sys.exit(f"unknown table(s) {unknown}; choose from {list(PARTS)}")
    for part in wanted:
        table, fn = PARTS[part]
        print(save_table(table, fn(), meta))


if __name__ == "__main__":
    main()
