"""Jacobian, sensitivities and AC solver against independent references (IEEE-33)."""
import warnings

import numpy as np
import pytest

from cdhim import load_grid, make_study, sample_eps
from cdhim.core import MethodConfig
from cdhim.powerflow import legacy_flat_jacobian, mismatch, polar_jacobian, stack_mismatch

# Values printed by the unmodified finalFinalCode.py (commit 4750ce2) on IEEE-33, delta = 15 %.
LEGACY_33 = dict(branches=31, ellipse=16, cri=86.70, top_branch=(5, 6), top_rho=0.9990, top_gain=96.45)


@pytest.fixture(scope="module")
def grid():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid("case33bw")


def finite_difference_jacobian(g, V, h=1e-6):
    va, vm = np.angle(V), np.abs(V)
    F = lambda a, m: stack_mismatch(mismatch(g.Ybus, m * np.exp(1j * a), 0.0), g.pvpq, g.pq)
    cols = []
    for b in g.pvpq:
        d = np.zeros(g.n_bus); d[b] = h
        cols.append((F(va + d, vm) - F(va - d, vm)) / (2 * h))
    for b in g.pq:
        d = np.zeros(g.n_bus); d[b] = h
        cols.append((F(va, vm + d) - F(va, vm - d)) / (2 * h))
    return np.array(cols).T


def test_exact_jacobian_matches_finite_differences(grid):
    J = polar_jacobian(grid.Ybus, grid.V0, grid.pvpq, grid.pq).toarray()
    np.testing.assert_allclose(J, finite_difference_jacobian(grid, grid.V0), atol=2e-6)


def test_legacy_jacobian_has_the_sign_error_in_H(grid):
    """Documents the defect of finalFinalCode.py:115 -- off-diagonal dP/dtheta is exactly negated."""
    v = grid.flat_voltage()
    n = len(grid.pvpq)
    legacy = legacy_flat_jacobian(grid.Ybus, np.abs(v), grid.pvpq, grid.pq).toarray()
    exact = finite_difference_jacobian(grid, v)
    off = ~np.eye(n, dtype=bool)
    np.testing.assert_allclose(legacy[:n, :n][off], -exact[:n, :n][off], atol=1e-6)        # H: wrong sign
    assert np.abs(exact[:n, :n][off]).max() > 1.0
    np.testing.assert_allclose(legacy[:n, n:], exact[:n, n:], atol=1e-6)                   # N, M, L agree
    np.testing.assert_allclose(legacy[n:], exact[n:], atol=1e-6)


def test_sensitivities_equal_dense_inverse_block(grid):
    st = make_study(grid)
    Jinv = np.linalg.inv(st.lin.J.toarray())
    ref = Jinv[len(grid.pvpq):][:, st.src.row] * st.src.unit
    np.testing.assert_allclose(st.C1, ref, rtol=1e-9, atol=1e-13)
    assert st.src.n == 64 and (~st.src.is_q).sum() == 32                  # separate P and Q symbols per load


def test_legacy_mode_reproduces_the_submitted_figures(grid):
    st = make_study(grid, point="flat", jacobian="legacy_signbug")
    res = st.run(0.15, MethodConfig(min_sigma_product=1e-12, min_radius_pu=0.0))
    assert res.n_branches == LEGACY_33["branches"]
    assert res.n_ellipse == LEGACY_33["ellipse"]
    assert round(100 * res.cri_area(), 2) == LEGACY_33["cri"]
    ok = res.valid & (np.abs(res.rho) < 0.999)
    k = np.flatnonzero(ok)[np.argmax(res.gain[ok])]                       # branch of the manuscript's Fig. 4
    assert (int(grid.pq[res.fi[k]]), int(grid.pq[res.ti[k]])) == LEGACY_33["top_branch"]
    assert round(float(res.rho[k]), 4) == LEGACY_33["top_rho"]
    assert round(100 * float(res.gain[k]), 2) == LEGACY_33["top_gain"]
    # gamma_th in {5, 10, 15, 20, 25} % selects the same branches: the threshold never binds
    for gamma in (0.05, 0.10, 0.20, 0.25):
        np.testing.assert_array_equal(st.run(0.15, MethodConfig(gamma_th=gamma, min_sigma_product=1e-12)).ellipse, res.ellipse)


def test_selection_is_independent_of_delta(grid):
    st = make_study(grid)
    ref = st.run(0.05)
    for delta in (0.10, 0.20, 0.30):
        res = st.run(delta)
        np.testing.assert_array_equal(res.ellipse, ref.ellipse)
        np.testing.assert_allclose(res.rho, ref.rho, rtol=1e-10)
        np.testing.assert_allclose(res.r_out, ref.r_out * delta / 0.05, rtol=1e-10)


def test_batch_solver_matches_pandapower(grid):
    """Every sample of the batched solver is the same AC solution pandapower's Newton-Raphson returns."""
    import pandapower as pp
    import pandapower.networks as pn

    st = make_study(grid)
    delta = 0.30
    eps = sample_eps(np.random.default_rng(3), st.src.n, 12)
    V, conv, stats = st.solver.solve(st.injections(eps, delta))
    assert conv.all() and stats.converged == 12
    res = np.abs(stack_mismatch(mismatch(grid.Ybus, V, st.injections(eps, delta)), grid.pvpq, grid.pq)).max()
    assert res < 1e-9
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net = pn.case33bw()
        p0, q0 = net.load.p_mw.values.copy(), net.load.q_mvar.values.copy()
        load_bus = net.load.bus.values
        pos = {int(b): k for k, b in enumerate(grid.bus_label)}
        ep = {int(grid.bus_label[b]): k for k, b in enumerate(st.src.bus) if not st.src.is_q[k]}
        eq = {int(grid.bus_label[b]): k for k, b in enumerate(st.src.bus) if st.src.is_q[k]}
        for j in range(12):
            net.load["p_mw"] = p0 * (1 + delta * np.array([eps[ep[int(b)], j] for b in load_bus]))
            net.load["q_mvar"] = q0 * (1 + delta * np.array([eps[eq[int(b)], j] for b in load_bus]))
            pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, tolerance_mva=1e-9)
            vm = np.array([net.res_bus.vm_pu.at[int(b)] for b in grid.bus_label])
            np.testing.assert_allclose(np.abs(V[:, j]), vm, atol=1e-8)


def test_linear_prediction_error_is_second_order(grid):
    """Halving delta must cut |dV_AC - dV_lin| by about four when J is evaluated at the AC solution."""
    st = make_study(grid)
    eps = sample_eps(np.random.default_rng(4), st.src.n, 200)
    errs = []
    for delta in (0.20, 0.10, 0.05):
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        assert conv.all()
        dv_ac = np.abs(V[grid.pq]) - st.vm_star[:, None]
        errs.append(np.abs(dv_ac - st.C1 @ eps * delta).max())
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.15)
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.15)


def test_load_factor_changes_the_operating_point(grid):
    st = make_study(grid, load_factor=1.4)
    res = np.abs(stack_mismatch(mismatch(grid.Ybus, st.V_star, grid.sbus(1.4)), grid.pvpq, grid.pq)).max()
    assert res < 1e-9
    assert np.abs(st.V_star).min() < np.abs(grid.V0).min()
    np.testing.assert_allclose(st.src.unit, 1.4 * make_study(grid).src.unit)
