"""The radial baseline (Sou, Khong, Sandberg 2026): geometry helpers, its published Table I, and its guarantees."""
import warnings

import numpy as np
import pytest

from cdhim import load_grid
from cdhim.distflow import (SCENARIOS, _half_disk_distance, _half_disk_reach, benchmark_ratings, feeder, net_injections,
                            worst_case_load_flow)
from cdhim.halfdisk import sample_half_disks
from cdhim.network import apply_matpower_units, matpower_unit_statements
from cdhim.powerflow import BatchACSolver


def _grid(case, prepare=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case, prepare=prepare)


def _close_ties(net):
    net.line["in_service"] = True


@pytest.fixture(scope="module")
def grid():
    return _grid("case33bw")


def _half_disk_points(rad, n=361):
    """Points covering the half disk x >= 0, |x + jy| <= rad: a fan of radii times a fan of angles."""
    ang = np.linspace(-np.pi / 2, np.pi / 2, n)
    r = np.linspace(0.0, rad, 81)[:, None]
    return (r * np.cos(ang)).ravel(), (r * np.sin(ang)).ravel()


def test_distance_between_rectangle_and_half_disk_matches_brute_force():
    rng = np.random.default_rng(0)
    for _ in range(40):
        rad = rng.uniform(0.2, 2.0)
        a_lo, b_lo = rng.uniform(-3, 3, 2)
        a_hi, b_hi = a_lo + rng.uniform(0, 2), b_lo + rng.uniform(0, 2)
        hx, hy = _half_disk_points(rad)
        ra, rb = np.meshgrid(np.linspace(a_lo, a_hi, 33), np.linspace(b_lo, b_hi, 33))
        brute = np.hypot(hx[None, :] - ra.reshape(-1, 1), hy[None, :] - rb.reshape(-1, 1)).min()
        got = float(_half_disk_distance(np.array(a_lo), np.array(a_hi), np.array(b_lo), np.array(b_hi), rad))
        assert got <= brute + 1e-12                                   # the grid can only over-estimate a distance
        assert brute - got < 0.08 * (1.0 + rad)                       # and by no more than its resolution


def test_reach_of_half_disk_matches_brute_force():
    rng = np.random.default_rng(1)
    ang = np.linspace(-np.pi / 2, np.pi / 2, 20001)
    for _ in range(200):
        rad = rng.uniform(0.1, 2.0)
        a, b = rng.uniform(-3, 3, 2)
        brute = np.hypot(a + rad * np.cos(ang), b + rad * np.sin(ang)).max()      # a convex function peaks on the arc
        got = float(_half_disk_reach(np.array(a), np.array(b), rad))
        assert got >= brute - 1e-12 and got - brute < 1e-6 * (1.0 + got)


def test_feeder_matrices_reproduce_the_load_flow(grid):
    """nu = nu0 + 2 R p + 2 X q + M l (their Eq. 2a) must hold for the exact AC solution of the base case."""
    fd = feeder(grid)
    V = grid.V0
    p, q = -grid.s_load[fd.nodes].real, -grid.s_load[fd.nodes].imag
    ell = np.abs((V[fd.parent] - V[fd.nodes]) / (fd.r + 1j * fd.x)) ** 2
    nu = fd.nu0 + 2 * fd.R @ p + 2 * fd.X @ q + fd.M @ ell
    np.testing.assert_allclose(nu, np.abs(V[fd.nodes]) ** 2, rtol=0, atol=1e-9)
    assert fd.M.max() <= 1e-15


def test_feeder_rejects_meshed_and_voltage_controlled_networks():
    with pytest.raises(ValueError, match="not radial"):
        feeder(_grid("case33bw", prepare=_close_ties))
    with pytest.raises(ValueError, match="voltage-controlled"):
        feeder(_grid("case118"))


def test_published_table_1_is_reproduced(grid):
    """33-bus feeder, five sources of Fig. 10, rho_D = 1, rho_G = 0.5, pf = 0.8: eps_ave and eps_max of (vU, vL)."""
    lab = np.asarray(grid.bus_label)
    rating = np.zeros(grid.n_bus)
    for bus, kva in ((9, 400), (17, 1000), (21, 500), (24, 800), (32, 1200)):     # paper bus n = bus label n (0 = substation)
        rating[np.flatnonzero(lab == bus)[0]] = kva / 1e3 / grid.base_mva
    assert abs(rating.sum() * grid.base_mva - 3.9) < 1e-12
    for rho_U, ave, top in ((0.5, "0.086", "0.17"), (1.0, "0.18", "0.34")):       # as printed in the paper
        s0, cap = net_injections(grid, rating, 1.0, 0.5, 0.8, rho_U)
        wc = worst_case_load_flow(grid, s0, cap, current_start_amp=5000.0)
        e = 100 * (wc.v_max_U - wc.v_max_L) / wc.v_max_L
        assert wc.converged
        assert f"{e.mean():.2g}" == ave and f"{e.max():.2g}" == top
    s0, cap = net_injections(grid, rating, 1.0, 0.5, 0.8, 0.0)                    # no uncertainty: the bounds must meet
    wc = worst_case_load_flow(grid, s0, cap, current_start_amp=5000.0)
    assert np.abs(wc.v_max_U - wc.v_max_L).max() < 1e-8 and np.abs(wc.v_min - wc.v_max_L).max() < 1e-12


@pytest.mark.parametrize("scen", ["F", "R"])
def test_bounds_hold_on_sampled_injections(grid, scen):
    sc = SCENARIOS[scen]
    rating = benchmark_ratings(grid)
    s0, cap = net_injections(grid, rating, sc["rho_D"], sc["rho_G"], sc["pf"], sc["rho_U"])
    wc = worst_case_load_flow(grid, s0, cap)
    dg = np.flatnonzero(cap > 0)
    assert len(dg) == 4 and wc.converged
    pts = sample_half_disks(np.random.default_rng(2), len(dg), 3000)
    S = np.repeat(s0[:, None], 3000, axis=1)
    S[dg] += cap[dg][:, None] * pts
    V, conv, _ = BatchACSolver(grid.Ybus, grid.pv, grid.pq, grid.V0).solve(S)
    assert conv.all()
    vm = np.abs(V[wc.nodes])
    assert np.all(vm.min(axis=1) >= wc.v_min - 1e-9)                              # minimum at (0, -C)
    assert np.all(vm.max(axis=1) <= wc.v_max_U + 1e-9)
    assert np.all(wc.v_max_L <= wc.v_max_U + 1e-12)


def test_result_does_not_depend_on_the_initial_current_bound(grid):
    sc = SCENARIOS["F"]
    s0, cap = net_injections(grid, benchmark_ratings(grid), sc["rho_D"], sc["rho_G"], sc["pf"], sc["rho_U"])
    ref = worst_case_load_flow(grid, s0, cap, current_start=3.0)
    assert ref.converged
    for kw in (dict(current_start=1.0), dict(current_start=1e3), dict(current_start_amp=5000.0), dict(current_start_amp=1e9)):
        other = worst_case_load_flow(grid, s0, cap, **kw)
        assert other.converged
        np.testing.assert_allclose(other.v_max_U, ref.v_max_U, rtol=0, atol=1e-8)


M_FILE = """
mpc.baseMVA = 10;
mpc.bus = [ 1 3 0 0 0 0 1 1 0 12.66 1 1 1; ];
Vbase = mpc.bus(1, BASE_KV) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);
mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;
"""


def test_unit_statements_of_matpower_feeders_are_applied():
    bus = np.zeros((2, 13))
    bus[:, 9] = 12.66
    bus[1, 2:4] = (100.0, 60.0)                                       # kW, kVAr
    branch = np.zeros((1, 13))
    branch[0, 2:4] = (0.0922, 0.0470)                                 # ohm
    ppc = {"baseMVA": 10.0, "bus": bus, "branch": branch}
    assert len(apply_matpower_units(ppc, M_FILE)) == 2
    z_base = 12.66e3 ** 2 / 10e6
    np.testing.assert_allclose(ppc["branch"][0, 2:4], np.array([0.0922, 0.0470]) / z_base)
    np.testing.assert_allclose(ppc["bus"][1, 2:4], [0.1, 0.06])
    assert matpower_unit_statements("mpc.bus = [ 1 2 3; ];\nmpc.gen = [ 1 0 0; ];") == []          # plain case file: nothing to do
    with pytest.raises(ValueError, match="does not know"):
        matpower_unit_statements(M_FILE + "mpc.bus(:, VM) = 1.05;\n")
