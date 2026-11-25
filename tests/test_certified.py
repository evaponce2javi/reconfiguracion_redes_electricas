"""Certified brackets of the exact voltage range: inner points, bisection, and their consistency."""
import warnings

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid, make_study, spectral_zones
from cdhim.certified import box_bracket, certify_points, half_disk_eps, half_disk_inner, in_half_disks, injection_matrix, matrix_study
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import injection_changes
from cdhim.halfdisk import attained_points, half_disk_study
from cdhim.verified import verified_enclosure
from cdhim.zones import laplacian_embedding


def _close_ties(net):
    net.line["in_service"] = True


@pytest.fixture(scope="module")
def grid():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid("case33bw", prepare=_close_ties)


@pytest.fixture(scope="module")
def half_disks(grid):
    sc = SCENARIOS["R"]
    s0, cap = net_injections(grid, benchmark_ratings(grid), sc["rho_D"], sc["rho_G"], sc["pf"], 0.5)
    return s0, cap, half_disk_study(grid, s0, cap)


def test_injection_matrix_is_the_one_of_the_enclosure_and_round_trips(grid, half_disks):
    s0, cap, st = half_disks
    M = injection_matrix(st)
    assert np.array_equal(M, injection_changes(SecondOrder(st, 1.0)))
    st2 = matrix_study(grid, st.sbus, M)
    assert np.array_equal(injection_matrix(st2), M) and np.array_equal(st2.sbus, st.sbus)
    zs = make_study(grid, zones=spectral_zones(grid, 4, laplacian_embedding(grid)))
    assert np.array_equal(injection_matrix(zs, 0.2), injection_changes(SecondOrder(zs, 0.2)))


def test_half_disk_membership_is_exact():
    assert in_half_disks(np.array([1.0, 0.0])) and in_half_disks(np.array([-1.0, 1.0])) and in_half_disks(np.array([-1.0, -1.0]))
    assert not in_half_disks(np.array([1.0, 1e-300])) and not in_half_disks(np.array([-1.0, np.nextafter(1.0, 2.0)]))
    assert not in_half_disks(np.array([np.nextafter(-1.0, -2.0), 0.0]))
    ang = np.linspace(-np.pi / 2, np.pi / 2, 181)
    pts = (np.cos(ang) + 1j * np.sin(ang))[:, None]                   # on the arc: rounding may leave them outside
    eps = half_disk_eps(pts)
    assert all(in_half_disks(eps[:, k]) for k in range(eps.shape[1]))
    assert np.abs(eps - np.vstack([2 * pts.real.T - 1, pts.imag.T])).max() < 1e-11


def test_certified_points_bracket_the_computed_solution(grid, half_disks):
    s0, cap, st = half_disks
    M = injection_matrix(st)
    eps = np.random.default_rng(0).uniform(-1, 1, (M.shape[1], 5))
    lo, hi = certify_points(grid, st.sbus, M, eps)
    V, conv, _ = st.solver.solve(st.sbus[:, None] + M @ eps)
    assert conv.all()
    vm = np.abs(V[grid.pq])
    assert np.all(lo <= vm + 1e-9) and np.all(vm <= hi + 1e-9) and (hi - lo).max() < 1e-7


def test_inner_bounds_over_half_disks_agree_with_the_enclosure(grid, half_disks):
    s0, cap, st = half_disks
    vmax, pmax = attained_points(grid, s0, cap, +1)
    vmin, pmin = attained_points(grid, s0, cap, -1)
    vmin_hi, vmax_lo = half_disk_inner(st, pmax, pmin)
    assert np.abs(vmin_hi - vmin).max() < 1e-7 and np.abs(vmax_lo - vmax).max() < 1e-7
    enc = verified_enclosure(SecondOrder(st, 1.0), half_disks=True)
    assert enc.verified and np.all(enc.lo <= vmin_hi) and np.all(enc.hi >= vmax_lo)          # two independent proofs


def test_bisection_brackets_the_exact_range(grid):
    st = make_study(grid, zones=spectral_zones(grid, 4, laplacian_embedding(grid)))
    M = injection_matrix(st, 0.3)
    br = box_bracket(grid, st.sbus, M, tol=1e-3, max_boxes=400)
    assert br.unverified == 0 and br.gap <= 1e-3
    assert np.all(br.vmin_lo <= br.vmin_hi) and np.all(br.vmax_lo <= br.vmax_hi)
    assert np.all(br.root_lo <= br.vmin_lo) and np.all(br.root_hi >= br.vmax_hi)             # refining can only tighten
    # every AC solution of a fine sample of the box lies inside the outer ends, and none beats the inner ends by much
    rng = np.random.default_rng(1)
    eps = np.hstack([rng.uniform(-1, 1, (4, 2000)), np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(16)] for j in range(4)])])
    V, conv, _ = st.solver.solve(st.injections(eps, 0.3))
    vm = np.abs(V[grid.pq][:, conv])
    assert np.all(vm.min(axis=1) >= br.vmin_lo - 1e-9) and np.all(vm.max(axis=1) <= br.vmax_hi + 1e-9)
    width = br.vmax_lo - br.vmin_hi
    assert np.all((br.vmax_hi - br.vmin_lo) <= width * 1.003 + 1e-12)
    assert np.all(vm.max(axis=1) <= br.vmax_lo + 2e-3 * width + 1e-9)


def test_subdivision_with_discs_is_valid_and_needs_more_boxes(grid):
    st = make_study(grid, zones=spectral_zones(grid, 4, laplacian_embedding(grid)))
    M = injection_matrix(st, 0.3)
    ours = box_bracket(grid, st.sbus, M, tol=2e-2, max_boxes=300)
    discs = box_bracket(grid, st.sbus, M, tol=2e-2, max_boxes=300, leaf="discs")
    assert ours.gap <= 2e-2 and discs.unverified == 0
    assert discs.evaluated > ours.evaluated
    # the two outer ends bracket the same certified inner values, and the root disc is the wider of the two
    assert np.all(discs.vmin_lo <= ours.vmin_hi) and np.all(discs.vmax_hi >= ours.vmax_lo)
    assert np.median((discs.root_hi - discs.root_lo) / (ours.root_hi - ours.root_lo)) > 1.05
    V, conv, _ = st.solver.solve(st.injections(np.random.default_rng(2).uniform(-1, 1, (4, 500)), 0.3))
    vm = np.abs(V[grid.pq][:, conv])
    assert np.all(vm.min(axis=1) >= discs.vmin_lo - 1e-9) and np.all(vm.max(axis=1) <= discs.vmax_hi + 1e-9)
