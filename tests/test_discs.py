"""The published disc bound for PQ networks (Wang et al. 2018): its printed example and its guarantee."""
import warnings

import numpy as np
import pytest

from cdhim import load_grid
from cdhim.discs import disc_bound, disc_radius
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.halfdisk import sample_half_disks
from cdhim.powerflow import BatchACSolver, newton_pf


def _close_ties(net):
    net.line["in_service"] = True


def test_radius_of_the_printed_example():
    """Sec. VI of the paper: u_min = 1.0050, xi(s^) = 0.5692, xi(s - s^) = 0.0164 give rho = 0.0412."""
    assert abs(disc_radius(1.0050, 0.5692, 0.0164) - 0.0412) < 1e-4
    assert disc_radius(1.0, 0.0, 0.25) is None and disc_radius(1.0, 1.0, 0.0) is None
    assert abs(disc_radius(1.0, 0.0, 0.09) - 0.1) < 1e-12                          # corollary: (1 - sqrt(1 - 4 xi)) / 2


@pytest.mark.parametrize("prepare", [None, _close_ties], ids=["radial", "meshed"])
def test_disc_contains_sampled_solutions(prepare):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = load_grid("case33bw", prepare=prepare)
    sc = SCENARIOS["F"]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], 0.2)
    v_hat, ok, _, _ = newton_pf(g.Ybus, s0, g.V0, g.pv, g.pq, tol=1e-10, max_iter=40)
    assert ok
    db = disc_bound(g, s0, v_hat, cap)                                             # the half disk of radius C lies in |ds| <= C
    assert db.ok and 0.0 < db.rho < db.u_min
    dg = np.flatnonzero(cap > 0)
    S = np.repeat(s0[:, None], 2000, axis=1)
    S[dg] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(7), len(dg), 2000)
    V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, v_hat).solve(S)
    assert conv.all()
    dev = np.abs(V[g.pq] - v_hat[g.pq][:, None]).max(axis=1)
    assert np.all(dev <= db.radius_bus + 1e-9) and np.all(db.radius_bus <= db.radius + 1e-15)
    assert dev.max() > 0.05 * db.radius_bus.max()                                  # not vacuous


def test_disc_bound_reports_failure_for_a_large_set():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = load_grid("case33bw")
    db = disc_bound(g, -g.s_load, g.V0, 50.0 * np.abs(g.s_load))
    assert not db.ok and np.isnan(db.radius).all()
