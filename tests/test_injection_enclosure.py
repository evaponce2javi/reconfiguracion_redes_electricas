"""The current-injection enclosure against exact AC solutions: every piece of its derivation, then the result."""
import warnings

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import Enclosure, InjectionEnclosure, injection_changes, voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks


def _close_ties(net):
    net.line["in_service"] = True


def _grid(case, prepare=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case, prepare=prepare)


def _add_shunts(net):
    """Line charging on every line and two capacitor banks: the DistFlow model of the baseline has neither."""
    import pandapower as pp

    net.line["c_nf_per_km"] = 300.0
    for bus in (17, 32):
        pp.create_shunt(net, bus, q_mvar=-0.3)


# a radial feeder, the same feeder with its five tie lines closed, and the same feeder with shunt elements
CASES = [("case33bw", None, "F", 0.5), ("case33bw", _close_ties, "R", 0.8), ("case33bw", _add_shunts, "F", 0.5)]


@pytest.fixture(scope="module", params=CASES, ids=["radial", "meshed", "shunts"])
def model(request):
    case, prepare, scen, rho_U = request.param
    g = _grid(case, prepare)
    sc = SCENARIOS[scen]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
    st = half_disk_study(g, s0, cap)
    return g, s0, cap, st, SecondOrder(st, 1.0)


def vertices(k):
    return np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** k)] for j in range(k)])


def solve(st, eps):
    V, conv, _ = st.solver.solve(st.injections(eps, 1.0))
    assert conv.all()
    return V


def test_study_is_centred_on_the_bounding_boxes(model):
    g, s0, cap, st, so = model
    dg = np.flatnonzero(cap > 0)
    n = len(dg)
    assert st.src.n == 2 * n <= 8
    lo = st.injections(-np.ones((2 * n, 1)), 1.0)[:, 0]
    hi = st.injections(np.ones((2 * n, 1)), 1.0)[:, 0]
    np.testing.assert_allclose(lo, s0 - 1j * cap, atol=1e-15)                      # (dp, dq) = (0, -C)
    np.testing.assert_allclose(hi, s0 + cap * (1 + 1j), atol=1e-15)                # (C, C): corner of the box


def test_exact_identity_of_the_current_injection_form(model):
    """J D = dS + V*^2 (-tau d + sig d^2) / (1 + d) for exact AC solutions inside the box."""
    g, s0, cap, st, so = model
    ns, m = g.pvpq, so.m
    eps = np.random.default_rng(3).uniform(-1, 1, (st.src.n, 20))
    Vs = st.V_star[ns][:, None]
    D = solve(st, eps)[ns] - Vs
    d = D / Vs
    dS = injection_changes(so)[ns] @ eps
    rhs = dS + Vs ** 2 * (-(dS / Vs ** 2) * d + (st.sbus[ns][:, None] / Vs ** 2) * d * d) / (1 + d)
    x = so.lu.solve(np.vstack([rhs.real, rhs.imag]))
    np.testing.assert_allclose(x[:m] + 1j * x[m:], D, rtol=0, atol=1e-8)


def test_partial_bounds_and_remainder_hold_at_every_vertex(model):
    g, s0, cap, st, so = model
    enc = InjectionEnclosure(so)
    fp = enc.solve(1.0)
    assert fp.converged and fp.rel_radius < 1.0
    ns, m, k = g.pvpq, so.m, st.src.n
    eps = vertices(k)
    Vs = st.V_star[ns][:, None]
    D = solve(st, eps)[ns] - Vs
    D1 = so.E[ns] @ eps
    r2 = enc.Hs.reshape(2 * m, -1) @ (eps[:, None, :] * eps[None, :, :]).reshape(k * k, -1)
    D2 = r2[:m] + 1j * r2[m:]
    tau = injection_changes(so)[ns] @ eps / Vs ** 2
    tol = 1e-12
    assert np.all(np.abs(D1 / Vs).max(axis=1) <= enc.a1 + tol)
    assert np.all(np.abs(D2 / Vs).max(axis=1) <= enc.a2 + tol)
    assert np.all(np.abs(tau).max(axis=1) <= enc.t + tol)
    s_true = np.abs(D - D1 - D2).max(axis=1)
    assert np.all(s_true <= fp.s_abs[ns] + 1e-10)                                  # 1e-10: tolerance of the AC solver
    assert s_true.max() > 0.1 * fp.s_abs[ns].max()                                 # and the bound is not vacuous


def test_remainder_is_of_third_order(model):
    g, s0, cap, st, so = model
    ns, m, k = g.pvpq, so.m, st.src.n
    Hs = InjectionEnclosure(so).Hs.reshape(2 * m, -1)
    eps = np.random.default_rng(4).choice([-1.0, 1.0], (k, 6))
    out = []
    for h in (0.5, 0.25):
        e = h * eps
        D = solve(st, e)[ns] - st.V_star[ns][:, None]
        r2 = Hs @ (e[:, None, :] * e[None, :, :]).reshape(k * k, -1)
        out.append(np.abs(D - so.E[ns] @ e - (r2[:m] + 1j * r2[m:])).max())
    assert 6.0 < out[0] / out[1] < 10.5


def test_first_order_range_over_half_disks_matches_brute_force():
    rng = np.random.default_rng(5)
    C = rng.normal(size=(6, 8))
    lo, hi = first_order_range(C)
    ang = np.linspace(-np.pi / 2, np.pi / 2, 20001)
    u, w = np.cos(ang), np.sin(ang)
    for i in range(6):
        best_hi = best_lo = 0.0
        for n in range(4):
            f = -C[i, n] + 2 * C[i, n] * u + C[i, 4 + n] * w          # on the arc; the linear form peaks on the boundary
            g0 = -C[i, n] + C[i, 4 + n] * np.linspace(-1, 1, 2001)    # and on the flat side u = 0
            best_hi += max(f.max(), g0.max())
            best_lo += min(f.min(), g0.min())
        assert hi[i] >= best_hi - 1e-12 and hi[i] - best_hi < 1e-6
        assert lo[i] <= best_lo + 1e-12 and best_lo - lo[i] < 1e-6


def test_enclosure_contains_the_ac_voltages(model):
    g, s0, cap, st, so = model
    enc = voltage_enclosure(so, first_order_range(so.c_lin))
    box = voltage_enclosure(so)
    assert enc.converged and box.converged
    dg = np.flatnonzero(cap > 0)
    tol = 1e-9
    vm = np.abs(solve(st, vertices(st.src.n))[g.pq])                               # the box, corners included
    assert np.all(vm.min(axis=1) >= box.lo - tol) and np.all(vm.max(axis=1) <= box.hi + tol)
    pts = sample_half_disks(np.random.default_rng(6), len(dg), 2000)              # the half disks
    S = np.repeat(s0[:, None], 2000, axis=1)
    S[dg] += cap[dg][:, None] * pts
    V, conv, _ = st.solver.solve(S)
    assert conv.all()
    vm = np.abs(V[g.pq])
    assert np.all(vm.min(axis=1) >= enc.lo - tol) and np.all(vm.max(axis=1) <= enc.hi + tol)
    for sign, bound in ((1, enc.hi), (-1, enc.lo)):                                # and the extremes found by local search
        ext = attained_extreme(g, s0, cap, sign)
        assert np.all(sign * (bound - ext) >= -tol)
    assert np.all(enc.hi <= box.hi + 1e-12) and np.all(enc.lo >= box.lo - 1e-12)   # half disks are inside their boxes
    width = attained_extreme(g, s0, cap, 1) - attained_extreme(g, s0, cap, -1)
    assert np.median((enc.hi - enc.lo) / width) < 1.1                              # tight: within 10 % of the attained range


def test_closes_where_the_rectangular_bound_does_not():
    g = _grid("case33bw")
    sc = SCENARIOS["F"]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], 0.8)
    so = SecondOrder(half_disk_study(g, s0, cap), 1.0)
    assert InjectionEnclosure(so).solve(1.0).converged
    assert not Enclosure(so).solve(1.0).converged


def test_voltage_controlled_buses_are_rejected():
    from cdhim import make_study, spectral_zones

    g = _grid("case118")
    so = SecondOrder(make_study(g, zones=spectral_zones(g, 8)), 0.02)
    with pytest.raises(ValueError, match="PQ"):
        InjectionEnclosure(so)
