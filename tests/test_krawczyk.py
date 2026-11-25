"""The Krawczyk baseline: every piece of the operator against sampled points, then the enclosure against AC solutions."""
import itertools
import warnings

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.krawczyk import Krawczyk, modulus_range


def _close_ties(net):
    net.line["in_service"] = True


def _add_shunts(net):
    import pandapower as pp

    net.line["c_nf_per_km"] = 300.0
    for bus in (17, 32):
        pp.create_shunt(net, bus, q_mvar=-0.3)


def _grid(case, prepare=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case, prepare=prepare)


def _model(prepare, scen, rho_U):
    g = _grid("case33bw", prepare)
    sc = SCENARIOS[scen]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
    st = half_disk_study(g, s0, cap)
    return g, s0, cap, st, SecondOrder(st, 1.0)


# a radial feeder, the same feeder with its five tie lines closed, and the same feeder with shunt elements
CASES = [(None, "F", 0.5), (_close_ties, "R", 0.8), (_add_shunts, "F", 0.5)]
VARIANTS = list(itertools.product(("power", "current"), (False, True), (False, True)))      # formulation, slope, rotate
IDS = [f"{f}-{'slope' if s else 'derivative'}-{'rotated' if r else 'plain'}" for f, s, r in VARIANTS]


@pytest.fixture(scope="module", params=CASES, ids=["radial", "meshed", "shunts"])
def model(request):
    return _model(*request.param)


def vertices(k):
    return np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** k)] for j in range(k)])


def solve(st, eps):
    V, conv, _ = st.solver.solve(st.injections(eps, 1.0))
    assert conv.all()
    return V


def state(kr, V, g):
    """AC solutions (n_bus, batch) in the coordinates of the variant: (2m, batch)."""
    U = np.conj(kr.rot)[:, None] * V[g.pvpq]
    return np.vstack([U.real, U.imag])


def test_modulus_range_matches_brute_force():
    rng = np.random.default_rng(1)
    for _ in range(200):
        c, h = rng.uniform(-1, 1, 2), rng.uniform(0, 1, 2) ** 3
        lo, hi = c - h, c + h
        u = np.r_[np.linspace(lo[0], hi[0], 201), np.clip(0.0, lo[0], hi[0])]
        w = np.r_[np.linspace(lo[1], hi[1], 201), np.clip(0.0, lo[1], hi[1])]
        mod = np.hypot(u[:, None], w[None, :])
        mn, mx = modulus_range(lo[0], hi[0], lo[1], hi[1])
        assert abs(mn - mod.min()) < 1e-12 and abs(mx - mod.max()) < 1e-12


@pytest.mark.parametrize("variant", VARIANTS, ids=IDS)
def test_first_order_term_is_the_one_of_the_second_order_model(model, variant):
    g, s0, cap, st, so = model
    kr = Krawczyk(so, *variant)
    L = kr.L[:kr.m] + 1j * kr.L[kr.m:]
    np.testing.assert_allclose(L, np.conj(kr.rot)[:, None] * so.E[g.pvpq], rtol=0, atol=1e-9)
    assert kr.b0.max() < 1e-8                                    # the centre is the AC solution


@pytest.mark.parametrize("variant", VARIANTS, ids=IDS)
def test_jacobian_and_slope_are_those_of_the_equations(model, variant):
    g, s0, cap, st, so = model
    kr = Krawczyk(so, *variant)
    rng = np.random.default_rng(2)
    x = kr.x_star + rng.uniform(-0.02, 0.02, 2 * kr.m)
    eps = rng.uniform(-1, 1, kr.K)
    J = kr.jacobian(x, eps)
    h = 1e-6
    fd = np.empty_like(J)
    for j in range(2 * kr.m):
        e = np.zeros(2 * kr.m)
        e[j] = h
        fd[:, j] = (kr.residual(x + e, eps) - kr.residual(x - e, eps)) / (2 * h)
    np.testing.assert_allclose(J, fd, rtol=0, atol=1e-6 * np.abs(J).max())
    lhs = kr.residual(x, eps) - kr.residual(kr.x_star, eps)
    np.testing.assert_allclose(kr.slope_matrix(x, eps) @ (x - kr.x_star), lhs, rtol=0, atol=1e-11 * np.abs(J).max())


def _radius(so, variant):
    """A box to test the bounds on: the one the current-injection operator ends with (it closes in these cases)."""
    res = Krawczyk(so, "current", variant[1], variant[2]).solve()
    assert res.closes
    return res.radius


@pytest.mark.parametrize("variant", VARIANTS, ids=IDS)
def test_interval_matrix_contains_the_matrices_at_sampled_points(model, variant):
    g, s0, cap, st, so = model
    kr = Krawczyk(so, *variant)
    r = _radius(so, variant)
    lo, hi = kr.matrix_interval(r)
    tol = 1e-12 * np.abs(lo).max()
    rng = np.random.default_rng(3)
    used = np.zeros(lo.shape)
    for t in range(60):
        corner = t % 2 == 0
        dx = r * (rng.choice([-1.0, 1.0], len(r)) if corner else rng.uniform(-1, 1, len(r)))
        eps = rng.choice([-1.0, 1.0], kr.K) if corner else rng.uniform(-1, 1, kr.K)
        A = (kr.slope_matrix if kr.slope else kr.jacobian)(kr.x_star + dx, eps)
        assert np.all(A >= lo - tol) and np.all(A <= hi + tol)
        used = np.maximum(used, np.abs(A - 0.5 * (lo + hi)))
    wide = (hi - lo) > 0
    assert (used[wide] / (0.5 * (hi - lo)[wide])).max() > 0.5     # and the interval is not far wider than needed


@pytest.mark.parametrize("variant", VARIANTS, ids=IDS)
def test_nonlinear_bound_holds_at_sampled_points(model, variant):
    """Phi(r) bounds |(I - C F'_x(xi, eps)) dx| for xi in the box and any |dx| <= r (derivative), and
    |(I - C S(x, eps)) (x - x*)| for x in the box (slope), for every eps of its box."""
    g, s0, cap, st, so = model
    kr = Krawczyk(so, *variant)
    r = _radius(so, variant)
    phi, a = kr.nonlinear_bound(r)
    assert phi is not None and a < 1.0 and np.all(np.isfinite(phi))
    rng = np.random.default_rng(4)
    worst = np.zeros(len(r))
    eye = np.eye(len(r))
    matrix = kr.slope_matrix if kr.slope else kr.jacobian
    for t in range(60):
        corner = t % 2 == 0
        dx = r * (rng.choice([-1.0, 1.0], len(r)) if corner else rng.uniform(-1, 1, len(r)))
        eps = rng.choice([-1.0, 1.0], kr.K) if corner else rng.uniform(-1, 1, kr.K)
        A = eye - kr.C @ matrix(kr.x_star + dx, eps)
        aligned = [r * np.sign(A[i]) for i in rng.integers(len(r), size=6)]          # vertices that are bad for one row
        if kr.slope:                                             # the matrix and the vector come from the same point
            pairs = [(A, dx)] + [(eye - kr.C @ matrix(kr.x_star + d, eps), d) for d in aligned]
        else:                                                    # any matrix of the box with any vector of the box
            pairs = [(A, d) for d in [dx, r * rng.choice([-1.0, 1.0], len(r))] + aligned]
        for B, d in pairs:
            val = np.abs(B @ d)
            assert np.all(val <= phi * (1 + 1e-12) + 1e-15)
            worst = np.maximum(worst, val)
    assert (worst / phi).max() > 0.75                             # and sampled values come close to the bound


@pytest.mark.parametrize("slope", [False, True], ids=["derivative", "slope"])
def test_power_formulation_interval_matrix_is_attained(model, slope):
    """The matrix is affine in x, so every entry of the interval is reached at a vertex of the box."""
    g, s0, cap, st, so = model
    kr = Krawczyk(so, "power", slope, True)
    r = _radius(so, ("power", slope, True))
    lo, hi = kr.matrix_interval(r)
    rng = np.random.default_rng(7)
    rows, cols = np.nonzero(hi - lo > 1e-9 * np.abs(hi).max())
    for k in rng.choice(len(rows), 25, replace=False):
        a, b = rows[k], cols[k]
        sgn = np.array([np.sign(A[a, b]) for A in kr._unit])
        sgn[sgn == 0] = 1.0
        M = (kr.slope_matrix if slope else kr.jacobian)(kr.x_star + r * sgn, np.zeros(kr.K))
        assert abs(M[a, b] - hi[a, b]) <= 1e-10 * max(1.0, abs(hi[a, b]))


@pytest.mark.parametrize("slope", [False, True], ids=["derivative", "slope"])
def test_current_formulation_bounds_of_the_bus_terms_are_sharp(model, slope):
    """Bounds of Re and Im of kap - kap* against its values at the vertices of (dU, eps): never exceeded, nearly reached."""
    g, s0, cap, st, so = model
    kr = Krawczyk(so, "current", slope, True)
    r = _radius(so, ("current", slope, True))
    m = kr.m
    lo, hi = kr.matrix_interval(r)
    rad = 0.5 * (hi - lo)
    k_r, k_i = np.diag(rad[:m, :m]), np.diag(rad[:m, m:])
    live = np.flatnonzero(k_r + k_i > 0)
    assert len(live) > 0
    best = np.zeros((2, m))
    tau_max = np.zeros(m)
    for eps in vertices(kr.K).T:
        S = kr.s_star + kr.M @ eps
        tau_max = np.maximum(tau_max, np.abs(np.conj(kr.M @ eps) / np.conj(kr.u_star) ** 2))
        for su, sw in itertools.product((-1.0, 1.0), repeat=2):
            U = kr.u_star + su * r[:m] + 1j * sw * r[m:]
            kap = -np.conj(S) / (np.conj(U) * np.conj(kr.u_star) if slope else np.conj(U) ** 2)
            d = kap + np.conj(kr.s_star) / np.conj(kr.u_star) ** 2
            vals = np.abs(np.vstack([d.real, d.imag]))
            assert np.all(vals <= np.vstack([k_r, k_i]) * (1 + 1e-12) + 1e-15)
            best = np.maximum(best, vals)
    assert np.all(best[0, live] >= 0.85 * k_r[live]) and np.all(best[1, live] >= 0.85 * k_i[live])
    np.testing.assert_allclose(kr.t, tau_max, rtol=1e-12, atol=0)             # t is the exact maximum of |tau| over the box


@pytest.mark.parametrize("variant", [v for v in VARIANTS if v[0] == "current"], ids=[i for i, v in zip(IDS, VARIANTS) if v[0] == "current"])
def test_enclosure_contains_the_ac_solutions(model, variant):
    g, s0, cap, st, so = model
    kr = Krawczyk(so, *variant)
    box = kr.solve()
    hd = kr.solve(first_order_range(kr.L))
    assert box.closes and hd.closes and box.rel_radius < 1.0
    assert np.all(kr.b0 + np.abs(kr.L).sum(axis=1) + kr.nonlinear_bound(box.radius)[0] <= box.radius * (1 + 1e-9))   # K(X) inside X
    tol = 1e-9
    m = kr.m
    eps = np.hstack([vertices(st.src.n), np.random.default_rng(5).uniform(-1, 1, (st.src.n, 500))])     # the box, corners included
    V = solve(st, eps)
    x = state(kr, V, g)
    assert np.all(x >= box.lo[:, None] - tol) and np.all(x <= box.hi[:, None] + tol)
    assert np.all(np.abs(x - kr.x_star[:, None] - kr.L @ eps) <= box.remainder[:, None] + tol)          # the affine form
    vm = np.abs(V[g.pq])
    assert np.all(vm.min(axis=1) >= box.vm_lo - tol) and np.all(vm.max(axis=1) <= box.vm_hi + tol)
    dg = np.flatnonzero(cap > 0)                                                                       # the half disks
    S = np.repeat(s0[:, None], 2000, axis=1)
    S[dg] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(6), len(dg), 2000)
    V, conv, _ = st.solver.solve(S)
    assert conv.all()
    x = state(kr, V, g)
    assert np.all(x >= hd.lo[:, None] - tol) and np.all(x <= hd.hi[:, None] + tol)
    vm = np.abs(V[g.pq])
    assert np.all(vm.min(axis=1) >= hd.vm_lo - tol) and np.all(vm.max(axis=1) <= hd.vm_hi + tol)
    vmax, vmin = attained_extreme(g, s0, cap, 1), attained_extreme(g, s0, cap, -1)
    assert np.all(vmax <= hd.vm_hi + tol) and np.all(vmin >= hd.vm_lo - tol)
    assert np.all(hd.vm_hi <= box.vm_hi + 1e-12) and np.all(hd.vm_lo >= box.vm_lo - 1e-12)             # half disks are inside their boxes
    assert np.median((hd.vm_hi - hd.vm_lo) / (vmax - vmin)) < 1.3                                       # and the bound is not vacuous


def test_power_formulation_closes_only_on_very_small_boxes():
    """Its Jacobian carries the admittances, so the box of independent voltages makes it lose the contraction early."""
    g, s0, cap, st, so = _model(None, "F", 0.001)
    for slope, rotate in itertools.product((False, True), (False, True)):
        kr = Krawczyk(so, "power", slope, rotate)
        res = kr.solve()
        assert res.closes
        V = solve(st, vertices(st.src.n))
        x = state(kr, V, g)
        assert np.all(x >= res.lo[:, None] - 1e-10) and np.all(x <= res.hi[:, None] + 1e-10)
        vm = np.abs(V[g.pq])
        assert np.all(vm.min(axis=1) >= res.vm_lo - 1e-10) and np.all(vm.max(axis=1) <= res.vm_hi + 1e-10)
    g, s0, cap, st, so = _model(None, "F", 0.05)
    assert not Krawczyk(so, "power", True, True).solve().closes
    assert Krawczyk(so, "current", False, False).solve().closes


def test_box_and_half_disk_modes_of_the_enclosure_differ_by_the_first_order_term(model):
    """experiments/x9_krawczyk.py gets the box mode of the enclosure from the half-disk one with this identity."""
    g, s0, cap, st, so = model
    hd = voltage_enclosure(so, first_order_range(so.c_lin))
    box = voltage_enclosure(so)
    lo1, hi1 = first_order_range(so.c_lin)
    r = np.abs(so.c_lin).sum(axis=1)
    np.testing.assert_allclose(box.lo, hd.lo - lo1 - r, rtol=0, atol=1e-13)
    np.testing.assert_allclose(box.hi, hd.hi - hi1 + r, rtol=0, atol=1e-13)


def test_voltage_controlled_buses_are_rejected():
    from cdhim import make_study, spectral_zones

    g = _grid("case118")
    so = SecondOrder(make_study(g, zones=spectral_zones(g, 8)), 0.02)
    with pytest.raises(ValueError, match="PQ"):
        Krawczyk(so)
