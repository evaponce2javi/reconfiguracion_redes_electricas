"""The second-order model against the exact AC power flow (independent of how H_i is derived)."""
import warnings

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid, make_study, participation_factors, sample_eps, spectral_zones


def _grid(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case)


# (case, zones, generators share the imbalance, delta): a radial feeder without PV buses and a meshed
# transmission system with PV buses, with and without the generator response
CASES = [("case33bw", 8, False, 0.15), ("case118", 8, False, 0.05), ("case118", 8, True, 0.05)]


@pytest.fixture(scope="module", params=CASES, ids=lambda c: f"{c[0]}-K{c[1]}-{'participation' if c[2] else 'slack'}")
def model(request):
    case, k, share, delta = request.param
    g = _grid(case)
    st = make_study(g, zones=spectral_zones(g, k), participation=participation_factors(g) if share else None)
    return st, SecondOrder(st, delta), delta


def ac_dv(st, eps, delta):
    V, conv, _ = st.solver.solve(st.injections(eps, delta))
    assert conv.all()
    return np.abs(V[st.grid.pq]) - st.vm_star[:, None]


def test_rectangular_and_polar_sensitivities_agree(model):
    st, so, delta = model
    C = st.coefficients(delta)
    np.testing.assert_allclose(so.c_lin, C, rtol=1e-7, atol=1e-10 * np.abs(C).max())


def test_assembled_forms_equal_the_term_by_term_expression(model):
    st, so, delta = model
    H = so.forms()
    eps = sample_eps(np.random.default_rng(0), st.src.n, 6)
    direct = np.einsum("ks,ikl,ls->is", eps, H, eps)
    np.testing.assert_allclose(direct, so.quadratic(eps), rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(H, H.transpose(0, 2, 1))
    np.testing.assert_allclose(np.einsum("ikk->i", H), so.trace(), rtol=1e-9, atol=1e-15)


def test_form_is_the_second_derivative_of_the_ac_solution(model):
    """d|V|(h d) + d|V|(-h d) = 2 h^2 d'H d + O(h^4): H_i must be the Hessian of the exact solution map."""
    st, so, delta = model
    rng = np.random.default_rng(1)
    d = rng.uniform(-1, 1, (st.src.n, 5))
    exact = so.quadratic(d)                                   # d'H d for every bus and direction
    big = np.abs(exact) > 1e-3 * np.abs(exact).max()          # leave out entries that are numerically zero
    err = []
    for h in (0.4, 0.2):
        fd = (ac_dv(st, h * d, delta) + ac_dv(st, -h * d, delta)) / (2.0 * h * h)
        err.append(np.abs(fd - exact)[big].max() / np.abs(exact)[big].max())
    # measured: 4e-6 to 3e-5 at h = 0.4 and about 1e-6 at h = 0.2 (the gap is the O(h^4) term of the solution;
    # below that the accuracy of the AC solver takes over, so the h^2 rate is not asserted)
    assert err[0] < 2e-4 and err[1] < 2e-5 and err[1] < err[0]


def test_error_after_the_quadratic_term_is_third_order(model):
    """Halving delta must cut |dV_AC - c'eps - eps'H eps| by about eight."""
    st, _, delta = model
    eps = sample_eps(np.random.default_rng(2), st.src.n, 100)
    errs = []
    for dl in (2 * delta, delta, delta / 2):
        so = SecondOrder(st, dl)
        errs.append(np.abs(ac_dv(st, eps, dl) - so.c_lin @ eps - so.quadratic(eps)).max())
    assert errs[0] / errs[1] == pytest.approx(8.0, rel=0.3)
    assert errs[1] / errs[2] == pytest.approx(8.0, rel=0.3)


@pytest.mark.parametrize("case", ["case33bw", "case118"])
def test_branch_admittances_and_second_order_state(case):
    import scipy.sparse as sp
    g = _grid(case)
    nb, n = g.Yf.shape[0], g.n_bus
    Cf = sp.csr_matrix((np.ones(nb), (np.arange(nb), g.br_f)), shape=(nb, n))
    Ct = sp.csr_matrix((np.ones(nb), (np.arange(nb), g.br_t)), shape=(nb, n))
    D = (Cf.T @ g.Yf + Ct.T @ g.Yt - g.Ybus).tocoo()                 # only bus shunts may remain
    assert np.all(np.abs(D.data[D.row != D.col]) < 1e-9)
    st = make_study(g)
    delta = 0.1 if case == "case33bw" else 0.05
    so = SecondOrder(st, delta)
    eps = sample_eps(np.random.default_rng(3), st.src.n, 50)
    dv = so.state(eps)
    ac = ac_dv(st, eps, delta)
    V, _, _ = st.solver.solve(st.injections(eps, delta))
    i_lin = np.abs(g.Yf @ (st.V_star[:, None] + so.E @ eps))
    i_so = np.abs(g.Yf @ (st.V_star[:, None] + dv))
    i_ac = np.abs(g.Yf @ V)
    assert np.abs(i_so - i_ac).max() < 0.05 * np.abs(i_lin - i_ac).max()  # the second-order state removes most of the error
    vm = np.abs(st.V_star[g.pq, None] + dv[g.pq]) - st.vm_star[:, None]
    assert np.abs(vm - ac).max() < 0.2 * np.abs(so.c_lin @ eps - ac).max()
