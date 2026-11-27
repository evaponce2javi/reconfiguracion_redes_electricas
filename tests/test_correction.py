"""Second-order correction of the CD-HIM bounds (cdhim.correction)."""
import warnings

import numpy as np
import pytest

from cdhim import load_grid, make_study, method_config, load_config, sample_eps
from cdhim import geometry as geo
from cdhim.correction import principal_axes, second_order_correction, vertex_values


@pytest.fixture(scope="module")
def setup():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = load_grid("case33bw")
    st = make_study(g)
    delta = 0.15
    res = st.run(delta, method_config(load_config()))
    cor = second_order_correction(st, delta, res, n_samples=4000, seed=3)
    return st, delta, res, cor


def test_exact_mean_matches_the_sampled_mean(setup):
    st, delta, res, cor = setup
    q = cor.so.quadratic(sample_eps(np.random.default_rng(11), st.src.n, 20000))
    se = q.std(axis=1) / np.sqrt(q.shape[1])
    assert np.all(np.abs(q.mean(axis=1) - cor.mu) < 5 * se + 1e-15)


def test_sampled_variance_matches_the_exact_one(setup):
    # Var(eps' H eps) for iid U[-1, 1]: 2 s^4 |H|_F^2 + (mu4 - 3 s^4) sum H_kk^2, s^2 = 1/3, mu4 = 1/5
    st, delta, res, cor = setup
    H = cor.so.forms()
    exact = 2.0 / 9.0 * np.einsum("ikl,ikl->i", H, H) - 2.0 / 15.0 * np.einsum("ikk,ikk->i", H, H)
    assert np.all(np.abs(cor.var_q - exact) <= 0.2 * exact + 1e-30)


def test_principal_axes_diagonalise_the_covariance():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(50, 2, 2))
    Q = a @ a.transpose(0, 2, 1)
    ux, uy, vx, vy = principal_axes(Q[:, 0, 0], Q[:, 1, 1], Q[:, 0, 1])
    lmin, lmax = geo.eig2x2(Q[:, 0, 0], Q[:, 1, 1], Q[:, 0, 1])
    u, v = np.stack([ux, uy], 1), np.stack([vx, vy], 1)
    np.testing.assert_allclose(np.einsum("bi,bij,bj->b", u, Q, u), lmax, rtol=1e-10)
    np.testing.assert_allclose(np.einsum("bi,bij,bj->b", v, Q, v), lmin, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(np.einsum("bi,bi->b", u, v), 0.0, atol=1e-12)


def test_minkowski_ellipse_contains_ellipse_plus_range_box(setup):
    st, delta, res, cor = setup
    m = res.ellipse
    cf, ct, q11, q22, q12 = cor.ellipse(res, "minkowski")
    ux, uy, vx, vy = principal_axes(res.q11[m], res.q22[m], res.q12[m])
    th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    a, b = np.sqrt(res.chi2 * res.lam_max[m]), np.sqrt(res.chi2 * np.maximum(res.lam_min[m], 0))
    for su in (cor.pu_min[m], cor.pu_max[m]):            # every corner of the range box
        for sv in (cor.pv_min[m], cor.pv_max[m]):
            pu = a[:, None] * np.cos(th) + su[:, None]
            pv = b[:, None] * np.sin(th) + sv[:, None]
            x = pu * ux[:, None] + pv * vx[:, None]
            y = pu * uy[:, None] + pv * vy[:, None]
            d2 = geo.mahalanobis2(x - cf[:, None], y - ct[:, None], q11, q22, q12)
            assert np.all(d2 <= res.chi2 * (1 + 1e-9))


def test_vertex_values_are_the_quadratic_at_the_extreme_vertices(setup):
    st, delta, res, cor = setup
    i = int(np.argmax(np.abs(cor.so.c_lin).sum(axis=1)))
    s = np.sign(cor.so.c_lin[i])[:, None]
    assert cor.q_vhi[i] == pytest.approx(cor.so.quadratic(s)[i, 0], rel=1e-12)
    assert cor.q_vlo[i] == pytest.approx(cor.so.quadratic(-s)[i, 0], rel=1e-12)
    mask = np.zeros(res.r_affine.size, bool)
    mask[i] = True
    vlo, vhi = vertex_values(cor.so, mask)
    assert vhi[i] == pytest.approx(cor.q_vhi[i], rel=1e-12) and np.count_nonzero(vhi) == 1


def test_no_correction_returns_the_linear_bounds(setup):
    st, delta, res, cor = setup
    lo, hi = cor.interval(res.r_out, res, "none")
    np.testing.assert_array_equal(hi, res.r_out)
    np.testing.assert_array_equal(lo, -res.r_out)
    lo, hi = cor.cdhim(res, "range")
    assert np.all(lo <= -res.r_out + 1e-15) and np.all(hi >= res.r_out - 1e-15)
