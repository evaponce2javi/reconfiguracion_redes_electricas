"""Targets built from the affine forms (cdhim.affine_forms)."""
import warnings

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid, make_study, sample_eps, spectral_zones
from cdhim.affine_forms import (TargetModel, affine_radius, bernstein_radius, bus_targets, drop_targets, sigma)


@pytest.fixture(scope="module")
def few_symbols():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = load_grid("case33bw")
    st = make_study(g, zones=spectral_zones(g, 8))
    so = SecondOrder(st, 0.3)
    fi, ti, _ = g.pq_pq_branches()
    return so, so.forms(), fi, ti


@pytest.mark.parametrize("kind", ["bus", "drop"])
def test_matrix_free_products_match_the_assembled_forms(few_symbols, kind):
    so, H, fi, ti = few_symbols
    n = H.shape[0]
    W = bus_targets(n) if kind == "bus" else drop_targets(n, fi, ti)
    tm = TargetModel(so, W)
    Ht = np.einsum("it,ikl->tkl", W.toarray(), H)
    idx = np.arange(W.shape[1])
    X = np.sign(np.random.default_rng(0).normal(size=(H.shape[1], idx.size)))
    scale = np.abs(Ht).max()
    np.testing.assert_allclose(tm.hx(X, idx), np.einsum("tkl,lt->kt", Ht, X), atol=1e-12 * scale)
    np.testing.assert_allclose(tm.hdiag(), np.einsum("tkk->tk", Ht), atol=1e-12 * scale)
    exact = np.einsum("tk,kt->t", tm.c, X) + np.einsum("kt,tkl,lt->t", X, Ht, X)
    np.testing.assert_allclose(tm.value(X, idx), exact, atol=1e-12 * scale)


@pytest.mark.parametrize("side", [1, -1])
def test_search_reaches_the_vertex_optimum_on_a_small_case(few_symbols, side):
    so, H, fi, ti = few_symbols
    tm = TargetModel(so, drop_targets(H.shape[0], fi, ti))
    K = H.shape[1]
    V = np.array(np.meshgrid(*[[-1.0, 1.0]] * K)).reshape(K, -1)
    Ht = np.einsum("it,ikl->tkl", tm.W.toarray(), H)
    full = tm.c @ V + np.einsum("kv,tkl,lv->tv", V, Ht, V)
    best = full.max(axis=1) if side > 0 else full.min(axis=1)
    X, f, f_lin, _ = tm.worst_vertices(side)
    assert np.all(side * (f - f_lin) >= -1e-15)                  # never worse than the first-order vertex
    np.testing.assert_allclose(f, best, atol=1e-12)


def test_bernstein_radius_covers_and_is_capped():
    rng = np.random.default_rng(1)
    C = np.vstack([rng.normal(size=(5, 400)), np.r_[1.0, np.zeros(399)][None, :]])
    alpha = 0.0027
    r = bernstein_radius(C, alpha)
    assert np.all(r <= affine_radius(C) + 1e-15)
    assert r[-1] == pytest.approx(1.0)                           # a single symbol: the cap is the support
    s = C @ sample_eps(rng, 400, 200000)
    assert np.all((np.abs(s) > r[:, None]).mean(axis=1) <= alpha)
    assert np.all(r[:5] > 3 * sigma(C)[:5])                       # distribution-free costs width over 3 sigma


def test_factor_model_correlation_support_and_rewrite():
    from cdhim.affine_forms import FactorModel
    zone = np.r_[0, 0, 0, 1, 1]
    for rho in (0.155, 0.5, 0.845):
        fm = FactorModel(zone, FactorModel.w_for_rho(rho))
        assert fm.rho == pytest.approx(rho)
        e = fm.sample(np.random.default_rng(0), 200000)
        assert np.abs(e).max() <= 1.0
        R = np.corrcoef(e)
        assert R[0, 1] == pytest.approx(rho, abs=0.01) and R[0, 3] == pytest.approx(0.0, abs=0.01)
    # the affine forms on the independent symbols give the same values
    fm = FactorModel(zone, 0.6)
    rng = np.random.default_rng(1)
    C = rng.normal(size=(4, 5))
    z, h = rng.uniform(-1, 1, (2, 7)), rng.uniform(-1, 1, (5, 7))
    np.testing.assert_allclose(C @ (fm.w * z[zone] + (1 - fm.w) * h), fm.coefficients(C) @ np.vstack([z, h]), atol=1e-14)
    assert np.all(affine_radius(fm.coefficients(C)) <= affine_radius(C) + 1e-14)


def test_order_statistic_bound_of_the_sampled_range():
    from cdhim.affine_forms import certified_split
    rng = np.random.default_rng(5)
    n, trials = 50, 200000
    q = rng.exponential(size=(trials, n + 1))                    # any continuous law
    p = np.mean(q[:, -1] > q[:, :-1].max(axis=1))
    assert p == pytest.approx(1.0 / (n + 1), abs=4 * np.sqrt(p / trials))
    assert certified_split(0.0027, 4000) == pytest.approx(0.0027 - 2 / 4001)
    with pytest.raises(ValueError):
        certified_split(0.0027, 500)
