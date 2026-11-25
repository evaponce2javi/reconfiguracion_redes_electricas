"""Support, mean and variance of every marginal of experiment C3, and the copula sampler."""
import numpy as np
import pytest

from cdhim.uncertainty import (UNIFORM_VARIANCE, copula_to_pearson, copula_to_pearson_uniform,
                               estimate_source_covariance, marginal_catalogue, sample_eps)

CAT = marginal_catalogue()
N = 400_000

# name -> (support lo, support hi, variance); None = not checked analytically here
EXPECTED = {
    "uniform": (-1.0, 1.0, 1 / 3),
    "beta25_support": (-0.4, 1.0, (10 / 392) / (5 / 7) ** 2),            # Var(Beta(2,5)) = 10/392
    "beta25_var13": (None, None, 1 / 3),
    "beta25_paper": (-1.14, 2.86, 16 * 10 / 392),
    "beta52_support": (-1.0, 0.4, (10 / 392) / (5 / 7) ** 2),
    "beta0505_ushape": (-1.0, 1.0, 0.5),
    "triangular": (-1.0, 1.0, 1 / 6),
    "student_t3": (-np.inf, np.inf, 1 / 3),
    "lognormal": (None, np.inf, 1 / 3),
}


@pytest.mark.parametrize("name", sorted(CAT))
def test_sampled_moments_match_declared(name):
    m = CAT[name]
    x = m.sample(np.random.default_rng(11), N)
    lo, hi = m.support
    assert x.min() >= lo - 1e-9 and x.max() <= hi + 1e-9
    heavy = name in ("student_t3", "lognormal")
    assert x.mean() == pytest.approx(m.mean, abs=0.02 if heavy else 4e-3)
    assert x.var() == pytest.approx(m.variance, rel=0.25 if heavy else 0.02)
    assert np.mean(np.abs(x) > 1.0) == pytest.approx(m.fraction_outside_unit, abs=2e-3)
    assert m.in_model == (lo >= -1 - 1e-9 and hi <= 1 + 1e-9)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_declared_support_and_variance(name):
    m = CAT[name]
    lo, hi, var = EXPECTED[name]
    if lo is not None:
        assert m.support[0] == pytest.approx(lo, abs=1e-9)
    if hi is not None:
        assert m.support[1] == pytest.approx(hi, abs=1e-9)
    assert m.variance == pytest.approx(var, rel=1e-6)


def test_zero_mean_except_the_paper_transform():
    for name, m in CAT.items():
        if name == "beta25_paper":
            assert m.mean == pytest.approx((2 / 7 - 0.285) * 4.0, abs=1e-12)     # +0.0029 bias in the original code
        else:
            assert m.mean == pytest.approx(0.0, abs=1e-7), name


def test_in_model_flags():
    inside = {n for n, m in CAT.items() if m.in_model}
    assert inside == {"uniform", "beta25_support", "beta52_support", "beta0505_ushape", "triangular",
                      "truncnorm_3s", "student_t3_trunc", "lognormal_trunc"}
    assert CAT["beta25_paper"].fraction_outside_unit > 0.05                      # the manuscript's Beta test left [-1, 1]
    assert CAT["beta25_var13"].fraction_outside_unit > 0.02


def test_truncated_normal_variance():
    from scipy import stats

    s = 1 / 3
    exact = s**2 * (1 - 2 * 3 * stats.norm.pdf(3) / (2 * stats.norm.cdf(3) - 1))
    assert CAT["truncnorm_3s"].variance == pytest.approx(exact, rel=1e-6)
    assert CAT["truncnorm_3s"].support == pytest.approx((-1.0, 1.0))


def test_uniform_sampler_and_batching_invariance():
    a = sample_eps(np.random.default_rng(5), 7, 1000)
    rng = np.random.default_rng(5)
    b = np.hstack([sample_eps(rng, 7, 300), sample_eps(rng, 7, 700)])
    np.testing.assert_array_equal(a, b)                    # same stream however the samples are batched
    assert a.shape == (7, 1000) and a.min() >= -1 and a.max() <= 1
    big = sample_eps(np.random.default_rng(6), 4, N)
    assert big.var(axis=1) == pytest.approx(np.full(4, UNIFORM_VARIANCE), rel=0.01)
    assert np.abs(np.corrcoef(big) - np.eye(4)).max() < 0.01


@pytest.mark.parametrize("r", [0.3, 0.6, 0.9])
def test_gaussian_copula_keeps_marginals_and_sets_correlation(r):
    eps = sample_eps(np.random.default_rng(7), 5, N, copula_r=r)
    assert eps.min() >= -1 and eps.max() <= 1
    assert eps.var(axis=1) == pytest.approx(np.full(5, UNIFORM_VARIANCE), rel=0.01)
    off = np.corrcoef(eps)[~np.eye(5, dtype=bool)]
    assert off == pytest.approx(np.full(20, copula_to_pearson_uniform(r)), abs=0.01)
    assert copula_to_pearson(CAT["uniform"], r) == pytest.approx(copula_to_pearson_uniform(r), abs=2e-4)


def test_copula_with_skewed_marginal():
    m = CAT["beta25_support"]
    eps = sample_eps(np.random.default_rng(8), 3, 200_000, marginal=m, copula_r=0.6)
    assert eps.min() >= -0.4 - 1e-9 and eps.max() <= 1 + 1e-9
    assert eps.var(axis=1) == pytest.approx(np.full(3, m.variance), rel=0.02)
    off = np.corrcoef(eps)[~np.eye(3, dtype=bool)]
    assert off == pytest.approx(np.full(6, copula_to_pearson(m, 0.6)), abs=0.01)


def test_estimate_source_covariance_from_series():
    rng = np.random.default_rng(9)
    L = np.array([[1.0, 0.0, 0.0], [0.7, 0.5, 0.0], [0.2, 0.1, 0.9]])
    x = 10.0 + (rng.uniform(-1, 1, (50_000, 3)) @ L.T) * np.array([2.0, 0.5, 4.0])
    x0, h, mean_eps, sigma = estimate_source_covariance(x)
    eps = (x - x0) / h
    assert eps.min() == pytest.approx(-1.0) and eps.max() == pytest.approx(1.0)
    np.testing.assert_allclose(sigma.matrix(3), np.cov(eps, rowvar=False), rtol=1e-12)
    assert np.all(np.linalg.eigvalsh(sigma.matrix(3)) > 0)
    corr = sigma.matrix(3) / np.sqrt(np.outer(np.diag(sigma.matrix(3)), np.diag(sigma.matrix(3))))
    np.testing.assert_allclose(corr, np.corrcoef(x, rowvar=False), atol=1e-10)
