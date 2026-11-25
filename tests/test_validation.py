"""Monte Carlo validation loop on IEEE-33: invariants that must hold whatever the random draw."""
import warnings

import numpy as np
import pytest

from cdhim import SourceCovariance, load_grid, make_study, marginal_catalogue, monte_carlo, standard_bounds
from cdhim.core import MethodConfig
from cdhim.uncertainty import UNIFORM_VARIANCE, copula_to_pearson_uniform


@pytest.fixture(scope="module")
def st():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_study(load_grid("case33bw"))


def test_affine_box_encloses_every_linear_sample(st):
    res = st.run(0.15)
    cs = monte_carlo(st, 0.15, res, standard_bounds(res), 3000, seed=1, ac=False)
    assert cs.n == 3000 and cs.n_ac == 0
    assert cs.bus_lin["affine_box"].sum() == 0                       # deterministic enclosure of the linear model
    assert cs.any_lin["affine_box"] == 0
    assert cs.bus_lin["cltm_3sigma"].sum() > 0                       # the 3-sigma box is only probabilistic
    assert cs.d2_lin_max.max() > 0.5


def test_ac_and_linear_statistics_are_consistent(st):
    res = st.run(0.15)
    res_c = st.run(0.15, MethodConfig(radius="coordinate"))
    cs = monte_carlo(st, 0.15, res, standard_bounds(res, res_c), 4000, seed=2)
    assert cs.n_ac == 4000 and cs.solver.converged == 4000
    s = cs.bus_summary("cdhim_lmax", "ac")
    c = cs.bus_summary("cdhim_coord", "ac")
    assert s["max"] <= c["max"] + 1e-12                              # lambda_max intervals contain the projections
    assert s["system"] <= c["system"] + 1e-12
    assert s["system"] <= s["expected_per_sample"] + 1e-12           # P(any) <= E[count]
    b = cs.branch_summary("lin")
    assert b["n"] == res.n_ellipse and b["mean"] < 0.0027            # sums of uniforms: lighter tails than Gaussian
    assert cs.err_abs.max < 1e-4                                     # second-order error at delta = 15 %
    assert 0 < cs.err_abs.quantile(0.99) <= cs.err_abs.max


def test_same_seed_gives_same_statistics_whatever_the_batch(st):
    res = st.run(0.15)
    a = monte_carlo(st, 0.15, res, standard_bounds(res), 1500, seed=3, batch=1500, ac=False)
    b = monte_carlo(st, 0.15, res, standard_bounds(res), 1500, seed=3, batch=400, ac=False)
    np.testing.assert_array_equal(a.br_lin, b.br_lin)
    np.testing.assert_array_equal(a.bus_lin["cltm_3sigma"], b.bus_lin["cltm_3sigma"])


def test_scaled_ellipses_are_nested(st):
    res = st.run(0.15)
    cs = monte_carlo(st, 0.15, res, {"cdhim_lmax": res.r_out}, 3000, seed=4, ac=False,
                     marginal=marginal_catalogue()["beta0505_ushape"], ellipse_scales={"half": 0.5, "double": 2.0})
    base = cs.br_lin
    assert np.all(cs.br_scaled[("half", "lin")][0] >= base)
    assert np.all(cs.br_scaled[("double", "lin")][0] <= base)
    assert cs.branch_summary("lin", "half")["mean"] > cs.branch_summary("lin")["mean"]


def test_correlated_sources_need_the_full_covariance(st):
    """With strongly correlated symbols, Q = (1/3) C C^T under-covers and Q = C Sigma_eps C^T does not."""
    r = 0.9
    sigma = SourceCovariance(var=UNIFORM_VARIANCE, rho_c=float(copula_to_pearson_uniform(r)))
    wrong, right = st.run(0.15), st.run(0.15, sigma_eps=sigma)
    assert np.all(right.q_ii > wrong.q_ii)                           # loads moving together widen every bus
    cw = monte_carlo(st, 0.15, wrong, {"cltm_3sigma": wrong.r_cltm}, 4000, seed=5, copula_r=r, ac=False)
    cr = monte_carlo(st, 0.15, right, {"cltm_3sigma": right.r_cltm}, 4000, seed=5, copula_r=r, ac=False)
    assert cw.bus_summary("cltm_3sigma", "lin")["mean"] > 0.05       # badly under-covered
    assert cr.bus_summary("cltm_3sigma", "lin")["mean"] < 0.01
