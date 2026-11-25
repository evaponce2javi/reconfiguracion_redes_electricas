"""Closed-form checks of the uncertainty-set geometry."""
import numpy as np
import pytest

from cdhim import geometry as geo
from cdhim.core import MethodConfig, branch_covariances, run_cdhim
from cdhim.uncertainty import SourceCovariance

CHI2 = geo.chi2_quantile(0.9973, 2)


def random_pairs(rng, n_branch=200, n_u=40):
    """Coefficient rows for ``n_branch`` independent bus pairs with a wide range of correlations."""
    C = rng.standard_normal((2 * n_branch, n_u))
    mix = rng.uniform(0.0, 1.0, n_branch)[:, None]
    C[1::2] = mix * C[0::2] + (1.0 - mix) * C[1::2] * rng.uniform(0.01, 3.0, (n_branch, 1))
    fi, ti = np.arange(0, 2 * n_branch, 2), np.arange(1, 2 * n_branch, 2)
    return C, fi, ti


def test_chi2_quantile_and_bonferroni():
    assert CHI2 == pytest.approx(11.829, abs=1e-3)
    assert CHI2 == pytest.approx(-2.0 * np.log(1 - 0.9973), rel=1e-12)           # chi2_2 is exponential
    assert geo.chi2_quantile_bonferroni(0.9973, 1) == pytest.approx(CHI2, rel=1e-12)
    m = 2545
    assert geo.chi2_quantile_bonferroni(0.9973, m) == pytest.approx(geo.chi2_quantile(1 - 0.0027 / m, 2), rel=1e-9)


def test_eig2x2_matches_lapack():
    rng = np.random.default_rng(0)
    C, fi, ti = random_pairs(rng)
    q_ii, q12 = branch_covariances(C, fi, ti)
    lam_min, lam_max = geo.eig2x2(q_ii[fi], q_ii[ti], q12)
    for k in range(len(fi)):
        ref = np.linalg.eigvalsh([[q_ii[fi[k]], q12[k]], [q12[k], q_ii[ti[k]]]])
        assert lam_min[k] == pytest.approx(ref[0], rel=1e-9, abs=1e-14)
        assert lam_max[k] == pytest.approx(ref[1], rel=1e-12)


@pytest.mark.parametrize("gain_box,kappa", [("covariance_1sigma", np.pi / 4), ("covariance_3sigma", np.pi * CHI2 / 36)])
def test_gain_equals_closed_form(gain_box, kappa):
    """G = 1 - kappa * sqrt(1 - rho^2) with a branch-independent kappa for both covariance boxes."""
    rng = np.random.default_rng(1)
    C, fi, ti = random_pairs(rng)
    res = run_cdhim(C, fi, ti, MethodConfig(gain_box=gain_box))
    assert geo.kappa_constant(gain_box, CHI2) == pytest.approx(kappa)
    np.testing.assert_allclose(res.gain, geo.gain_closed_form(res.rho, kappa), atol=1e-10)
    stable = np.abs(res.rho) < 0.999
    np.testing.assert_allclose(res.kappa[stable], kappa, rtol=1e-7)


def test_gain_with_affine_box_is_not_a_function_of_rho():
    rng = np.random.default_rng(2)
    C, fi, ti = random_pairs(rng)
    res = run_cdhim(C, fi, ti, MethodConfig(gain_box="affine"))
    r = np.abs(C).sum(axis=1)
    lam_min, lam_max = geo.eig2x2(res.q11, res.q22, res.q12)
    expected = np.maximum(0.0, 1.0 - np.pi * CHI2 * np.sqrt(lam_min * lam_max) / (4.0 * r[fi] * r[ti]))
    np.testing.assert_allclose(res.gain, expected, atol=1e-12)
    assert geo.kappa_constant("affine", CHI2) is None
    assert np.nanstd(res.kappa) > 1e-3                      # kappa changes from branch to branch


def test_paper_threshold_is_inactive_with_the_original_box():
    """With kappa = pi/4, G >= 1 - pi/4 = 21.46 % for every rho, so gamma_th = 15 % never binds."""
    assert geo.gain_closed_form(0.0, np.pi / 4) == pytest.approx(0.2146, abs=1e-4)
    assert geo.rho_threshold(0.15, np.pi / 4) == 0.0
    assert geo.rho_threshold(0.25, np.pi / 4) == pytest.approx(0.2968, abs=1e-4)
    assert geo.rho_threshold(0.15, np.pi * CHI2 / 36) == pytest.approx(0.5674, abs=1e-3)


@pytest.mark.parametrize("gain_box", ["covariance_1sigma", "covariance_3sigma"])
@pytest.mark.parametrize("gamma", [0.05, 0.15, 0.25, 0.40])
def test_gain_rule_equals_rho_rule_for_constant_kappa(gain_box, gamma):
    rng = np.random.default_rng(3)
    C, fi, ti = random_pairs(rng, n_branch=500)
    kappa = geo.kappa_constant(gain_box, CHI2)
    by_gain = run_cdhim(C, fi, ti, MethodConfig(gain_box=gain_box, gamma_th=gamma))
    by_rho = run_cdhim(C, fi, ti, MethodConfig(gain_box=gain_box, selection_rule="rho", rho_th=geo.rho_threshold(gamma, kappa)))
    np.testing.assert_array_equal(by_gain.ellipse, by_rho.ellipse)


def test_ellipse_projection_against_analytic_maximum():
    """max |x_i| over the ellipse x'Q^-1x = chi2 is sqrt(chi2 Q_ii); its largest norm is sqrt(chi2 lam_max)."""
    rng = np.random.default_rng(4)
    t = np.linspace(0.0, 2.0 * np.pi, 200001)
    circle = np.vstack([np.cos(t), np.sin(t)])
    for _ in range(20):
        A = rng.standard_normal((2, 2)) * rng.uniform(0.1, 5.0)
        Q = A @ A.T
        boundary = np.linalg.cholesky(Q) @ circle * np.sqrt(CHI2)
        lam_min, lam_max = geo.eig2x2(Q[0, 0], Q[1, 1], Q[0, 1])
        assert np.abs(boundary[0]).max() == pytest.approx(geo.radius_coordinate(Q[0, 0], CHI2), rel=1e-8)
        assert np.abs(boundary[1]).max() == pytest.approx(geo.radius_coordinate(Q[1, 1], CHI2), rel=1e-8)
        assert np.linalg.norm(boundary, axis=0).max() == pytest.approx(geo.radius_lambda_max(lam_max, CHI2), rel=1e-8)
        assert geo.radius_lambda_max(lam_max, CHI2) >= geo.radius_coordinate(max(Q[0, 0], Q[1, 1]), CHI2) - 1e-12
        d2 = geo.mahalanobis2(boundary[:1], boundary[1:], np.array([Q[0, 0]]), np.array([Q[1, 1]]), np.array([Q[0, 1]]))
        np.testing.assert_allclose(d2, CHI2, rtol=1e-8)


def test_default_source_covariance_is_one_third_identity():
    """Q = C Sigma_eps C^T with Sigma_eps = I/3 must equal (1/3) C C^T, however Sigma_eps is supplied."""
    rng = np.random.default_rng(5)
    C, fi, ti = random_pairs(rng, n_branch=50, n_u=30)
    ref = C @ C.T / 3.0
    for sigma in (None, SourceCovariance(), SourceCovariance(var=np.full(30, 1 / 3)), SourceCovariance.dense(np.eye(30) / 3)):
        q_ii, q_ij = branch_covariances(C, fi, ti, sigma)
        np.testing.assert_allclose(q_ii, np.diag(ref), rtol=1e-12)
        np.testing.assert_allclose(q_ij, ref[fi, ti], rtol=1e-11, atol=1e-13)


def test_equicorrelated_source_covariance_matches_dense():
    rng = np.random.default_rng(6)
    C, fi, ti = random_pairs(rng, n_branch=30, n_u=20)
    var = rng.uniform(0.05, 0.5, 20)
    s = np.sqrt(var)
    dense = 0.4 * np.outer(s, s) + 0.6 * np.diag(var)
    a = branch_covariances(C, fi, ti, SourceCovariance(var=var, rho_c=0.4))
    b = branch_covariances(C, fi, ti, SourceCovariance.dense(dense))
    np.testing.assert_allclose(a[0], b[0], rtol=1e-11)
    np.testing.assert_allclose(a[1], b[1], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(SourceCovariance(var=var, rho_c=0.4).matrix(20), dense, rtol=1e-12)


def test_mahalanobis_is_invariant_to_the_uncertainty_level():
    """delta scales every coefficient alike, so d2 of linear samples cannot depend on delta (Table I)."""
    rng = np.random.default_rng(7)
    C, fi, ti = random_pairs(rng, n_branch=20, n_u=25)
    eps = rng.uniform(-1, 1, (25, 500))
    ref = None
    for delta in (0.05, 0.10, 0.20, 0.30):
        Cd = C * delta
        q_ii, q12 = branch_covariances(Cd, fi, ti)
        x = Cd @ eps
        d2 = geo.mahalanobis2(x[fi], x[ti], q_ii[fi], q_ii[ti], q12)
        ref = d2 if ref is None else ref
        np.testing.assert_allclose(d2, ref, rtol=1e-9)


def test_output_radius_modes_and_combination():
    rng = np.random.default_rng(8)
    C = rng.standard_normal((3, 12))
    C[1] = 0.8 * C[0] + 0.2 * C[1]
    fi, ti = np.array([0, 1]), np.array([1, 2])
    res = run_cdhim(C, fi, ti, MethodConfig(delta_stab=1e-6))
    assert res.ellipse.all() and res.bus_ellipse.all()
    assert res.r_out[1] == pytest.approx(res.r_lmax.min())            # shared bus keeps the smaller radius
    res_max = run_cdhim(C, fi, ti, MethodConfig(delta_stab=1e-6, bus_combine="max"))
    assert res_max.r_out[1] == pytest.approx(res.r_lmax.max())
    coord = run_cdhim(C, fi, ti, MethodConfig(delta_stab=1e-6, radius="coordinate"))
    np.testing.assert_allclose(coord.r_out, np.sqrt(CHI2 * res.q_ii))
    assert np.all(coord.r_out <= res.r_out + 1e-15)                   # lambda_max never beats the projection
    np.testing.assert_allclose(res.overconservatism().min() >= 1.0, True)
    clipped = run_cdhim(C, fi, ti, MethodConfig(delta_stab=1e-6, clip_to_affine_box=True))
    assert np.all(clipped.r_out <= res.r_affine + 1e-15)


def test_stability_guard_sends_near_singular_pairs_to_box():
    C = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0001], [3.0, -1.0, 0.5]])
    res = run_cdhim(C, np.array([0, 0]), np.array([1, 2]), MethodConfig())
    assert abs(res.rho[0]) > 0.999 and not res.ellipse[0]
    assert res.ellipse[1]
    assert res.r_out[1] == pytest.approx(np.abs(C[1]).sum())          # bus 1 only touches a Box branch
