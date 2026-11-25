"""Uncertainty-set geometry for a pair of terminal voltages.

Every function is vectorised over branches and has no side effects. The three "boxes" that the
manuscript used interchangeably are kept apart, each under its own name:

* covariance box  : half-widths  k * sqrt(Q_ii)            -> :func:`cov_box_halfwidth`
* CLTM 3-sigma box: the covariance box with k = 3          -> :func:`cltm_box_halfwidth`
* affine box      : half-widths  sum_k |c_ik| (L1 norm)    -> :func:`affine_box_halfwidth`
"""
from __future__ import annotations

import numpy as np
from scipy import stats

GAIN_BOXES = ("covariance_1sigma", "covariance_3sigma", "affine")


# --------------------------------------------------------------------------- quantiles
def chi2_quantile(confidence: float, dof: int = 2) -> float:
    """Chi-square quantile; 11.829 for confidence 0.9973 and 2 degrees of freedom."""
    return float(stats.chi2.ppf(confidence, dof))


def bonferroni_confidence(confidence: float, m: int) -> float:
    """Per-set confidence 1 - alpha/m that gives simultaneous confidence >= 1 - alpha over m sets."""
    return 1.0 - (1.0 - confidence) / max(int(m), 1)


def chi2_quantile_bonferroni(confidence: float, m: int) -> float:
    """chi2_2 quantile at level 1 - alpha/m, computed as -2 ln(alpha/m) to keep precision."""
    return float(-2.0 * np.log((1.0 - confidence) / max(int(m), 1)))


# --------------------------------------------------------------------------- 2x2 spectra
def eig2x2(q11, q22, q12):
    """Eigenvalues (lam_min, lam_max) of symmetric 2x2 matrices [[q11, q12], [q12, q22]]."""
    q11, q22, q12 = np.asarray(q11, float), np.asarray(q22, float), np.asarray(q12, float)
    half_tr = 0.5 * (q11 + q22)
    rad = np.hypot(0.5 * (q11 - q22), q12)
    lam_max = half_tr + rad
    # lam_min from the determinant: avoids the cancellation of half_tr - rad when |rho| -> 1
    det = np.maximum(q11 * q22 - q12 * q12, 0.0)
    lam_min = np.divide(det, lam_max, out=np.zeros_like(lam_max), where=lam_max > 0)
    return lam_min, lam_max


def correlation(q11, q22, q12):
    den = np.sqrt(np.asarray(q11, float) * np.asarray(q22, float))
    return np.divide(q12, den, out=np.zeros_like(den), where=den > 0)


# --------------------------------------------------------------------------- boxes
def cov_box_halfwidth(qii, k: float = 1.0):
    """Covariance box: k standard deviations per coordinate."""
    return k * np.sqrt(np.asarray(qii, float))


def cltm_box_halfwidth(qii):
    """3-sigma probabilistic box ("CLTM Box (3 sigma)" of the manuscript, Fig. 4)."""
    return cov_box_halfwidth(qii, 3.0)


def affine_box_halfwidth(C):
    """Deterministic affine-arithmetic enclosure: L1 norm of each row of the coefficient matrix."""
    return np.abs(C).sum(axis=1)


def box_area(h_i, h_j):
    return 4.0 * np.asarray(h_i, float) * np.asarray(h_j, float)


# --------------------------------------------------------------------------- ellipse
def ellipse_area(lam_min, lam_max, scale2: float = 1.0):
    """Area of {x : x' Q^-1 x <= scale2}: pi * scale2 * sqrt(lam_min * lam_max)."""
    return np.pi * scale2 * np.sqrt(np.asarray(lam_min, float) * np.asarray(lam_max, float))


def radius_lambda_max(lam_max, chi2: float):
    """Manuscript Eq. (11): radius of the circle circumscribing the confidence ellipse."""
    return np.sqrt(chi2 * np.asarray(lam_max, float))


def radius_coordinate(qii, chi2: float):
    """Exact projection of the confidence ellipse on one coordinate axis."""
    return np.sqrt(chi2 * np.asarray(qii, float))


def mahalanobis2(x, y, q11, q22, q12):
    """Squared Mahalanobis distance of points (x, y) for per-branch 2x2 covariances.

    ``x`` and ``y`` have shape (n_branch, n_samples); the covariance entries have shape (n_branch,).
    Computed in standardised coordinates, d2 = (u^2 - 2 rho u v + v^2) / (1 - rho^2). No diagonal
    regularisation is added: an absolute jitter is not scale invariant and inflates the ellipse.
    """
    sx, sy = np.sqrt(q11)[:, None], np.sqrt(q22)[:, None]
    rho = correlation(q11, q22, q12)[:, None]
    u, v = x / sx, y / sy
    return (u * u - 2.0 * rho * u * v + v * v) / (1.0 - rho * rho)


# --------------------------------------------------------------------------- coupling gain
def coupling_gain(area_ell, area_box):
    """G = max(0, 1 - A_ell / A_box), as a fraction in [0, 1]."""
    ratio = np.divide(area_ell, area_box, out=np.ones_like(np.asarray(area_box, float)), where=np.asarray(area_box) > 0)
    return np.maximum(0.0, 1.0 - ratio)


def kappa_constant(gain_box: str, chi2: float) -> float | None:
    """kappa such that G = 1 - kappa * sqrt(1 - rho^2); None when it is branch dependent."""
    if gain_box == "covariance_1sigma":
        return np.pi / 4.0
    if gain_box == "covariance_3sigma":
        return np.pi * chi2 / 36.0
    if gain_box == "affine":
        return None
    raise ValueError(f"unknown gain_box {gain_box!r}; expected one of {GAIN_BOXES}")


def gain_closed_form(rho, kappa: float):
    return np.maximum(0.0, 1.0 - kappa * np.sqrt(np.maximum(1.0 - np.asarray(rho, float) ** 2, 0.0)))


def rho_threshold(gamma_th: float, kappa: float) -> float:
    """|rho| above which G >= gamma_th. Returns 0.0 when every branch satisfies the gain test."""
    x = (1.0 - gamma_th) / kappa
    return 0.0 if x >= 1.0 else float(np.sqrt(1.0 - x * x))
