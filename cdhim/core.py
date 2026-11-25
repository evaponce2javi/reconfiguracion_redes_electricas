"""CD-HIM: branch-wise Ellipse/Box selection and the interval bounds it returns (Algorithm 1).

Three objects are produced for a branch in Ellipse mode and they are NOT the same set:

1. the covariance ellipse        {dz : dz' Q^-1 dz <= chi2}                      (joint, rotated)
2. the circumscribed radius      sqrt(chi2 * lambda_max)       (Eq. 11, same for both terminals)
3. the output interval per bus   V0_i +- r_i                   (what Algorithm 1 returns)

``radius="coordinate"`` replaces (2) by the exact axis projection sqrt(chi2 * Q_ii).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from . import geometry as geo
from .uncertainty import SourceCovariance


@dataclass(frozen=True)
class MethodConfig:
    gamma_th: float = 0.15                 # gain threshold (fraction)
    delta_stab: float = 1e-3               # Ellipse only if |rho| < 1 - delta_stab
    confidence: float = 0.9973
    radius: str = "lambda_max"             # "lambda_max" (Eq. 11) | "coordinate" (exact projection)
    gain_box: str = "covariance_1sigma"    # box inside G: see geometry.GAIN_BOXES
    selection_rule: str = "gain"           # "gain": G >= gamma_th | "rho": |rho| >= rho_th
    rho_th: float | None = None            # threshold of the "rho" rule
    bus_combine: str = "min"               # a bus in several Ellipse branches keeps the min/max radius
    clip_to_affine_box: bool = False       # never return an interval wider than the affine box
    bonferroni_m: int | None = None        # use the chi2 quantile at level 1 - alpha/m
    min_radius_pu: float = 1e-10           # buses whose affine radius is below this are numerically
                                           # deterministic (round-off level); their branches are skipped
    min_sigma_product: float = 0.0         # original code: skip branches with sx*sy <= 1e-12 (an absolute
                                           # filter that moves with delta); 0 disables it

    def chi2(self) -> float:
        if self.bonferroni_m:
            return geo.chi2_quantile_bonferroni(self.confidence, self.bonferroni_m)
        return geo.chi2_quantile(self.confidence, 2)


@dataclass
class CDHIMResult:
    cfg: MethodConfig
    chi2: float
    # ---- per PQ bus (index = position in grid.pq)
    q_ii: np.ndarray               # variance of the linearised voltage deviation
    r_affine: np.ndarray           # affine box  sum|c_k|
    r_cltm: np.ndarray             # 3-sigma box
    r_coord: np.ndarray            # sqrt(chi2 Q_ii)
    r_out: np.ndarray              # interval half-width returned by the method
    bus_ellipse: np.ndarray        # True if the bus takes its bound from an Ellipse-mode branch
    # ---- per PQ-PQ branch
    fi: np.ndarray                 # positions of the two terminals in grid.pq
    ti: np.ndarray
    q11: np.ndarray
    q22: np.ndarray
    q12: np.ndarray
    rho: np.ndarray
    lam_min: np.ndarray
    lam_max: np.ndarray
    gain: np.ndarray               # G as a fraction
    kappa: np.ndarray              # (1 - A_ell/A_box) = 1 - kappa sqrt(1 - rho^2)
    ellipse: np.ndarray            # decision
    r_lmax: np.ndarray             # sqrt(chi2 lambda_max)
    valid: np.ndarray              # False for degenerate branches (zero variance), excluded like the original code
    extra: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ counts and CRI definitions
    @property
    def n_branches(self) -> int:
        return int(self.valid.sum())

    @property
    def n_ellipse(self) -> int:
        return int(self.ellipse.sum())

    def cri_area(self, over: str = "ellipse") -> float:
        """Manuscript CRI: mean coupling gain G (area of ellipse vs. the box named in ``gain_box``)."""
        mask = self.ellipse if over == "ellipse" else self.valid
        return float(self.gain[mask].mean()) if mask.any() else float("nan")

    def width_reduction(self, reference: str = "affine", over: str = "ellipse") -> float:
        """CRI on what the method returns: mean of 1 - r_out / r_reference over buses."""
        ref = {"affine": self.r_affine, "cltm": self.r_cltm}[reference]
        mask = self.bus_ellipse if over == "ellipse" else np.ones_like(self.bus_ellipse)
        mask = mask & (ref > 0)
        return float(np.mean(1.0 - self.r_out[mask] / ref[mask])) if mask.any() else float("nan")

    def output_area_reduction(self, reference: str = "affine") -> float:
        """Mean over Ellipse branches of 1 - (output rectangle area) / (reference box area)."""
        ref = {"affine": self.r_affine, "cltm": self.r_cltm}[reference]
        m = self.ellipse
        if not m.any():
            return float("nan")
        out = self.r_out[self.fi[m]] * self.r_out[self.ti[m]]
        return float(np.mean(1.0 - out / (ref[self.fi[m]] * ref[self.ti[m]])))

    def overconservatism(self):
        """sqrt(lambda_max / Q_ii) for every (Ellipse branch, terminal) pair."""
        m = self.ellipse
        return np.concatenate([np.sqrt(self.lam_max[m] / self.q11[m]), np.sqrt(self.lam_max[m] / self.q22[m])])


def branch_covariances(C, fi, ti, sigma_eps: SourceCovariance | None = None, chunk: int = 2048):
    """Q_ii for every row of ``C`` and Q_ij for every branch (fi, ti), with Q = C Sigma_eps C^T."""
    sigma_eps = sigma_eps or SourceCovariance()
    structured = sigma_eps.full is None and not sigma_eps.rho_c
    if structured:                                   # diagonal Sigma_eps: no (n_bus x N_u) temporary
        var = np.broadcast_to(np.asarray(sigma_eps.var, float), (C.shape[1],))
        q_ii = np.einsum("ij,ij,j->i", C, C, var)
        q_ij = np.empty(len(fi))
        for s in range(0, len(fi), chunk):
            sl = slice(s, s + chunk)
            q_ij[sl] = np.einsum("ij,ij,j->i", C[fi[sl]], C[ti[sl]], var)
        return q_ii, q_ij
    M = sigma_eps.apply(C)
    q_ii = np.einsum("ij,ij->i", M, C)
    q_ij = np.empty(len(fi))
    for s in range(0, len(fi), chunk):
        sl = slice(s, s + chunk)
        q_ij[sl] = np.einsum("ij,ij->i", M[fi[sl]], C[ti[sl]])
    return q_ii, q_ij


def run_cdhim(C, fi, ti, cfg: MethodConfig = MethodConfig(), sigma_eps: SourceCovariance | None = None,
              r_affine=None) -> CDHIMResult:
    """Run the hybrid selection and build the output intervals.

    ``C``: (n_pq, N_u) affine coefficients of the voltage magnitudes of the PQ buses.
    ``fi``, ``ti``: positions (rows of ``C``) of the terminals of every PQ-PQ branch.
    """
    chi2 = cfg.chi2()
    q_ii, q12 = branch_covariances(C, fi, ti, sigma_eps)
    if r_affine is None:
        r_affine = geo.affine_box_halfwidth(C)
    q11, q22 = q_ii[fi], q_ii[ti]

    # branches touching a bus without (numerically meaningful) uncertainty are not analysed
    valid = ((r_affine[fi] > cfg.min_radius_pu) & (r_affine[ti] > cfg.min_radius_pu)
             & (np.sqrt(q11 * q22) > cfg.min_sigma_product))
    rho = np.where(valid, geo.correlation(q11, q22, q12), 0.0)
    lam_min, lam_max = geo.eig2x2(q11, q22, q12)

    if cfg.gain_box == "covariance_1sigma":          # ellipse and box at the same (1-sigma) scale
        a_ell = geo.ellipse_area(lam_min, lam_max, 1.0)
        a_box = geo.box_area(geo.cov_box_halfwidth(q11, 1.0), geo.cov_box_halfwidth(q22, 1.0))
    elif cfg.gain_box == "covariance_3sigma":        # chi2 ellipse against the 3-sigma (CLTM) box
        a_ell = geo.ellipse_area(lam_min, lam_max, chi2)
        a_box = geo.box_area(geo.cltm_box_halfwidth(q11), geo.cltm_box_halfwidth(q22))
    elif cfg.gain_box == "affine":                   # chi2 ellipse against the deterministic affine box
        a_ell = geo.ellipse_area(lam_min, lam_max, chi2)
        a_box = geo.box_area(r_affine[fi], r_affine[ti])
    else:
        raise ValueError(f"unknown gain_box {cfg.gain_box!r}")
    gain = np.where(valid, geo.coupling_gain(a_ell, a_box), 0.0)
    root = np.sqrt(np.maximum(1.0 - rho**2, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = np.where(valid & (root > 0), (a_ell / np.where(a_box > 0, a_box, 1.0)) / root, np.nan)

    stable = np.abs(rho) < 1.0 - cfg.delta_stab
    if cfg.selection_rule == "gain":
        ellipse = valid & (gain >= cfg.gamma_th) & stable
    elif cfg.selection_rule == "rho":
        if cfg.rho_th is None:
            raise ValueError("selection_rule='rho' needs rho_th")
        ellipse = valid & (np.abs(rho) >= cfg.rho_th) & stable
    else:
        raise ValueError(f"unknown selection_rule {cfg.selection_rule!r}")

    r_lmax = geo.radius_lambda_max(lam_max, chi2)
    r_coord = geo.radius_coordinate(q_ii, chi2)

    # ---- per-bus output interval
    n = C.shape[0]
    bus_ellipse = np.zeros(n, dtype=bool)
    bus_ellipse[fi[ellipse]] = True
    bus_ellipse[ti[ellipse]] = True
    if cfg.radius == "lambda_max":
        if cfg.bus_combine == "min":
            acc = np.full(n, np.inf)
            np.minimum.at(acc, fi[ellipse], r_lmax[ellipse])
            np.minimum.at(acc, ti[ellipse], r_lmax[ellipse])
        elif cfg.bus_combine == "max":
            acc = np.zeros(n)
            np.maximum.at(acc, fi[ellipse], r_lmax[ellipse])
            np.maximum.at(acc, ti[ellipse], r_lmax[ellipse])
        else:
            raise ValueError("bus_combine must be 'min' or 'max'")
        r_out = np.where(bus_ellipse, acc, r_affine)
    elif cfg.radius == "coordinate":
        r_out = np.where(bus_ellipse, r_coord, r_affine)
    else:
        raise ValueError("radius must be 'lambda_max' or 'coordinate'")
    if cfg.clip_to_affine_box:
        r_out = np.minimum(r_out, r_affine)

    return CDHIMResult(cfg, chi2, q_ii, r_affine, geo.cltm_box_halfwidth(q_ii), r_coord, r_out, bus_ellipse,
                       np.asarray(fi), np.asarray(ti), q11, q22, q12, rho, lam_min, lam_max, gain, kappa,
                       ellipse, r_lmax, valid)


def with_method(res_cfg: MethodConfig, **changes) -> MethodConfig:
    return replace(res_cfg, **changes)
