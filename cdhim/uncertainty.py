"""Noise-symbol model: marginal distributions on the affine support, source correlation, sampling.

A noise symbol ``eps_k`` multiplies the half-width of one uncertain injection. Affine arithmetic only
says ``eps_k in [-1, 1]``; treating the symbols as random variables with covariance ``Sigma_eps`` is an
additional modelling assumption, and it is made explicit here:  Q = C Sigma_eps C^T.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special, stats

UNIFORM_VARIANCE = 1.0 / 3.0


# --------------------------------------------------------------------------- marginals
@dataclass(frozen=True)
class Marginal:
    """A zero-mean marginal: eps = (X - shift) * scale, with X ~ ``base`` restricted to [x_lo, x_hi]."""

    name: str
    base: object                 # frozen scipy.stats distribution of X
    shift: float
    scale: float
    x_lo: float = -np.inf        # truncation of X (conditional distribution, not clipping)
    x_hi: float = np.inf
    note: str = ""

    # -- exact properties ---------------------------------------------------
    @property
    def _p_lo(self) -> float:
        return float(self.base.cdf(self.x_lo))

    @property
    def _p_hi(self) -> float:
        return float(self.base.cdf(self.x_hi))

    @property
    def support(self) -> tuple[float, float]:
        lo = max(self.x_lo, float(self.base.support()[0]))
        hi = min(self.x_hi, float(self.base.support()[1]))
        return (lo - self.shift) * self.scale, (hi - self.shift) * self.scale

    def _x_moment(self, order: int) -> float:
        lo, hi = max(self.x_lo, float(self.base.support()[0])), min(self.x_hi, float(self.base.support()[1]))
        if np.isinf(lo) and np.isinf(hi) or (self._p_lo == 0.0 and self._p_hi == 1.0):
            return float(self.base.moment(order))
        return float(self.base.expect(lambda x: x**order, lb=lo, ub=hi, conditional=True))

    @property
    def in_model(self) -> bool:
        """True when the support lies inside the affine interval [-1, 1]."""
        lo, hi = self.support
        return bool(lo >= -1.0 - 1e-9 and hi <= 1.0 + 1e-9)

    @property
    def mean(self) -> float:
        return (self._x_moment(1) - self.shift) * self.scale

    @property
    def variance(self) -> float:
        m1, m2 = self._x_moment(1), self._x_moment(2)
        return (m2 - m1 * m1) * self.scale**2

    @property
    def fraction_outside_unit(self) -> float:
        """P(|eps| > 1): mass that leaves the affine support [-1, 1]."""
        span = self._p_hi - self._p_lo
        p_in = self.base.cdf(np.clip(self.shift + 1.0 / self.scale, self.x_lo, self.x_hi)) - self.base.cdf(
            np.clip(self.shift - 1.0 / self.scale, self.x_lo, self.x_hi)
        )
        return float(max(0.0, 1.0 - p_in / span))

    # -- sampling -------------------------------------------------------------
    def from_uniform(self, u):
        """Inverse-CDF transform of U(0,1) variates (used with the Gaussian copula)."""
        p = self._p_lo + u * (self._p_hi - self._p_lo)
        return (self.base.ppf(p) - self.shift) * self.scale

    def sample(self, rng: np.random.Generator, size):
        """Independent draws. Uses the generator's native sampler (much faster than the inverse CDF)
        and rejection for truncated marginals."""
        x = np.asarray(self.base.rvs(size=size, random_state=rng), dtype=float)
        if np.isfinite(self.x_lo) or np.isfinite(self.x_hi):
            bad = (x < self.x_lo) | (x > self.x_hi)
            while bad.any():
                x[bad] = self.base.rvs(size=int(bad.sum()), random_state=rng)
                bad = (x < self.x_lo) | (x > self.x_hi)
        return (x - self.shift) * self.scale


def _centred(name, base, scale=None, target_var=None, x_lo=-np.inf, x_hi=np.inf, note=""):
    """Centre on the (possibly truncated) mean; scale to respect the support or to hit a variance."""
    tmp = Marginal(name, base, 0.0, 1.0, x_lo, x_hi)
    mu, var = tmp.mean, tmp.variance
    if target_var is not None:
        scale = float(np.sqrt(target_var / var))
    elif scale is None:                       # largest scale that keeps eps inside [-1, 1]
        lo, hi = tmp.support
        scale = 1.0 / max(mu - lo, hi - mu)
    return Marginal(name, base, mu, scale, x_lo, x_hi, note)


def marginal_catalogue() -> dict[str, Marginal]:
    """All marginals used in experiment C3, keyed by name."""
    b25, b52, bu = stats.beta(2, 5), stats.beta(5, 2), stats.beta(0.5, 0.5)
    t3, ln = stats.t(3), stats.lognorm(s=0.5)
    cat = [
        Marginal("uniform", stats.uniform(-1, 2), 0.0, 1.0, note="U[-1,1], the manuscript's assumption"),
        _centred("beta25_support", b25, note="Beta(2,5): (u-mu)/max(mu,1-mu), support respected"),
        _centred("beta25_var13", b25, target_var=UNIFORM_VARIANCE,
                 note="Beta(2,5) rescaled to variance 1/3; leaves [-1,1]"),
        Marginal("beta25_paper", b25, 0.285, 4.0, note="(u-0.285)*4 exactly as in finalFinalCode.py:580"),
        _centred("beta52_support", b52, note="Beta(5,2): opposite skew, support respected"),
        _centred("beta0505_ushape", bu, note="Beta(0.5,0.5): U-shaped, variance 1/2 > 1/3"),
        Marginal("triangular", stats.triang(c=0.5, loc=-1, scale=2), 0.0, 1.0, note="symmetric triangular on [-1,1]"),
        Marginal("truncnorm_3s", stats.norm(0, 1 / 3), 0.0, 1.0, -1.0, 1.0, note="N(0,1/3^2) truncated at +-3 sigma"),
        _centred("student_t3_trunc", t3, x_lo=float(t3.ppf(0.00135)), x_hi=float(t3.ppf(0.99865)),
                 note="Student t (3 dof) truncated to its central 99.73 %, rescaled into [-1,1]"),
        _centred("student_t3", t3, target_var=UNIFORM_VARIANCE,
                 note="Student t (3 dof), variance 1/3, unbounded: outside the model"),
        _centred("lognormal_trunc", ln, x_lo=0.0, x_hi=float(ln.ppf(0.9973)),
                 note="lognormal(s=0.5) truncated at its 99.73 % quantile, rescaled into [-1,1]"),
        _centred("lognormal", ln, target_var=UNIFORM_VARIANCE,
                 note="lognormal(s=0.5), variance 1/3, unbounded above: outside the model"),
    ]
    return {m.name: m for m in cat}


# --------------------------------------------------------------------------- source correlation
def copula_to_pearson_uniform(r):
    """Pearson correlation of two U[-1,1] variables joined by a Gaussian copula of parameter r."""
    return 6.0 / np.pi * np.arcsin(np.asarray(r, float) / 2.0)


def copula_to_pearson(marginal: Marginal, r: float, nodes: int = 96) -> float:
    """Pearson correlation of two ``marginal`` variables joined by a Gaussian copula of parameter r
    (Gauss-Hermite quadrature of E[g(z1) g(z2)], g = marginal quantile of Phi(z))."""
    if r == 0.0:
        return 0.0
    x, w = np.polynomial.hermite_e.hermegauss(nodes)
    w = w / np.sqrt(2.0 * np.pi)
    g = lambda z: marginal.from_uniform(np.clip(special.ndtr(z), 1e-15, 1 - 1e-15))
    z2 = r * x[:, None] + np.sqrt(1.0 - r * r) * x[None, :]
    e12 = float(np.sum(w[:, None] * w[None, :] * g(x)[:, None] * g(z2)))
    m1 = float(np.sum(w * g(x)))
    v1 = float(np.sum(w * g(x) ** 2)) - m1 * m1
    return (e12 - m1 * m1) / v1


@dataclass(frozen=True)
class SourceCovariance:
    """Sigma_eps = diag(d) + rho_c * (s s^T - diag(s^2)) -- diagonal plus one common factor.

    ``d`` holds the variances and ``s`` the standard deviations of the symbols, ``rho_c`` their common
    (equi-)correlation. ``rho_c = 0`` with ``d = 1/3`` is the manuscript's model, Sigma_eps = I/3.
    A full matrix can be supplied instead through :meth:`dense`.
    """

    var: np.ndarray | float = UNIFORM_VARIANCE
    rho_c: float = 0.0
    full: np.ndarray | None = None

    @staticmethod
    def dense(sigma: np.ndarray) -> "SourceCovariance":
        sigma = np.asarray(sigma, float)
        if sigma.shape[0] != sigma.shape[1] or not np.allclose(sigma, sigma.T, atol=1e-12):
            raise ValueError("Sigma_eps must be square and symmetric")
        return SourceCovariance(full=sigma)

    def apply(self, C: np.ndarray) -> np.ndarray:
        """Return C @ Sigma_eps without forming Sigma_eps when it is structured."""
        if self.full is not None:
            return C @ self.full
        var = np.broadcast_to(np.asarray(self.var, float), (C.shape[1],))
        out = C * var
        if self.rho_c:
            s = np.sqrt(var)
            out = out * (1.0 - self.rho_c) + self.rho_c * np.outer(C @ s, s)
        return out

    def matrix(self, n: int) -> np.ndarray:
        return self.apply(np.eye(n))


def estimate_source_covariance(series: np.ndarray, centre: str = "midrange"):
    """Affine half-widths and Sigma_eps from measured injections (rows = time steps, columns = sources).

    Each column is mapped to the affine form x = x0 + h * eps with eps in [-1, 1]:
    ``centre="midrange"`` uses x0 = (max + min)/2 and h = (max - min)/2 (support exactly [-1, 1], the
    symbols may then have non-zero mean); ``centre="mean"`` uses x0 = mean and h = max |x - mean|.
    Returns ``(x0, h, mean_eps, SourceCovariance)``; the covariance is taken about ``mean_eps``.
    """
    x = np.asarray(series, float)
    if x.ndim != 2 or x.shape[0] < 2:
        raise ValueError("series must be (n_time, n_sources) with at least two time steps")
    lo, hi = x.min(axis=0), x.max(axis=0)
    if centre == "midrange":
        x0, h = 0.5 * (hi + lo), 0.5 * (hi - lo)
    elif centre == "mean":
        x0 = x.mean(axis=0)
        h = np.max(np.abs(x - x0), axis=0)
    else:
        raise ValueError("centre must be 'midrange' or 'mean'")
    h = np.where(h > 0, h, 1.0)
    eps = (x - x0) / h
    return x0, h, eps.mean(axis=0), SourceCovariance.dense(np.cov(eps, rowvar=False))


# --------------------------------------------------------------------------- sampling
def sample_eps(rng: np.random.Generator, n_sources: int, n_samples: int, marginal: Marginal | None = None,
               copula_r: float = 0.0) -> np.ndarray:
    """Draw eps with shape (n_sources, n_samples).

    ``copula_r`` > 0 couples all sources through a one-factor Gaussian copula with parameter r
    (z_k = sqrt(r) z_0 + sqrt(1-r) e_k); the marginals are unchanged. Each sample consumes a fixed
    block of the random stream, so the result does not depend on how the samples are batched.
    """
    uniform = marginal is None or marginal.name == "uniform"
    if copula_r <= 0.0:
        if uniform:
            return 2.0 * rng.random((n_samples, n_sources)).T - 1.0
        return np.ascontiguousarray(marginal.sample(rng, (n_samples, n_sources)).T)
    g = rng.standard_normal((n_samples, n_sources + 1)).T
    u = special.ndtr(np.sqrt(1.0 - copula_r) * g[1:] + np.sqrt(copula_r) * g[:1])
    return 2.0 * u - 1.0 if uniform else marginal.from_uniform(u)
