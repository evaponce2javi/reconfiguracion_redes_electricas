"""Monte Carlo validation of interval bounds against the linear model and the nonlinear AC power flow.

One pass over the samples measures, with the same noise draws:

(a) the linearisation error  |dV_AC - dV_lin|  per bus;
(b) the miscoverage of the covariance ellipse inside the linear model (per branch);
(c) the total violation: AC samples outside the bounds (per bus, per branch, and system-wide).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import geometry as geo
from .core import CDHIMResult
from .powerflow import BatchStats
from .study import Study
from .uncertainty import Marginal, sample_eps

VIOLATION_ATOL = 1e-8      # p.u.; an excursion below the AC solver accuracy is not counted


class LogHistogram:
    """Streaming quantiles of positive values (relative resolution about 1 %)."""

    def __init__(self, lo=-16.0, hi=4.0, bins=4000):
        self.edges = np.linspace(lo, hi, bins + 1)
        self.counts = np.zeros(bins + 2, dtype=np.int64)
        self.max = 0.0

    def add(self, values):
        v = np.asarray(values).ravel()
        if v.size == 0:
            return
        self.max = max(self.max, float(v.max()))
        with np.errstate(divide="ignore"):
            idx = np.searchsorted(self.edges, np.log10(v), side="right")
        self.counts += np.bincount(idx, minlength=self.counts.size)

    def quantile(self, q):
        total = self.counts.sum()
        if total == 0:
            return float("nan")
        k = int(np.searchsorted(np.cumsum(self.counts), q * total))
        if k == 0:
            return 0.0
        return min(self.max, float(10.0 ** self.edges[min(k, len(self.edges) - 1)]))


@dataclass
class CellStats:
    n: int = 0                       # samples drawn
    n_ac: int = 0                    # samples with a converged AC solution
    solver: BatchStats = field(default_factory=BatchStats)
    bus_lin: dict = field(default_factory=dict)      # name -> per-bus violation counts (linear model)
    bus_ac: dict = field(default_factory=dict)
    any_lin: dict = field(default_factory=dict)      # name -> samples with at least one bus outside
    any_ac: dict = field(default_factory=dict)
    br_lin: np.ndarray | None = None                 # per Ellipse-branch violation counts (Mahalanobis)
    br_ac: np.ndarray | None = None
    any_br_lin: int = 0
    any_br_ac: int = 0
    br_scaled: dict = field(default_factory=dict)    # (scale name, model) -> [per-branch counts, any-count]
    err_abs: LogHistogram = field(default_factory=LogHistogram)   # |dV_AC - dV_lin| [p.u.]
    err_rel: LogHistogram = field(default_factory=LogHistogram)   # same, over the output radius
    err_bus_max: np.ndarray | None = None
    d2_lin_max: np.ndarray | None = None             # per-branch max of d2 / chi2 (linear samples)

    # ------------------------------------------------------------------ summaries
    @staticmethod
    def _rates(counts, n, mask=None):
        r = counts / max(n, 1)
        if mask is not None:
            r = r[mask]
        if r.size == 0:
            return dict(mean=np.nan, p99=np.nan, max=np.nan, n_over=0, n=0)
        return dict(mean=float(r.mean()), p99=float(np.quantile(r, 0.99)), max=float(r.max()),
                    n_over=int((r > 0.0027).sum()), n=int(r.size))

    def bus_summary(self, name, model="ac", mask=None):
        counts, n = (self.bus_ac[name], self.n_ac) if model == "ac" else (self.bus_lin[name], self.n)
        out = self._rates(counts, n, mask)
        any_ = (self.any_ac if model == "ac" else self.any_lin)[name]
        out["system"] = any_ / max(n, 1)
        out["expected_per_sample"] = float(counts.sum() / max(n, 1))
        return out

    def branch_summary(self, model="ac", scale=None):
        n = self.n_ac if model == "ac" else self.n
        if scale is None:
            counts, any_ = (self.br_ac, self.any_br_ac) if model == "ac" else (self.br_lin, self.any_br_lin)
        else:
            counts, any_ = self.br_scaled[(scale, model)]
        out = self._rates(counts, n)
        out["system"] = any_ / max(n, 1)
        out["expected_per_sample"] = float(counts.sum() / max(n, 1)) if counts.size else 0.0
        return out


def monte_carlo(study: Study, delta: float, res: CDHIMResult, bounds: dict, n_samples: int, seed: int,
                marginal: Marginal | None = None, copula_r: float = 0.0, batch: int = 1000, ac: bool = True,
                ellipse_scales: dict | None = None) -> CellStats:
    """Draw ``n_samples`` noise vectors and accumulate all violation statistics.

    ``bounds``: name -> half-width per PQ bus; every interval is centred on the converged AC solution.
    The ellipses tested are those of the Ellipse-mode branches of ``res`` at ``res.chi2``.
    ``ellipse_scales``: name -> s; the same ellipses with covariance s * Q are tested as well (a sample
    violates them when d2 > s * chi2), e.g. s = (1/3)/Var(eps) for the uniform-variance assumption.
    """
    rng = np.random.default_rng(seed)
    n_pq = study.C1.shape[0]
    m = res.ellipse
    fi, ti = res.fi[m], res.ti[m]
    q11, q22, q12 = res.q11[m], res.q22[m], res.q12[m]
    st = CellStats(br_lin=np.zeros(int(m.sum()), np.int64), br_ac=np.zeros(int(m.sum()), np.int64),
                   err_bus_max=np.zeros(n_pq), d2_lin_max=np.zeros(int(m.sum())))
    for name in bounds:
        st.bus_lin[name] = np.zeros(n_pq, np.int64)
        st.bus_ac[name] = np.zeros(n_pq, np.int64)
        st.any_lin[name] = st.any_ac[name] = 0
    scales = ellipse_scales or {}
    for name in scales:
        for model in ("lin", "ac"):
            st.br_scaled[(name, model)] = [np.zeros(int(m.sum()), np.int64), 0]
    r_ref = np.where(res.r_out > max(res.cfg.min_radius_pu, 1e-10), res.r_out, np.inf)[:, None]   # skip round-off radii
    C = study.C1 * delta

    for start in range(0, n_samples, batch):
        b = min(batch, n_samples - start)
        eps = sample_eps(rng, study.src.n, b, marginal, copula_r)
        dv_lin = C @ eps
        st.n += b
        for name, r in bounds.items():
            out = np.abs(dv_lin) > r[:, None] + VIOLATION_ATOL
            st.bus_lin[name] += out.sum(axis=1)
            st.any_lin[name] += int(out.any(axis=0).sum())
        if fi.size:
            d2 = geo.mahalanobis2(dv_lin[fi], dv_lin[ti], q11, q22, q12)
            out = d2 > res.chi2
            st.br_lin += out.sum(axis=1)
            st.any_br_lin += int(out.any(axis=0).sum())
            np.maximum(st.d2_lin_max, d2.max(axis=1) / res.chi2, out=st.d2_lin_max)
            for name, sc in scales.items():
                out = d2 > sc * res.chi2
                st.br_scaled[(name, "lin")][0] += out.sum(axis=1)
                st.br_scaled[(name, "lin")][1] += int(out.any(axis=0).sum())
        if not ac:
            continue
        V, conv, bs = study.solver.solve(study.injections(eps, delta))
        st.solver.add(bs)
        if not conv.any():
            continue
        dv_ac = np.abs(V[study.grid.pq][:, conv]) - study.vm_star[:, None]
        st.n_ac += int(conv.sum())
        for name, r in bounds.items():
            out = np.abs(dv_ac) > r[:, None] + VIOLATION_ATOL
            st.bus_ac[name] += out.sum(axis=1)
            st.any_ac[name] += int(out.any(axis=0).sum())
        if fi.size:
            d2 = geo.mahalanobis2(dv_ac[fi], dv_ac[ti], q11, q22, q12)
            out = d2 > res.chi2
            st.br_ac += out.sum(axis=1)
            st.any_br_ac += int(out.any(axis=0).sum())
            for name, sc in scales.items():
                out = d2 > sc * res.chi2
                st.br_scaled[(name, "ac")][0] += out.sum(axis=1)
                st.br_scaled[(name, "ac")][1] += int(out.any(axis=0).sum())
        err = np.abs(dv_ac - dv_lin[:, conv])
        st.err_abs.add(err)
        st.err_rel.add(err / r_ref)
        np.maximum(st.err_bus_max, err.max(axis=1), out=st.err_bus_max)
    return st


def standard_bounds(res: CDHIMResult, res_coord: CDHIMResult | None = None) -> dict:
    """The interval families compared throughout the experiments (half-width per PQ bus)."""
    out = {"cdhim_lmax": res.r_out, "affine_box": res.r_affine, "cltm_3sigma": res.r_cltm}
    if res_coord is not None:
        out["cdhim_coord"] = res_coord.r_out
    return out
