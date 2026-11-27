"""Second-order correction of the CD-HIM bounds for the linearisation error.

The voltage deviation of every PQ bus is, to second order in delta,
    d|V_i| = c_i' eps + q_i(eps),        q_i(eps) = eps' H_i eps        (cdhim.second_order.SecondOrder),
and q_i is what separates the AC samples from the linear model at the design point (RESULTS.md, B1).
For independent, zero-mean, symmetric symbols the linear and quadratic parts are uncorrelated
(E[eps_k eps_l eps_m] = 0), so the first two moments of the second-order model are
    mean:      mu_i = E[q_i] = Var(eps) tr(H_i)              (exact, no H_i is formed)
    variance:  Q_ii + Var(q_i)
    branch:    Q_ij + Cov(q_i, q_j)
Var(q_i) and Cov(q_i, q_j) are estimated from samples of the quadratic surrogate (cheap: no power flow).

Three corrections are built from them, all relative to the converged AC solution:
    shift      every interval and ellipse is re-centred on mu (radius unchanged);
    moment     re-centred on mu and every probabilistic radius scaled by sqrt(1 + Var(q_i) / Q_ii);
               every ellipse uses the covariance Q + Cov(q);
    box2       the deterministic affine box becomes [-r_aff + lo(q_i), r_aff + hi(q_i)], with
               [lo, hi] an enclosure of q_i on the unit box (SecondOrder.term_bounds, "binomial");
               rigorous for the second-order model, but loose when many symbols act on the bus;
    box2v      the affine box with q evaluated at the two vertices where the linear part is extreme,
               eps = +-sign(c_i): [-r_aff + min(0, q_i(-s_i)), r_aff + max(0, q_i(s_i))];
    minkowski  every ellipse replaced by an outer ellipse of  E(Q) + B,  where B is the box spanned by the
               sampled range of (q_i, q_j) projected on the principal axes of Q; the quadratic term along the minor axis of a thin ellipse is
               one-sided and heavy-tailed, which a moment correction cannot represent.
The third-order remainder is not bounded: the AC check measures what is left.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import CDHIMResult
from .second_order import SecondOrder
from .uncertainty import sample_eps

VAR_UNIFORM = 1.0 / 3.0


@dataclass
class Correction:
    mu: np.ndarray            # E[q_i] per PQ bus
    var_q: np.ndarray         # Var(q_i) per PQ bus (sampled)
    cov_q: np.ndarray         # Cov(q_fi, q_ti) per branch of the result (sampled)
    q_lo: np.ndarray          # enclosure of q_i on [-1, 1]^N_u
    q_hi: np.ndarray
    q_min: np.ndarray         # sampled range of q_i
    q_max: np.ndarray
    q_vlo: np.ndarray         # q_i at eps = -sign(c_i) and at eps = +sign(c_i)
    q_vhi: np.ndarray
    pu_min: np.ndarray        # sampled range of (q_fi, q_ti) projected on the principal axes of each branch
    pu_max: np.ndarray
    pv_min: np.ndarray
    pv_max: np.ndarray
    n_samples: int
    seed: int
    mean_check: float         # max |sampled mean - mu|, a consistency check of the exact mean
    so: SecondOrder

    def scale(self, res: CDHIMResult) -> np.ndarray:
        """Radius factor sqrt(1 + Var(q_i)/Q_ii) (1 where the linear variance is zero)."""
        q = res.q_ii
        return np.sqrt(1.0 + np.divide(self.var_q, q, out=np.zeros_like(q), where=q > 0))

    def interval(self, r: np.ndarray, res: CDHIMResult, kind: str):
        """(lo, hi) offsets from the AC solution for a symmetric linear radius ``r``."""
        if kind == "none":
            return -r, r
        if kind == "shift":
            return self.mu - r, self.mu + r
        if kind == "moment":
            s = self.scale(res) * r
            return self.mu - s, self.mu + s
        if kind == "box2":
            return -r + self.q_lo, r + self.q_hi
        if kind == "box2v":
            return -r + np.minimum(0.0, self.q_vlo), r + np.maximum(0.0, self.q_vhi)
        if kind in ("minkowski", "range"):
            return -r + np.minimum(0.0, self.q_min), r + np.maximum(0.0, self.q_max)
        raise ValueError(kind)

    def cdhim(self, res: CDHIMResult, kind: str, box: str = "box2v"):
        """Output interval of Algorithm 1 with the correction: Ellipse buses take ``kind``,
        Box buses take the second-order affine box ``box``."""
        lo, hi = self.interval(res.r_out, res, kind)
        blo, bhi = self.interval(res.r_affine, res, box)
        box = ~res.bus_ellipse
        return np.where(box, blo, lo), np.where(box, bhi, hi)

    def ellipse(self, res: CDHIMResult, kind: str):
        """(centre_f, centre_t, q11, q22, q12) of every Ellipse-mode branch."""
        m = res.ellipse
        f, t = res.fi[m], res.ti[m]
        if kind == "none":
            return np.zeros(f.size), np.zeros(f.size), res.q11[m], res.q22[m], res.q12[m]
        cf, ct = self.mu[f], self.mu[t]
        if kind == "shift":
            return cf, ct, res.q11[m], res.q22[m], res.q12[m]
        if kind == "moment":
            return cf, ct, res.q11[m] + self.var_q[f], res.q22[m] + self.var_q[t], res.q12[m] + self.cov_q[m]
        if kind == "minkowski":
            # Principal axes of Q: u (major) and v (minor), squared semi-axes a2 = chi2 lam_max, b2 = chi2 lam_min.
            # The sampled range of q projected on (u, v) is a box with centre (cu, cv) and half-widths (hu, hv),
            # contained in the aligned ellipse with squared semi-axes (2 hu^2, 2 hv^2). For any p > 0,
            #     E(Q1) + E(Q2)  is inside  E((1 + 1/p) Q1 + (1 + p) Q2),
            # and with aligned ellipses the result stays aligned; p is chosen per branch to minimise the area.
            chi2 = res.chi2
            ux, uy, vx, vy = principal_axes(res.q11[m], res.q22[m], res.q12[m])
            cu, cv = 0.5 * (self.pu_min[m] + self.pu_max[m]), 0.5 * (self.pv_min[m] + self.pv_max[m])
            hu, hv = 0.5 * (self.pu_max[m] - self.pu_min[m]), 0.5 * (self.pv_max[m] - self.pv_min[m])
            a2, b2 = chi2 * res.lam_max[m], chi2 * np.maximum(res.lam_min[m], 0.0)
            A2, B2 = outer_sum_axes(a2, b2, 2.0 * hu ** 2, 2.0 * hv ** 2)
            l1, l2 = A2 / chi2, B2 / chi2
            return (cu * ux + cv * vx, cu * uy + cv * vy,
                    l1 * ux * ux + l2 * vx * vx, l1 * uy * uy + l2 * vy * vy, l1 * ux * uy + l2 * vx * vy)
        raise ValueError(kind)


def second_order_correction(st, delta: float, res: CDHIMResult, n_samples: int = 4000, seed: int = 0,
                            batch: int = 500, var_eps: float = VAR_UNIFORM, vertex_buses=None,
                            enclosure: bool = True) -> Correction:
    """Moments, sampled range and enclosures of the quadratic term at every PQ bus, for the branches of ``res``.
    ``vertex_buses``: mask of the PQ buses where the vertex values are computed (default: all).
    ``enclosure``: also compute the rigorous (loose) enclosure used by "box2"; not needed by the method."""
    so = SecondOrder(st, delta)
    mu = var_eps * so.trace()
    rng = np.random.default_rng(seed)
    n_pq = so.c_lin.shape[0]
    s1, s2 = np.zeros(n_pq), np.zeros(n_pq)
    qmin, qmax = np.full(n_pq, np.inf), np.full(n_pq, -np.inf)
    sft = np.zeros(res.fi.size)
    ux, uy, vx, vy = principal_axes(res.q11, res.q22, res.q12)
    nb = res.fi.size
    pu_min, pv_min = np.full(nb, np.inf), np.full(nb, np.inf)
    pu_max, pv_max = np.full(nb, -np.inf), np.full(nb, -np.inf)
    for start in range(0, n_samples, batch):
        b = min(batch, n_samples - start)
        q = so.quadratic(sample_eps(rng, st.src.n, b))
        s1 += q.sum(axis=1)
        s2 += (q * q).sum(axis=1)
        sft += (q[res.fi] * q[res.ti]).sum(axis=1)
        np.minimum(qmin, q.min(axis=1), out=qmin)
        np.maximum(qmax, q.max(axis=1), out=qmax)
        qf, qt = q[res.fi], q[res.ti]
        pu = ux[:, None] * qf + uy[:, None] * qt
        pv = vx[:, None] * qf + vy[:, None] * qt
        np.minimum(pu_min, pu.min(axis=1), out=pu_min)
        np.maximum(pu_max, pu.max(axis=1), out=pu_max)
        np.minimum(pv_min, pv.min(axis=1), out=pv_min)
        np.maximum(pv_max, pv.max(axis=1), out=pv_max)
    m1 = s1 / n_samples
    var_q = np.maximum(s2 / n_samples - m1 * m1, 0.0) * n_samples / (n_samples - 1)
    cov_q = (sft / n_samples - m1[res.fi] * m1[res.ti]) * n_samples / (n_samples - 1)
    lo, hi = so.term_bounds()["binomial"] if enclosure else (np.full(n_pq, -np.inf), np.full(n_pq, np.inf))
    vlo, vhi = vertex_values(so, vertex_buses)
    return Correction(mu, var_q, cov_q, lo, hi, qmin, qmax, vlo, vhi, pu_min, pu_max, pv_min, pv_max, n_samples, seed, float(np.abs(m1 - mu).max()), so)


def outer_sum_axes(a2, b2, c2, d2, n_grid: int = 400):
    """Squared semi-axes of an outer ellipse of  E(a, b) + E(c, d)  (aligned ellipses):
    ((1 + 1/p) a2 + (1 + p) c2, (1 + 1/p) b2 + (1 + p) d2) with the p of minimum area on a log grid."""
    a2, b2, c2, d2 = (np.asarray(x, float)[:, None] for x in (a2, b2, c2, d2))
    p = np.logspace(-8, 8, n_grid)[None, :]
    A2 = (1.0 + 1.0 / p) * a2 + (1.0 + p) * c2
    B2 = (1.0 + 1.0 / p) * b2 + (1.0 + p) * d2
    k = np.argmin(A2 * B2, axis=1)
    r = np.arange(A2.shape[0])
    none = (c2[:, 0] == 0) & (d2[:, 0] == 0)                 # nothing to add: keep the ellipse
    return np.where(none, a2[:, 0], A2[r, k]), np.where(none, b2[:, 0], B2[r, k])


def principal_axes(q11, q22, q12):
    """Unit vectors of the major (u) and minor (v) axes of the 2x2 covariances [[q11, q12], [q12, q22]]."""
    th = 0.5 * np.arctan2(2.0 * q12, q11 - q22)
    ux, uy = np.cos(th), np.sin(th)
    return ux, uy, -uy, ux


def vertex_values(so: SecondOrder, buses=None, chunk: int = 256):
    """q_i at eps = -sign(c_i) and eps = +sign(c_i), the vertices where the linear part of bus i is extreme.
    ``buses``: positions in grid.pq to evaluate (default all); the others get 0."""
    n_pq = so.c_lin.shape[0]
    idx = np.arange(n_pq) if buses is None else np.flatnonzero(buses)
    vlo, vhi = np.zeros(n_pq), np.zeros(n_pq)
    for s in range(0, idx.size, chunk):
        b = idx[s:s + chunk]
        sgn = np.sign(so.c_lin[b]).T                       # (N_u, chunk)
        cols = np.arange(b.size)
        vhi[b] = so.quadratic(sgn)[b, cols]
        vlo[b] = so.quadratic(-sgn)[b, cols]
    return vlo, vhi
