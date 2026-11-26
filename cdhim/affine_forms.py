"""Outputs built directly from the affine forms, without the ellipse: per-bus intervals and branch quantities.

A *target* is a linear combination of PQ-bus voltage magnitudes, t = sum_i w_i |V_i|: a single bus (w = e_i)
or the voltage drop of a branch (w = e_i - e_j). Its first- and second-order model is
    dt = c_t' eps + eps' H_t eps,     c_t = sum_i w_i c_i,   H_t = sum_i w_i H_i,
so the dependence between buses is kept exactly by the shared noise symbols: no ellipse or decision rule.

Two kinds of half-width are given for the linear part:
    affine      sum_k |c_tk|  (worst case over the noise box; deterministic for the linear model);
    bernstein   Bernstein's inequality for a sum of independent, zero-mean, bounded terms |c_tk eps_k| <= |c_tk|
                with variance sigma_t^2:  P(|c_t' eps| >= r) <= alpha  for
                r = M L / 3 + sqrt((M L / 3)^2 + 2 sigma_t^2 L),  L = ln(2 / alpha),  M = max_k |c_tk|;
                capped by the affine radius. Distribution-free given independence and the support [-1, 1].

Worst vertices: the extreme of the second-order model over the vertices of the noise box is searched from the
first-order vertex sign(c_t) by flips guided by the exact gain of each single flip,
    f(x - 2 x_k e_k) - f(x) = -2 x_k (c_tk + 2 (H_t x)_k) + 4 (H_t)_kk,
with H_t x and diag(H_t) computed from the factors of cdhim.second_order (no H_t is formed).
The result is a point that is attained, not a bound: it is used to test bounds with the exact AC power flow.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .second_order import SecondOrder

VAR_UNIFORM = 1.0 / 3.0


def bus_targets(n_pq: int) -> sp.csc_matrix:
    return sp.identity(n_pq, format="csc")


def drop_targets(n_pq: int, fi, ti) -> sp.csc_matrix:
    """Column b: |V_fi| - |V_ti|."""
    nb = len(fi)
    rows = np.concatenate([fi, ti])
    cols = np.concatenate([np.arange(nb), np.arange(nb)])
    vals = np.concatenate([np.ones(nb), -np.ones(nb)])
    return sp.csc_matrix((vals, (rows, cols)), shape=(n_pq, nb))


def target_coefficients(C: np.ndarray, W) -> np.ndarray:
    """(n_target, N_u) coefficients of the linear part."""
    return np.asarray((W.T @ C))


def affine_radius(Ct: np.ndarray) -> np.ndarray:
    return np.abs(Ct).sum(axis=1)


def bernstein_radius(Ct: np.ndarray, alpha: float, var_eps: float = VAR_UNIFORM) -> np.ndarray:
    L = np.log(2.0 / alpha)
    M = np.abs(Ct).max(axis=1)
    s2 = var_eps * (Ct ** 2).sum(axis=1)
    a = M * L / 3.0
    return np.minimum(a + np.sqrt(a * a + 2.0 * s2 * L), affine_radius(Ct))


def sigma(Ct: np.ndarray, var_eps: float = VAR_UNIFORM) -> np.ndarray:
    return np.sqrt(var_eps * (Ct ** 2).sum(axis=1))


# ----------------------------------------------------------------------------- second-order model of a target
class TargetModel:
    """Second-order model of the targets W (n_pq x T) built on a cdhim.second_order.SecondOrder."""

    def __init__(self, so: SecondOrder, W):
        self.so = so
        self.W = sp.csc_matrix(W)
        self.c = np.asarray(self.W.T @ so.c_lin)                    # (T, N_u)
        self.YH = so.Y.conj().T.tocsr()
        self._diag = None

    def _cols(self, idx):
        return self.W[:, idx].toarray()                             # (n_pq, B)

    def value(self, X, idx):
        """f_t(x) = c_t'x + x'H_t x for target idx[b] at column b of X."""
        Wb = self._cols(idx)
        q = self.so.quadratic(X)
        return np.einsum("tk,kt->t", self.c[idx], X) + (Wb * q).sum(axis=0)

    def hx(self, X, idx):
        """H_t x for target idx[b] at column b of X; shape (N_u, B)."""
        so, Wb = self.so, self._cols(idx)
        G = so.gam @ Wb                                             # (n, B) weight of q_m in each target
        BV = so.bv @ Wb if so.has_pv else None
        A = so.E @ X                                                # dV of every column
        Z = 0.5 * (self.YH @ (G * A) + np.conj(G) * (so.Y @ A))
        if BV is not None:
            Z += BV * A
        out = -np.real(so.E.conj().T @ Z)
        P = so.c_perp @ X                                           # (n_pq, B)
        out += so.c_perp.T @ (Wb * P / (2.0 * so.vm[:, None]))
        return out

    def hdiag(self):
        """diag(H_t) for every target, shape (T, N_u)."""
        if self._diag is None:
            so = self.so
            T = so.E * np.conj(so.F)                                # (n, N_u): dV_m conj(dI_m) per symbol
            D = -(so.gam.T @ T).real                                # (n_pq, N_u)
            if so.has_pv:
                D -= so.bv.T @ (np.abs(so.E) ** 2)
            D += so.c_perp ** 2 / (2.0 * so.vm[:, None])
            self._diag = np.asarray(self.W.T @ D)
        return self._diag

    def worst_vertices(self, side: int, iters: int = 8, chunk: int = 128):
        """Vertex maximising side * f_t for every target (side = +1 maximum, -1 minimum), by flip ascent from
        side * sign(c_t). Returns (X (N_u, T), f_t at X, f_t at the first-order vertex, flips per target)."""
        so = self.so
        nt, nu = self.c.shape
        Xout = np.empty((nu, nt))
        f_out, f_lin, n_flip = np.empty(nt), np.empty(nt), np.zeros(nt, int)
        Hd = self.hdiag()
        for s in range(0, nt, chunk):
            idx = np.arange(s, min(s + chunk, nt))
            X = side * np.sign(self.c[idx]).T
            X[X == 0] = side
            f = side * self.value(X, idx)
            f_lin[idx] = side * f
            for _ in range(iters):
                gain = side * (-2.0 * X * (self.c[idx].T + 2.0 * self.hx(X, idx)) + 4.0 * Hd[idx].T)
                pos = gain > 1e-15 * np.maximum(np.abs(f), 1e-12)[None, :]
                if not pos.any():
                    break
                # all positive single-flip gains at once; fall back to the best single flip where that loses
                Xa = np.where(pos, -X, X)
                fa = side * self.value(Xa, idx)
                best = np.argmax(gain, axis=0)
                cols = np.arange(len(idx))
                Xb = X.copy()
                Xb[best, cols] = -Xb[best, cols]
                fb = side * self.value(Xb, idx)
                use_a = fa >= fb
                Xn = np.where(use_a[None, :], Xa, Xb)
                fn = np.where(use_a, fa, fb)
                better = (fn > f) & pos.any(axis=0)
                if not better.any():
                    break
                n_flip[idx[better]] += np.where(use_a, pos.sum(axis=0), 1)[better]
                X[:, better] = Xn[:, better]
                f = np.where(better, fn, f)
            Xout[:, idx] = X
            f_out[idx] = side * f
        return Xout, f_out, f_lin, n_flip


# ----------------------------------------------------------------------------- correlated sources as shared symbols
class FactorModel:
    """Correlated sources written with independent bounded symbols (the affine-arithmetic way):
        eps_k = w * zeta_{z(k)} + (1 - w) * eta_k,     zeta_z, eta_k independent U[-1, 1],
    so eps_k stays in [-1, 1] and two sources of the same zone have correlation w^2 / (w^2 + (1 - w)^2).
    ``zone``: zone index of every source (all zeros: one global factor)."""

    def __init__(self, zone, w: float):
        self.zone = np.asarray(zone, int)
        self.k = int(self.zone.max()) + 1
        self.w = float(w)
        self.Z = sp.csr_matrix((np.ones(self.zone.size), (np.arange(self.zone.size), self.zone)), shape=(self.zone.size, self.k))

    @property
    def rho(self) -> float:
        w = self.w
        return w * w / (w * w + (1 - w) ** 2)

    def coefficients(self, Ct: np.ndarray) -> np.ndarray:
        """Coefficients of the targets on the independent symbols (zeta_1..zeta_K, eta_1..eta_N)."""
        return np.hstack([self.w * np.asarray(Ct @ self.Z.toarray()), (1.0 - self.w) * Ct])

    def sample(self, rng, n: int) -> np.ndarray:
        from .uncertainty import sample_eps
        zeta = sample_eps(rng, self.k, n)
        eta = sample_eps(rng, self.zone.size, n)
        return self.w * zeta[self.zone] + (1.0 - self.w) * eta

    @staticmethod
    def w_for_rho(rho: float) -> float:
        """Factor weight giving correlation rho between two sources of the same zone."""
        if rho <= 0:
            return 0.0
        a = np.sqrt(rho / (1.0 - rho)) if rho < 1 else np.inf
        return float(a / (1.0 + a)) if np.isfinite(a) else 1.0


# ----------------------------------------------------------------------------- sampled range of the quadratic part
def quadratic_ranges(so: SecondOrder, Ws: dict, n: int, seed: int, batch: int = 500, sampler=None):
    """Sampled (min, max) of the quadratic part of every target, for several target sets sharing the samples.
    Double precision is needed: q is a sum of large terms that nearly cancel (single precision errs by
    about 1e-7 p.u. on PEGASE-2869, the size of the ranges being estimated at many buses).
    ``sampler(rng, n)`` -> eps (N_u, n): the joint law of the symbols (default: independent U[-1, 1])."""
    from .uncertainty import sample_eps
    rng = np.random.default_rng(seed)
    WT = {k: sp.csr_matrix(W.T) for k, W in Ws.items()}
    out = {k: (np.full(W.shape[1], np.inf), np.full(W.shape[1], -np.inf)) for k, W in Ws.items()}
    for s in range(0, n, batch):
        b = min(batch, n - s)
        q = so.quadratic(sampler(rng, b) if sampler is not None else sample_eps(rng, so.c_lin.shape[1], b))
        for k, Wt in WT.items():
            v = Wt @ q
            np.minimum(out[k][0], v.min(axis=1), out=out[k][0])
            np.maximum(out[k][1], v.max(axis=1), out=out[k][1])
    return out


# ----------------------------------------------------------------------------- probability guarantee
def certified_split(alpha: float, n_samples: int) -> float:
    """Level left for the linear part when the range of the quadratic term comes from ``n_samples`` samples.

    For a target Y = L + Q of the second-order model, with the bound [-r + min(0, q_min), r + max(0, q_max)] and
    q_min, q_max the extremes of n independent samples of Q drawn from the law of the symbols,
        P(Y outside) <= P(|L| > r) + P(Q > q_max) + P(Q < q_min) <= alpha_L + 2 / (n + 1),
    because a new exchangeable sample exceeds the maximum (or falls below the minimum) of n others with
    probability at most 1 / (n + 1), whatever its distribution (order statistics; the probability is over both
    the new draw and the n samples). So alpha_L = alpha - 2 / (n + 1) gives P(Y outside) <= alpha for the
    second-order model; the third and higher orders are not covered by this argument.
    """
    a = alpha - 2.0 / (n_samples + 1)
    if a <= 0:
        raise ValueError(f"{n_samples} samples cannot certify level {alpha}: need more than {2 / alpha - 1:.0f}")
    return a


def certified_bernstein(Ct: np.ndarray, alpha: float, n_samples: int, var_eps: float = VAR_UNIFORM) -> np.ndarray:
    """Bernstein radius at the level that, with the sampled range of the quadratic term, certifies ``alpha``."""
    return bernstein_radius(Ct, certified_split(alpha, n_samples), var_eps)
