"""Second-order model of the bus-voltage magnitudes under bounded injection uncertainty.

In rectangular coordinates x = (e, f) the power-flow equations are exactly quadratic,
    g(x* + D) = g(x*) + J D + q(D),     q_m(D) = dV_m conj((Y dV)_m)   (|dV_m|^2 at PV buses),
so the deviation of the state solves  D = C eps - J^-1 q(D)  and
    d|V_i| = c_i' eps + eps' H_i eps + O(delta^3),
    eps' H_i eps = -Re sum_m gam_im dV_m conj(dI_m) - sum_pv bv_im |dV_m|^2 + (c_perp_i' eps)^2 / (2 |V_i|),
with dV = E eps, dI = Y dV, and (gam_i, bv_i) = J^-T (d|V_i|/dx): one transposed solve per monitored bus.
H_i is the Hessian of the AC solution map (tests/test_second_order.py checks it by finite differences);
the same quantity appears, for the state, as d2x/dxi2 in Liao et al. (2017), obtained there with one
linear solve per pair of perturbations.

Also here: enclosures of eps' H eps on [-1, 1]^n used as references (term by term, spectral blocks).
The group's RAA-K operator is in cdhim.knomial.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from numba import njit
from scipy.sparse.linalg import splu

from .powerflow import solve_blocked
from .study import drive_matrix


def l1(x):
    return np.abs(x).sum(axis=1)


def pos(x):
    return np.maximum(x, 0.0)


def neg(x):
    return np.maximum(-x, 0.0)


# ----------------------------------------------------------------------------- second-order model
class SecondOrder:
    """First- and second-order terms of d|V| at every PQ bus, kept as sums of products of affine forms."""

    def __init__(self, st, delta):
        g = st.grid
        Y, V = g.Ybus.tocsr(), st.V_star
        n, ns, pq = g.n_bus, g.pvpq, g.pq
        m, npq = len(ns), len(pq)
        col = np.full(n, -1)
        col[ns] = np.arange(m)
        cur = Y @ V
        dS_de = (sp.diags(V) @ Y.conj() + sp.diags(cur.conj())).tocsr()
        dS_df = (1j * (sp.diags(cur.conj()) - sp.diags(V) @ Y.conj())).tocsr()
        blocks = [sp.hstack([dS_de[ns][:, ns].real, dS_df[ns][:, ns].real]),         # P at pv + pq
                  sp.hstack([dS_de[pq][:, ns].imag, dS_df[pq][:, ns].imag])]         # Q at pq
        if len(g.pv):                                                                # |V|^2 at pv
            k = np.arange(len(g.pv))
            blocks.append(sp.hstack([sp.csr_matrix((2 * V[g.pv].real, (k, col[g.pv])), shape=(len(k), m)),
                                     sp.csr_matrix((2 * V[g.pv].imag, (k, col[g.pv])), shape=(len(k), m))]))
        lu = splu(sp.vstack(blocks).tocsc())

        src = st.src
        row_p, row_q = np.full(n, -1), np.full(n, -1)
        row_p[ns] = np.arange(m)
        row_q[pq] = m + np.arange(npq)
        rhs = (drive_matrix(src, 2 * m, row_p, np.where(src.is_q, row_q[src.bus], row_p[src.bus])) * delta).toarray()
        C = solve_blocked(lu, rhs)                                   # (de, df) of the non-slack buses
        E = np.zeros((n, src.n), complex)                            # dV = E eps (zero at the slack)
        E[ns] = C[:m] + 1j * C[m:]
        F = Y @ E                                                    # d(injected current) = F eps

        vm = np.abs(V[pq])
        ce, cf = C[col[pq]], C[m + col[pq]]
        er, ei = (V[pq].real / vm)[:, None], (V[pq].imag / vm)[:, None]
        self.c_lin = er * ce + ei * cf                               # d|V| first order
        self.c_perp = er * cf - ei * ce                              # component normal to V
        alpha = np.zeros((2 * m, npq))
        alpha[col[pq], np.arange(npq)] = er[:, 0]
        alpha[m + col[pq], np.arange(npq)] = ei[:, 0]
        beta = np.empty_like(alpha)                                  # J^-T alpha, one column per PQ bus
        for s in range(0, npq, 32):
            beta[:, s:s + 32] = lu.solve(np.ascontiguousarray(alpha[:, s:s + 32]), trans="T")
        gam = np.zeros((n, npq), complex)                            # weight of q_m in d|V_i|
        gam[ns] = beta[:m]
        gam[pq] -= 1j * beta[m:m + npq]
        bv = np.zeros((n, npq))
        bv[g.pv] = beta[m + npq:]
        self.st, self.Y, self.E, self.F, self.gam, self.bv, self.vm = st, Y, E, F, gam, bv, vm
        self.lu, self.C, self.col, self.m, self.delta = lu, C, col, m, delta
        self.has_pv = bool(len(g.pv))

    def quadratic(self, eps):
        """eps' H_i eps for every PQ bus and every column of ``eps``, without forming any H_i."""
        a = self.E @ eps
        z = a * np.conj(self.Y @ a)
        out = -(self.gam.T @ z).real + (self.c_perp @ eps) ** 2 / (2.0 * self.vm[:, None])
        if self.has_pv:
            out -= self.bv.T @ (a.real ** 2 + a.imag ** 2)
        return out

    def trace(self):
        s = np.einsum("mk,mk->m", self.E, self.F.conj())
        e2 = np.einsum("mk,mk->m", self.E, self.E.conj()).real
        return (-(self.gam.real.T @ s.real - self.gam.imag.T @ s.imag) - self.bv.T @ e2
                + (self.c_perp ** 2).sum(axis=1) / (2.0 * self.vm))

    def term_bounds(self):
        """Enclosures of eps' H_i eps for all buses from the terms  -Re(gam_m dV_m conj(dI_m)).

        natural:  each real product (a'eps)(b'eps) in [-|a|_1 |b|_1, |a|_1 |b|_1]  (standard AA product);
        binomial: (a'eps)(b'eps) = ((t a + b/t)'eps)^2/4 - ((t a - b/t)'eps)^2/4 with t^2 = |b|_1/|a|_1,
                  i.e. two k-nomial blocks per product; never wider than the natural bound.
        """
        Er, Ei, Fr, Fi = self.E.real, self.E.imag, self.F.real, self.F.imag
        Ar, Ai, Br, Bi = l1(Er), l1(Ei), l1(Fr), l1(Fi)

        def binomial(a, b, A, B):
            ok = (A > 0) & (B > 0)
            t = np.ones_like(A)
            t[ok] = np.sqrt(B[ok] / A[ok])
            hi = 0.25 * l1(t[:, None] * a + b / t[:, None]) ** 2
            lo = -0.25 * l1(t[:, None] * a - b / t[:, None]) ** 2
            hi[~ok] = lo[~ok] = 0.0
            return lo, hi

        gr, gi = self.gam.real, self.gam.imag
        sq = Ar ** 2 + Ai ** 2                                        # |dV_m|^2 <= sq_m (PV term)
        pv_lo, pv_hi = -pos(self.bv).T @ sq, neg(self.bv).T @ sq
        curv = l1(self.c_perp) ** 2 / (2.0 * self.vm)                 # (c_perp'eps)^2/(2|V|) in [0, curv]
        nat = np.abs(gr).T @ (Ar * Br + Ai * Bi) + np.abs(gi).T @ (Ai * Br + Ar * Bi)
        # products and their coefficients: -gr (ar br), -gr (ai bi), +gi (ai br), -gi (ar bi)
        l1_, h1 = binomial(Er, Fr, Ar, Br)
        l2_, h2 = binomial(Ei, Fi, Ai, Bi)
        l3_, h3 = binomial(Ei, Fr, Ai, Br)
        l4_, h4 = binomial(Er, Fi, Ar, Bi)
        hi = neg(gr).T @ (h1 + h2) - pos(gr).T @ (l1_ + l2_) + pos(gi).T @ h3 - neg(gi).T @ l3_ + neg(gi).T @ h4 - pos(gi).T @ l4_
        lo = neg(gr).T @ (l1_ + l2_) - pos(gr).T @ (h1 + h2) + pos(gi).T @ l3_ - neg(gi).T @ h3 + neg(gi).T @ l4_ - pos(gi).T @ h4
        self.natural_parts = {"products": 2.0 * nat, "pv_squares": pv_hi - pv_lo, "curvature": curv}
        return {"natural": (-nat + pv_lo, nat + pv_hi + curv), "binomial": (lo + pv_lo, hi + pv_hi + curv)}

    def forms(self):
        """H_i for every PQ bus, shape (n_pq, K, K). Memory grows as n_bus K^2: meant for few symbols."""
        E, F = self.E, self.F
        n, k = E.shape
        T = (E[:, :, None] * np.conj(F)[:, None, :]).reshape(n, k * k)
        H = -(self.gam.T @ T).real
        if self.has_pv:
            H -= self.bv.T @ (np.conj(E)[:, :, None] * E[:, None, :]).real.reshape(n, k * k)
        H = H.reshape(-1, k, k) + self.c_perp[:, :, None] * self.c_perp[:, None, :] / (2.0 * self.vm)[:, None, None]
        return 0.5 * (H + H.transpose(0, 2, 1))

    def assemble(self, k):
        """Dense symmetric H_i of PQ bus number ``k`` (N_u x N_u)."""
        M = (self.Y.conj().T @ sp.diags(self.gam[:, k]) + sp.diags(self.bv[:, k])).tocsr()
        Mh = 0.5 * (M + M.conj().T)
        H = -np.real(self.E.conj().T @ (Mh @ self.E))
        H += np.outer(self.c_perp[k], self.c_perp[k]) / (2.0 * self.vm[k])
        return 0.5 * (H + H.T)


# ----------------------------------------------------------------------------- bounds on eps' H eps
def bound_monomial(H):
    """k = 1 (AF2 / RAA term by term): sum H_kk [0, 1] + sum_{k<l} 2 H_kl [-1, 1]."""
    d = np.diag(H)
    off = np.abs(H).sum() - np.abs(d).sum()
    return float(neg(d).sum() * -1.0 - off), float(pos(d).sum() + off)


@njit(cache=True)
def _greedy_pass(Q, oi, oj):
    lo = hi = 0.0
    for k in range(len(oi)):
        i, j = oi[k], oj[k]
        di, dj = Q[i, i], Q[j, j]
        if di * dj <= 0.0:
            continue
        a = min(abs(Q[i, j]), abs(di), abs(dj))
        if a <= 0.0:
            continue
        r = 1.0 if di > 0.0 else -1.0
        if r > 0.0:
            hi += 4.0 * a                      # a (eps_i +- eps_j)^2 in a [0, 4]
        else:
            lo -= 4.0 * a
        s = a if Q[i, j] > 0.0 else -a
        Q[i, j] -= s
        Q[j, i] -= s
        Q[i, i] -= r * a
        Q[j, j] -= r * a
    return lo, hi


def bound_greedy_binomial(H):
    """k = 2: move a|eps_i eps_j| cross terms into a(eps_i +- eps_j)^2 while the diagonals can pay for it.

    One pass over the pairs in decreasing |H_ij|; a simplified reading of the diagonal-budgeting idea,
    not the authors' implementation. The residual is bounded monomial-wise.
    """
    Q = H.copy()
    iu, ju = np.triu_indices(len(Q), 1)
    order = np.argsort(-np.abs(Q[iu, ju]), kind="stable")
    lo, hi = _greedy_pass(Q, iu[order], ju[order])
    l2, h2 = bound_monomial(Q)
    return lo + l2, hi + h2


def bound_spectral(lam, vec):
    """k = n, no residual: H = sum lam_k v_k v_k'  ->  sum lam_k [0, |v_k|_1^2]."""
    s2 = np.abs(vec).sum(axis=0) ** 2
    return float(-(neg(lam) @ s2)), float(pos(lam) @ s2)


def bound_knomial(H, lam, vec, r):
    """k-nomial theorem with the r leading spectral blocks: sum sigma_l (w_l'eps)^2 + eps' R eps, R monomial-wise."""
    top = np.argsort(-np.abs(lam))[:r]
    W = vec[:, top]
    lo, hi = bound_monomial(H - (W * lam[top]) @ W.T)
    s2 = np.abs(W).sum(axis=0) ** 2
    return float(lo - neg(lam[top]) @ s2), float(hi + pos(lam[top]) @ s2)
