"""Parametric Krawczyk operator for the power-flow equations: the interval baseline of the voltage enclosure.

Network with PQ buses and reference bus(es). Unknown x = (u, w), real and imaginary parts of the voltages
U = u + j w of the PQ buses; parameters eps in [-1, 1]^K, on which the specified injections depend affinely,
S(eps) = S* + M eps (the same symbols as ``SecondOrder``). x* is the AC solution at eps = 0.

Operator. For a square system F(x, eps) = 0, a non-singular real matrix C and a box X = x* +- r, define
g(x) = x - C F(x, eps). If for every x in X and eps in the box  F(x, eps) - F(x*, eps) = S(x, eps) (x - x*)
with a matrix S(x, eps) (the slope; by the mean value theorem, row by row, S can be taken as the Jacobian F'_x
at points of X), then

    g(x) = x* - C F(x*, eps) + (I - C S(x, eps)) (x - x*)   in   K(X) = x* - C F(x*, P) + (I - C S(X, P)) (X - x*).

If K(X) lies in the interior of X, g maps X into itself for every eps, so by Brouwer's theorem it has a fixed
point there, which solves F(x, eps) = 0 because C is non-singular. With S = F'_x over all of X the solution is
unique in X (Krawczyk 1969; Moore 1977); with slopes only existence follows. Every solution of X lies in K(X),
so replacing X by its intersection with K(X) keeps them. The branch that continues x* cannot leave X either: wherever it is in X it is in
K(X), i.e. in the interior of X.

The pieces, with C = F'_x(x*, 0)^-1 (the inverse of the midpoint Jacobian, as in the published methods; it is
non-singular because the row sums of |I - C J*| are below one, which is checked):

    first order   F is affine in eps, so -C (F(x*, eps) - F(x*, 0)) = L eps with L = -C dF/deps. Its range is
                  taken exactly, component by component (row sums of |L| on the box; any other parameter set
                  through ``first_order``), so no dependency between the symbols is lost.
    centre        b0 = |C F(x*, 0)|: zero up to the tolerance of the load flow.
    non-linear    Phi(r) >= |(I - C S(x, eps)) (x - x*)| for x in x* +- r and eps in the box (below).

Then K(X) = x* + [lo1, hi1] +- (b0 + Phi(r)) and K(X) is inside X when  b + Phi(r) < r,  b = b0 + max(-lo1, hi1).
``solve`` iterates r <- b + Phi(r) from r = b (Phi is monotone, so this reaches the least fixed point whenever
some box x* +- r is mapped into itself), inflates the limit, checks the strict inequality, and then replaces
X by its intersection with K(X) until the box stops shrinking. The result also encloses the solution as an
affine form: |x(eps) - x* - L eps| <= b0 + Phi(r).

Formulation "power" (complex power mismatch, quadratic in x):
    F = U conj(Y U + c) - S(eps),   c: currents fed by the reference buses.
    F'_x(x) = J* + sum_j (x_j - x*_j) L_j is affine in x (L_j: Jacobian of the unit deviation j), so
    |C (F'_x(x) - J*)| <= sum_j r_j |C L_j| entrywise, exactly, with the cancellations inside C L_j kept.
    Slope: F is quadratic, so S(x) = F'_x((x + x*)/2) and the same bound holds with r_j / 2.
    Phi(r) = |I - C J*| r + theta sum_j r_j |C L_j| r,   theta = 1 (derivative) or 1/2 (slope).

Formulation "current" (current mismatch; the non-linearity is local to every bus):
    F = Y U + c - conj(S(eps)) / conj(U).
    d[conj(S)/conj(U)] = kap conj(dU),  kap = -conj(S)/conj(U)^2, so in real form F'_x = Yr - H(kap) with
    H(kap) = [[diag kap_r, diag kap_i], [diag kap_i, -diag kap_r]] and I - C F'_x = (I - C J*) + C H(kap - kap*).
    Row a of C H(dkap) dx is the sum over the buses of  Re( conj(c_am) dkap_m conj(dU_m) ),  c_am = C_a,um + j C_a,wm.
    Per bus, with sig = conj(S*)/conj(U*)^2, tau(eps) = conj(M eps)/conj(U*)^2, dxi = xi - U* for a point xi of
    the box, a >= |dxi|/|U*|, a < 1, and t >= |tau|:
        derivative   kap(xi)   - kap* = -tau + g conj(dxi) + rem,  g = 2 sig/conj(U*),
                                        |rem| <= 2 t a + (|sig| + t)(3 a^2 + 2 a^3)/(1 - a)^2
        slope        kaps(U)   - kap* = -tau + g conj(dU)  + rem,  g =   sig/conj(U*),   |rem| <= (t a + |sig| a^2)/(1 - a)
    (slope: conj(S)/conj(U) - conj(S)/conj(U*) = kaps conj(U - U*), kaps = -conj(S)/(conj(U) conj(U*))).
    The three terms of the contribution of a bus are bounded one by one, each as sharply as its form allows:
        -tau conj(dU)           bilinear in the symbols of the bus and dU. At a vertex of the symbols its maximum
                                over the rectangle of dU is |Re w| r_u + |Im w| r_w, w = conj(c) tau(eps); that is
                                convex in eps, so the maximum over the box is at a vertex and all are enumerated;
        g conj(dxi) conj(dU)    |w0| rho^2 with w0 = conj(c) g, rho = hypot(r_u, r_w), or the exact maximum over
                                two independent rectangles, whichever is smaller; with slopes dxi = dU and
                                |Re w0| max(r_u, r_w)^2 + 2 |Im w0| r_u r_w also holds;
        rem conj(dU)            |c| rem rho.
    t is the exact maximum of |tau(eps)| over the box. Buses without injection have sig = tau = 0 and do not
    contribute. The admittances do not enter the non-linear term.

``rotate``: every bus is written in axes turned by the angle of its voltage at x* (U = exp(-j theta*) V, Y
replaced by D^H Y D, D = diag exp(j theta*)). The equations keep their form, |V| = |U|, and to first order |V|
is the first coordinate, so the box loses little when it is turned into bounds of |V|. Without it the box is
aligned with the global axes and the rectangle of a bus with a large angle over-estimates |V|. (Rodrigues
Junior, Melo and Nepomuceno, IJEPES 2022, rotate all phasors together and intersect, for the same reason.)

|V| is bounded by the exact range of |u + j w| over the rectangle of each bus.

Valid in exact arithmetic. Nothing is rounded outward, so this is not a verified enclosure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .enclosure import injection_changes
from .second_order import SecondOrder

INFLATIONS = (1e-9, 1e-7, 1e-5, 1e-3, 1e-2, 1e-1)     # relative inflations tried for the strict inclusion
MAX_ENUMERATED_SYMBOLS = 12                           # symbols per bus up to which the vertices are enumerated


def modulus_range(u_lo, u_hi, w_lo, w_hi):
    """Exact (min, max) of |u + j w| over the rectangles [u_lo, u_hi] x [w_lo, w_hi]."""
    def near(lo, hi):
        return np.where((lo <= 0.0) & (hi >= 0.0), 0.0, np.minimum(np.abs(lo), np.abs(hi)))

    def far(lo, hi):
        return np.maximum(np.abs(lo), np.abs(hi))

    return np.hypot(near(u_lo, u_hi), near(w_lo, w_hi)), np.hypot(far(u_lo, u_hi), far(w_lo, w_hi))


def _hblock(kap):
    """H(kap): real form of dU -> kap conj(dU); rows (Re, Im), columns (u, w)."""
    r, i = np.diag(kap.real), np.diag(kap.imag)
    return np.block([[r, i], [i, -r]])


def _sign_patterns(k):
    """The 2^(k-1) vertices of [-1, 1]^k with first coordinate +1 (one of each pair v, -v): (patterns, k)."""
    idx = np.arange(2 ** (k - 1))
    return np.hstack([np.ones((len(idx), 1)), 2.0 * ((idx[:, None] >> np.arange(k - 1)[None, :]) & 1) - 1.0])


@dataclass
class KrawczykResult:
    closes: bool                # K(X) in the interior of X was shown: lo, hi, vm_lo, vm_hi hold over the parameter set
    iterations: int             # evaluations of Phi
    radius: np.ndarray          # (2m,) half-widths r of the final box X = x* +- r
    remainder: np.ndarray       # (2m,) |x(eps) - x* - L eps| <= remainder
    lo: np.ndarray              # (2m,) enclosure of (u, w) in the coordinates of the variant (NaN if it does not close)
    hi: np.ndarray
    vm_lo: np.ndarray           # (n_pq,) enclosure of |V| (NaN if it does not close)
    vm_hi: np.ndarray
    rel_radius: float           # largest |U - U*| / |U*| allowed by the box
    inflation: float            # relative inflation with which the strict inclusion held (NaN if none)


class Krawczyk:
    """Krawczyk operator at the centre of ``so`` (study and half-width); see the module docstring.

    ``formulation``: "power" or "current". ``slope``: slopes instead of derivatives (existence, not uniqueness).
    ``rotate``: bus-wise rotated coordinates. State order: u of the PQ buses (order of ``grid.pq``), then w.
    Only ``so.st`` and ``so.delta`` are used.
    """

    def __init__(self, so: SecondOrder, formulation: str = "current", slope: bool = False, rotate: bool = False):
        g, st = so.st.grid, so.st
        if len(g.pv):
            raise ValueError("the Krawczyk baseline needs all non-slack buses to be PQ")
        if formulation not in ("power", "current"):
            raise ValueError("formulation must be 'power' or 'current'")
        self.formulation, self.slope, self.rotate = formulation, bool(slope), bool(rotate)
        ns = g.pvpq
        fixed = np.setdiff1d(np.arange(g.n_bus), ns)
        m = self.m = len(ns)
        V = st.V_star
        rot = np.exp(1j * np.angle(V)) if rotate else np.ones(g.n_bus, complex)
        Y = (sp.diags(rot.conj()) @ g.Ybus.tocsr() @ sp.diags(rot)).tocsr()
        U = rot.conj() * V
        self.rot = rot[ns]
        self.Ynn = Y[ns][:, ns].tocsr()
        self.c = Y[ns][:, fixed] @ U[fixed]
        self.u_star = U[ns]
        self.x_star = np.r_[U[ns].real, U[ns].imag]
        self.s_star = np.asarray(st.sbus[ns], dtype=complex)
        self.M = injection_changes(so)[ns]                                   # (m, K) complex
        self.K = self.M.shape[1]
        Yd = self.Ynn.toarray()
        self.Yr = np.block([[Yd.real, -Yd.imag], [Yd.imag, Yd.real]])
        J = self.jacobian(self.x_star, np.zeros(self.K))
        self.C = np.linalg.inv(J)
        self.E0 = np.abs(np.eye(2 * m) - self.C @ J)
        self.regular = bool(self.E0.sum(axis=1).max() < 1.0)                 # then C (and J*) are non-singular
        self.b0 = np.abs(self.C @ self.residual(self.x_star, np.zeros(self.K)))
        B = self.M if formulation == "power" else np.conj(self.M) / np.conj(self.u_star)[:, None]
        self.L = self.C @ np.vstack([B.real, B.imag])                        # x(eps) - x* = L eps + higher orders
        self.theta = 0.5 if slope else 1.0
        if formulation == "power":
            self._unit = [self._linear_jacobian(j) for j in range(2 * m)]
            self._quad = []
            for A in self._unit:                                             # |C L_j|, stored by its non-zero columns
                cols = np.flatnonzero(np.diff(A.indptr))
                self._quad.append((cols, np.abs(self.C @ A[:, cols].toarray())))
            return
        self.ua = np.abs(self.u_star)
        self.sig = np.conj(self.s_star) / np.conj(self.u_star) ** 2
        tau = np.conj(self.M) / (np.conj(self.u_star) ** 2)[:, None]         # (m, K): tau(eps) = tau eps
        self.t_r, self.t_i = np.abs(tau.real).sum(axis=1), np.abs(tau.imag).sum(axis=1)
        self.g = (1.0 if slope else 2.0) * self.sig / np.conj(self.u_star)
        Cc = self.C[:, :m] + 1j * self.C[:, m:]                              # c_am
        self.c_abs = np.abs(Cc)
        self.w0 = np.conj(Cc) * self.g[None, :]
        self.t = np.zeros(m)
        self._bilinear = []                                                  # per bus with symbols: (bus, |Re w|, |Im w|) at the vertices
        for i in np.flatnonzero(np.abs(tau).sum(axis=1) > 0.0):
            ks = np.flatnonzero(tau[i])
            if len(ks) > MAX_ENUMERATED_SYMBOLS:                             # too many symbols on one bus: moduli are added
                self.t[i] = np.abs(tau[i]).sum()
                self._bilinear.append((i, None, None))
                continue
            at_vertex = _sign_patterns(len(ks)) @ tau[i, ks]                 # tau(eps) at the vertices, up to sign
            self.t[i] = np.abs(at_vertex).max()
            w = np.conj(Cc[:, i])[:, None] * at_vertex[None, :]
            self._bilinear.append((i, np.abs(w.real), np.abs(w.imag)))

    # ------------------------------------------------------------------ the equations, point by point
    def _split(self, x):
        return x[:self.m] + 1j * x[self.m:]

    def residual(self, x, eps):
        """F(x, eps) as (Re, Im)."""
        U = self._split(np.asarray(x, dtype=float))
        S = self.s_star + self.M @ np.asarray(eps, dtype=float)
        F = U * np.conj(self.Ynn @ U + self.c) - S if self.formulation == "power" else self.Ynn @ U + self.c - np.conj(S) / np.conj(U)
        return np.r_[F.real, F.imag]

    def _power_jacobian(self, U, cur):
        """Real form of dS = dU conj(cur) + U conj(Y dU): linear in (U, cur)."""
        Yc = self.Ynn.conj()
        Ae = sp.diags(U) @ Yc + sp.diags(np.conj(cur))
        Af = 1j * (sp.diags(np.conj(cur)) - sp.diags(U) @ Yc)
        return sp.bmat([[Ae.real, Af.real], [Ae.imag, Af.imag]]).tocsc()

    def _linear_jacobian(self, j):
        """L_j = F'_x(x* + e_j) - F'_x(x*) of the power formulation."""
        a = np.zeros(self.m, complex)
        a[j % self.m] = 1.0 if j < self.m else 1j
        return self._power_jacobian(a, self.Ynn @ a)

    def jacobian(self, x, eps):
        """F'_x(x, eps), dense."""
        U = self._split(np.asarray(x, dtype=float))
        if self.formulation == "power":
            return self._power_jacobian(U, self.Ynn @ U + self.c).toarray()
        S = self.s_star + self.M @ np.asarray(eps, dtype=float)
        return self.Yr - _hblock(-np.conj(S) / np.conj(U) ** 2)

    def slope_matrix(self, x, eps):
        """S(x, eps) with F(x, eps) - F(x*, eps) = S (x - x*), dense."""
        x = np.asarray(x, dtype=float)
        if self.formulation == "power":
            return self.jacobian(0.5 * (x + self.x_star), eps)
        S = self.s_star + self.M @ np.asarray(eps, dtype=float)
        return self.Yr - _hblock(-np.conj(S) / (np.conj(self._split(x)) * np.conj(self.u_star)))

    # ------------------------------------------------------------------ bounds over X = x* +- r and the box of eps
    def _remainder(self, a):
        """Bound of |rem| in kap - kap* (derivative) or kaps - kap* (slope), per bus, for |dxi|/|U*| <= a < 1."""
        sa, t = np.abs(self.sig), self.t
        with np.errstate(divide="ignore", invalid="ignore"):
            if self.slope:
                return (t * a + sa * a * a) / (1.0 - a)
            return 2.0 * t * a + (sa + t) * (3.0 * a * a + 2.0 * a ** 3) / (1.0 - a) ** 2

    def nonlinear_bound(self, r):
        """(Phi(r), largest |U - U*|/|U*|); Phi is None where the bound breaks down (some a >= 1)."""
        r = np.asarray(r, dtype=float)
        m = self.m
        ru, rw = r[:m], r[m:]
        rho = np.hypot(ru, rw)
        a = rho / np.abs(self.u_star)
        a_max = float(a.max())
        out = self.E0 @ r
        if self.formulation == "power":
            for j, (cols, P) in enumerate(self._quad):
                if r[j] > 0.0:
                    out = out + self.theta * r[j] * (P @ r[cols])
            return out, a_max
        if a_max >= 1.0:
            return None, a_max
        re, im = self.w0.real, self.w0.imag
        quad = np.minimum(np.abs(self.w0) * rho ** 2,                        # g conj(dxi) conj(dU): disc, or two independent rectangles
                          np.maximum(np.abs(re * ru + im * rw) * ru + np.abs(im * ru - re * rw) * rw,
                                     np.abs(re * ru - im * rw) * ru + np.abs(im * ru + re * rw) * rw))
        if self.slope:                                                       # dxi = dU: a0 (du^2 - dw^2) + 2 b0 du dw
            quad = np.minimum(quad, np.abs(re) * np.maximum(ru, rw) ** 2 + 2.0 * np.abs(im) * ru * rw)
        out = out + (quad + self.c_abs * (self._remainder(a) * rho)).sum(axis=1)
        for i, w_re, w_im in self._bilinear:                                 # -tau conj(dU)
            if w_re is None:
                out = out + self.c_abs[:, i] * self.t[i] * rho[i]
            else:
                out = out + (w_re * ru[i] + w_im * rw[i]).max(axis=1)
        return out, a_max

    def matrix_interval(self, r):
        """(lo, hi): entrywise enclosure of F'_x(x, eps) (of the slope matrix if ``slope``) over X = x* +- r and the box."""
        J = self.jacobian(self.x_star, np.zeros(self.K))
        r = np.asarray(r, dtype=float)
        if self.formulation == "power":
            rad = np.zeros_like(J)
            for j, A in enumerate(self._unit):
                rad += self.theta * r[j] * np.abs(A.toarray())
            return J - rad, J + rad
        ru, rw = r[:self.m], r[self.m:]
        a = np.hypot(ru, rw) / self.ua
        if a.max() >= 1.0:
            raise ValueError("the box reaches U = 0 at some bus")
        rem = self._remainder(a)
        gr, gi = np.abs(self.g.real), np.abs(self.g.imag)
        kr = self.t_r + gr * ru + gi * rw + rem                              # |Re(kap - kap*)|, |Im(kap - kap*)|
        ki = self.t_i + gi * ru + gr * rw + rem
        rad = np.block([[np.diag(kr), np.diag(ki)], [np.diag(ki), np.diag(kr)]])
        return J - rad, J + rad

    # ------------------------------------------------------------------ the iteration
    def solve(self, first_order=None, max_iter: int = 5000, tol: float = 1e-13) -> KrawczykResult:
        """Krawczyk enclosure over the box of eps, or over a subset of it whose range of L eps is ``first_order``.

        ``first_order``: (lo, hi), per state component, of L eps over the parameter set (e.g.
        ``cdhim.halfdisk.first_order_range(self.L)`` for half disks inside their bounding boxes).
        """
        m = self.m
        if first_order is None:
            hi1 = np.abs(self.L).sum(axis=1)
            lo1 = -hi1
        else:
            lo1, hi1 = (np.asarray(v, dtype=float) for v in first_order)
        b = self.b0 + np.maximum(-lo1, hi1)
        nan = np.full(2 * m, np.nan)
        r, n_eval, a = b.copy(), 0, 0.0
        for _ in range(max_iter):                                # least fixed point of r = b + Phi(r), from below
            phi, a = self.nonlinear_bound(r)
            n_eval += 1
            if not self.regular or phi is None or not np.all(np.isfinite(phi)) or phi.max() > 1e3:
                return KrawczykResult(False, n_eval, r, nan, nan, nan, nan[:m], nan[:m], a, float("nan"))
            new = b + phi
            done = np.abs(new - r).max() <= tol * max(new.max(), 1e-300)
            r = new
            if done:
                break
        used = float("nan")
        for eta in INFLATIONS:                                   # epsilon-inflation: K(X) strictly inside X
            cand = r * (1.0 + eta) + 1e-15 * r.max()
            phi, a_c = self.nonlinear_bound(cand)
            n_eval += 1
            if phi is not None and np.all(np.isfinite(phi)) and np.all(b + phi < cand):
                r, a, used = cand, a_c, eta
                break
        if not np.isfinite(used):
            return KrawczykResult(False, n_eval, r, nan, nan, nan, nan[:m], nan[:m], a, used)
        for _ in range(max_iter):                                # X <- intersection of K(X) and X until it stops shrinking
            phi, a_c = self.nonlinear_bound(r)
            n_eval += 1
            new = np.minimum(r, b + phi)
            shrink = (r - new).max()
            r, a = new, a_c
            if shrink <= 1e-12 * r.max():
                break
        phi, a = self.nonlinear_bound(r)
        rem = self.b0 + phi
        lo, hi = self.x_star + lo1 - rem, self.x_star + hi1 + rem
        vm_lo, vm_hi = modulus_range(lo[:m], hi[:m], lo[m:], hi[m:])
        return KrawczykResult(True, n_eval + 1, r, rem, lo, hi, vm_lo, vm_hi, a, used)
