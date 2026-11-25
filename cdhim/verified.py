"""Floating-point rigorous enclosure of the voltage magnitudes of a network with PQ buses (a posteriori).

Problem. Non-slack buses are PQ; the injections are S(eps) = S0 + M eps, eps in [-1, 1]^K. In current form,

    F(V, eps) = Y V + c - conj(S(eps)) / conj(V) = 0,        c = Y_(L,slack) V_slack.

Model. Any floating-point data V* (m,), E (m, K), B (m, K, K) define m(eps) = V* + E eps + eps' B eps; they are
meant to be the AC solution at eps = 0 and its first- and second-order coefficients, but nothing below assumes
that they are exact. Likewise C is any real 2m x 2m matrix, meant to be the inverse of the Jacobian of F.

Operator. With kap = conj(S0)/conj(V*)^2 and the real-linear map J0[s] = Y s + kap conj(s),

    F(m + s) = F(m) + J0[s] + w conj(s),     w = conj(S) / (conj(m) conj(m + s)) - kap      (per bus),

so Phi(s) = s - C F(m + s) = (I - C J0)[s] - C F(m) - C [w conj(s)]. For |s_j| <= z_j and every eps of the box,

    |Phi(s)_i| <= sum_j R_ij z_j + sum_j T_ij (phi_j + omega_j(z) z_j),

    R_ij, T_ij >= spectral norm of the 2x2 real block (i, j) of I - C J0 and of C,
    phi_j     >= sup |F_j(m(eps), eps)|     (defect of the model),
    omega_j   >= sup |w_j| = [tau_j + |S0_j| (|v|(2 dl + z) + dl (dl + z)) / |v|^2] / ((|v| - dl)(|v| - dl - z)),
                 v = V*_j, tau_j = sum_k |M_jk|, dl_j >= |m_j - v|.

If the right-hand side is <= z_i for every i, Phi maps the product of discs into itself and Brouwer's theorem
gives s with C F(m + s) = 0; R z < z makes C non-singular, so F(m + s) = 0: for every eps there is a power-flow
solution with |V_j - m_j(eps)| <= z_j. If moreover (R + T diag(lam)) z < z, lam being the Lipschitz bound of
s -> w conj(s), Phi is a contraction and that solution is the only one in the tube.

Defect. conj(m_j) F_j(m) = G_j = conj(m_j)(Y m + c)_j - conj(S_j) is a polynomial of degree 4 in eps:

    degree 0   g0 = conj(V*) I* - conj(S0),                         I* = Y V* + c
    degree 1   g1 = conj(V*) (Y E) + conj(E) I* - conj(M)
    degree 2   g2 = conj(V*) (Y B) + conj(B) I* + sym[conj(E)_k (Y E)_l]
    degree 3+  conj(a)(Y b) + conj(b)(Y a) + conj(b)(Y b),         a = E eps, b = eps' B eps.

g0, g1 and g2 vanish for exact data; here they are evaluated with error bounds and kept. The rest is bounded
with norms of the computed Y E and Y B (the changes of the injected currents), so large admittances do not enter.

Magnitudes. |V_i| is within z_i of |v + a + b|. With p, q the components of a + b along and across v and
A = |v| + p > 0,  A + q^2/(2A) - q^4/(8A^3) <= sqrt(A^2 + q^2) <= A + q^2/(2A), hence

    |v + a + b| = |v| + c_lin' eps + eps' H eps + rem,     H = Re(conj(v) B)/|v| + c_perp c_perp'/(2|v|),
    |rem| <= Q1^2 Pm/(2|v| den) + (2 Q1 Q2 + Q2^2)/(2 den) + (Q1 + Q2)^4/(8 den^3),     den = |v| - Pm,

with Q1, Q2 >= |q| at first and second order and Pm >= |p|.

Range of c' eps + eps' H eps over the set S (the box, or the half disks of cdhim.halfdisk). Two bounds, the
better of which is kept at every bus and on every side:
    separate   range of c' eps over S (exact) plus the range of the form over the box with the RAA-K blocks;
    joint      for any vectors u_b and numbers a_b,  -(u_b' eps)^2 <= a_b^2 - 2 a_b u_b' eps  (a tangent of a
               concave function lies above it), so with H = Hrest - sum_b u_b u_b',
                   c' eps + eps' H eps <= max_S (c - 2 sum_b a_b u_b)' eps + sum_b a_b^2 + max_box eps' Hrest eps.
               The u_b are the eigenvectors of the negative eigenvalues of H (scaled), the a_b their values at
               the maximiser of the linear part, refined a few times; Hrest is bounded term by term. The lower
               bound is the mirror image with the positive eigenvalues. Neither the eigenvectors nor the a_b
               need be exact: Hrest is formed from them with an error bound.
The separate bound loses the second-order distance between the maximisers of the two terms; the joint one does
not when the form is concave (upper bound) or convex (lower bound) along the directions that matter.

Rounding. Every bound is an upper bound in real arithmetic of the stated quantity for the stored floating-point
data: after each elementary operation the result is moved to the next floating-point number in the safe
direction; sums and matrix products are corrected with the a priori bound gamma_n = n u / (1 - n u) of a sum of
n terms (Higham, Accuracy and Stability of Numerical Algorithms, 2002, Sec. 3.1), which holds for any order of
evaluation; quantities with cancellation (the defects) are computed in floating point and their error is bounded
from the magnitudes of their terms. This assumes IEEE-754 double arithmetic rounded to nearest with correctly
rounded +, -, *, / and sqrt, and no underflow (the magnitudes involved are far from it). The network data
(Ybus, S0, M) are taken as exact: the statement is about the equations defined by the stored numbers.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .enclosure import injection_changes, state_forms
from .knomial import knomial_rewrite
from .second_order import SecondOrder

U = 2.0 ** -53                               # unit roundoff of IEEE double, round to nearest


# ----------------------------------------------------------------------------- directed bounds
def up(x):
    return np.nextafter(x, np.inf)


def down(x):
    return np.nextafter(x, -np.inf)


def add_up(*xs):
    """Upper bound of the sum of upper bounds."""
    out = xs[0]
    for x in xs[1:]:
        out = up(out + x)
    return out


def mul_up(a, b):
    """Upper bound of the product of non-negative upper bounds."""
    return up(a * b)


def div_up(a, b_lo):
    """Upper bound of a / b for a non-negative upper bound ``a`` and a positive lower bound ``b_lo``."""
    return up(a / b_lo)


def sub_down(a_lo, b_up):
    """Lower bound of a - b from a lower bound of a and an upper bound of b."""
    return down(a_lo - b_up)


def _nterms(x, axis):
    if axis is None:
        return max(x.size, 1)
    axis = (axis,) if isinstance(axis, int) else axis
    return max(int(np.prod([x.shape[a] for a in axis])), 1)


def sum_up(x, axis=None):
    """Upper bound of the exact sum of non-negative terms: fl(sum) >= sum (1 - gamma_(n-1))."""
    return up(np.sum(x, axis=axis) * (1.0 + 4.0 * _nterms(x, axis) * U))


def ssum_up(x, axis=None):
    """Upper bound of the exact sum of signed terms: |fl(sum) - sum| <= gamma_(n-1) sum |x|."""
    return up(np.sum(x, axis=axis) + up(4.0 * _nterms(x, axis) * U * sum_up(np.abs(x), axis)))


def ssum_down(x, axis=None):
    return down(np.sum(x, axis=axis) - up(4.0 * _nterms(x, axis) * U * sum_up(np.abs(x), axis)))


def abs_up(z):
    """Upper bound of |z|, entrywise."""
    if not np.iscomplexobj(z):
        return np.abs(z)
    return up(np.sqrt(up(up(z.real * z.real) + up(z.imag * z.imag))))


def abs_down(z):
    """Lower bound of |z|, entrywise."""
    if not np.iscomplexobj(z):
        return np.abs(z)
    a, b = np.maximum(down(z.real * z.real), 0.0), np.maximum(down(z.imag * z.imag), 0.0)
    return np.maximum(down(np.sqrt(np.maximum(down(a + b), 0.0))), 0.0)


def matmul_up(A, B, p):
    """Upper bound of A B for non-negative A, B; ``p`` >= number of terms of every inner product."""
    if sp.issparse(B) and not sp.issparse(A):
        P = (B.T.tocsr() @ A.T).T
    else:
        P = A @ B
    P = P.toarray() if sp.issparse(P) else np.asarray(P)
    return up(P * (1.0 + 4.0 * (p + 1) * U))


def err_const(terms):
    """c with |fl(e) - e| <= c * (sum of the magnitudes of the terms) for a sum of ``terms`` real or complex
    products evaluated in any order: sqrt(2) gamma_(terms + 2) <= 4 (terms + 2) u."""
    return 4.0 * (terms + 2) * U


def _sparse_abs(A):
    A = A.tocsr()
    return sp.csr_matrix((abs_up(A.data), A.indices, A.indptr), shape=A.shape)


def csum_up(A, axis, rot=None):
    """Upper bound of sup |sum_k A_k eps_k| over the box for complex coefficients: the smallest of the sum of the
    moduli and of hypot(sum |Re A_k|, sum |Im A_k|) in the given axes and, with ``rot`` (one number of modulus
    close to 1 per row), in the axes turned by it."""
    def hyp(B):
        r, i = sum_up(np.abs(B.real), axis=axis), sum_up(np.abs(B.imag), axis=axis)
        return up(np.sqrt(up(up(r * r) + up(i * i))))

    l1 = sum_up(abs_up(A), axis=axis)
    out = np.minimum(l1, hyp(A))
    if rot is not None:
        shape = (-1,) + (1,) * (A.ndim - 1)
        turned = up(add_up(hyp(A * rot.reshape(shape)), up(8.0 * U * l1)) / down(abs_down(rot) * (1.0 - 8.0 * U)))
        out = np.minimum(out, turned)
    return out


def _form_sup(A):
    """Upper bound of sup |eps' A eps| on the box for a stack of real forms (n, K, K): diagonal in [0, 1]."""
    k = A.shape[1]
    d = np.einsum("ikk->ik", A)
    off = sum_up(np.abs(A) * (1.0 - np.eye(k)), axis=(1, 2))
    return add_up(off, np.maximum(sum_up(np.maximum(d, 0.0), axis=1), sum_up(np.maximum(-d, 0.0), axis=1)))


def _cform_sup(A, rot=None):
    """The same for complex forms: the smallest of the l1 norm and of the term-by-term bound of the real and
    imaginary parts, in the given axes and, with ``rot``, in the axes turned by it (see ``csum_up``)."""
    def parts(B):
        r, i = _form_sup(B.real), _form_sup(B.imag)
        return up(np.sqrt(up(up(r * r) + up(i * i))))

    l1 = sum_up(abs_up(A), axis=(1, 2))
    out = np.minimum(l1, parts(A))
    if rot is not None:
        turned = up(add_up(parts(A * rot[:, None, None]), up(8.0 * U * l1)) / down(abs_down(rot) * (1.0 - 8.0 * U)))
        out = np.minimum(out, turned)
    return out


# ----------------------------------------------------------------------------- the tube around the model
@dataclass
class Tube:
    exists: bool                # self-mapping verified: a solution within z of the model for every eps
    unique: bool                # contraction verified: it is the only solution within z of the model
    z: np.ndarray               # radius per non-slack bus [p.u.]
    iterations: int
    reason: str                 # why the verification failed ("" if it did not)
    rounding_share: float       # largest share of z explained by the defects of degree 0 to 2 (rounding level)
    aux: dict | None = None     # with ``keep``: C, R, T, phi and kap, for the tests of the derivation


def verify_tube(Ybus, ns, V, S0, M, E, B, s_err=None, inflate=(1e-9, 1e-6, 1e-3), max_iter: int = 5000, keep: bool = False,
                cache: dict | None = None) -> Tube:
    """Radius z of the tube around m(eps) = V[ns] + E eps + eps' B eps, verified with directed bounds.

    ``Ybus`` (n, n) sparse; ``ns`` the non-slack buses, all PQ; ``V`` (n,) the centre, slack voltages included;
    ``S0`` (m,), ``M`` (m, K) the injections of the non-slack buses; ``E`` (m, K), ``B`` (m, K, K) the model.
    ``s_err`` (m,): the true injections are S0 + M eps + e(eps) with |e_j| <= s_err_j (e.g. the rounding of a
    centre computed in floating point); it adds to the defect and to the bound of the injection change.
    ``cache``: a dictionary shared by several calls on the same network. Both C and the splitting point kap can
    be any numbers, so those of the first call are reused in the next ones together with the bounds R and T,
    which depend only on them; what changes with the centre is the distance from kap to conj(S0)/conj(V*)^2,
    which enters omega. The tube is then somewhat larger than with matrices computed for that centre.
    """
    Y = Ybus.tocsr()
    m, K = E.shape
    Yn = Y[ns].tocsr()
    YL = Yn[:, ns].tocsr()
    p = max(int(np.diff(Yn.indptr).max()), 1)
    aYn, aYL = _sparse_abs(Yn), _sparse_abs(YL)
    Vn, Vc = V[ns], np.conj(V[ns])
    vup, vlo = abs_up(Vn), abs_down(Vn)
    nan = np.full(m, np.nan)

    # currents of the model and their error bounds
    cur = Yn @ V
    e_cur = up(err_const(p) * matmul_up(aYn, abs_up(V), p))
    a_cur = abs_up(cur)
    YE = YL @ E
    aE = abs_up(E)
    e_YE = up(err_const(p) * matmul_up(aYL, aE, p))
    aYE = abs_up(YE)
    YB = (YL @ B.reshape(m, K * K)).reshape(m, K, K)
    aB = abs_up(B)
    e_YB = up(err_const(p) * matmul_up(aYL, aB.reshape(m, K * K), p)).reshape(m, K, K)
    aYB = abs_up(YB)
    aS0, aM = abs_up(S0), abs_up(M)

    # defects of degree 0, 1, 2 (zero for exact data) with the error of their evaluation
    g0 = Vc * cur - np.conj(S0)
    G0 = add_up(abs_up(g0), mul_up(vup, e_cur), up(err_const(2) * add_up(mul_up(vup, a_cur), aS0)))
    g1 = Vc[:, None] * YE + np.conj(E) * cur[:, None] - np.conj(M)
    e1 = add_up(mul_up(vup[:, None], e_YE), mul_up(aE, e_cur[:, None]),
                up(err_const(3) * add_up(mul_up(vup[:, None], aYE), mul_up(aE, a_cur[:, None]), aM)))
    G1 = sum_up(add_up(abs_up(g1), e1), axis=1)
    cross = np.conj(E)[:, :, None] * YE[:, None, :]
    g2 = Vc[:, None, None] * YB + np.conj(B) * cur[:, None, None] + 0.5 * (cross + cross.transpose(0, 2, 1))
    a_cross, e_cross = mul_up(aE[:, :, None], aYE[:, None, :]), mul_up(aE[:, :, None], e_YE[:, None, :])
    e2 = add_up(mul_up(vup[:, None, None], e_YB), mul_up(aB, e_cur[:, None, None]), 0.5 * add_up(e_cross, e_cross.transpose(0, 2, 1)),
                up(err_const(4) * add_up(mul_up(vup[:, None, None], aYB), mul_up(aB, a_cur[:, None, None]),
                                         0.5 * add_up(a_cross, a_cross.transpose(0, 2, 1)))))
    G2 = sum_up(add_up(abs_up(g2), e2), axis=(1, 2))
    if s_err is not None:
        G0 = add_up(G0, s_err)
    G012 = add_up(G0, G1, G2)

    # degree 3 and 4
    rot = Vc / np.sqrt(Vn.real * Vn.real + Vn.imag * Vn.imag)          # axes along and across V*
    alpha = csum_up(E, 1, rot)                                         # >= |a|
    beta = _cform_sup(B, rot)                                          # >= |b|
    YEn = add_up(csum_up(YE, 1, rot), sum_up(e_YE, axis=1))            # >= |Y a|
    YBn = add_up(_cform_sup(YB, rot), sum_up(e_YB, axis=(1, 2)))       # >= |Y b|
    G3 = add_up(mul_up(alpha, YBn), mul_up(beta, YEn), mul_up(beta, YBn))
    delta = add_up(alpha, beta)                                        # >= |m - v|
    den = sub_down(vlo, delta)                                         # <= |m|
    if not np.all(den > 0):
        return Tube(False, False, nan, 0, "the model reaches zero voltage", np.nan)
    phi = div_up(add_up(G012, G3), den)
    phi_round = div_up(G012, den)

    # splitting point kap (any floating-point value) and the distance to conj(S0)/conj(V*)^2
    tau = csum_up(M, 1) if s_err is None else add_up(csum_up(M, 1), s_err)
    shared = cache is not None and "R" in cache
    kap = cache["kap"] if shared else np.conj(S0) / (Vc * Vc)
    vlo2 = down(vlo * vlo)
    e_kap = div_up(add_up(abs_up(np.conj(S0) - kap * Vc * Vc), up(err_const(4) * add_up(aS0, mul_up(abs_up(kap), mul_up(vup, vup))))), vlo2)

    if shared:
        C, R, T = cache["C"], cache["R"], cache["T"]
    else:
        # C, an approximate inverse of J0 = Yr + Kr in real form, and bounds of I - C J0
        Gs, Bs = YL.real, YL.imag
        Yr = sp.bmat([[Gs, -Bs], [Bs, Gs]], format="csc")
        kr, ki = kap.real, kap.imag
        Kr = sp.bmat([[sp.diags(kr), sp.diags(ki)], [sp.diags(ki), sp.diags(-kr)]], format="csc")
        C = np.linalg.inv((Yr + Kr).toarray())
        aC = np.abs(C)
        pc = max(int(np.diff(Yr.indptr).max()), 1)
        P1 = (Yr.T.tocsr() @ C.T).T
        e_P1 = up(err_const(pc) * matmul_up(aC, abs(Yr), pc))
        P2 = np.hstack([C[:, :m] * kr + C[:, m:] * ki, C[:, :m] * ki - C[:, m:] * kr])
        akr, aki = np.abs(kr), np.abs(ki)
        e_P2 = up(err_const(2) * np.hstack([add_up(mul_up(aC[:, :m], akr), mul_up(aC[:, m:], aki)),
                                            add_up(mul_up(aC[:, :m], aki), mul_up(aC[:, m:], akr))]))
        eye = np.eye(2 * m)
        res = add_up(np.abs(eye - P1 - P2), e_P1, e_P2, up(err_const(3) * add_up(eye, np.abs(P1), np.abs(P2))))

        def blocks(A):
            return A[:m, :m], A[:m, m:], A[m:, :m], A[m:, m:]

        a, b, c, d = blocks(res)
        R = up(np.sqrt(add_up(mul_up(a, a), mul_up(b, b), mul_up(c, c), mul_up(d, d))))    # Frobenius >= spectral
        a, b, c, d = blocks(C)
        s = add_up(up(a * a), up(b * b), up(c * c), up(d * d))
        e_det = up(err_const(2) * add_up(up(np.abs(a * d)), up(np.abs(b * c))))
        det_lo = np.maximum(down(np.abs(a * d - b * c) - e_det), 0.0)
        disc = np.maximum(up(up(s * s) - 4.0 * np.maximum(down(det_lo * det_lo), 0.0)), 0.0)
        T = up(np.sqrt(up(0.5 * up(s + up(np.sqrt(disc))))))                               # 2x2 spectral norms
        if cache is not None:
            cache.update(kap=kap, C=C, R=R, T=T)

    def omega(z, lipschitz=False):
        d2 = sub_down(den, z)
        if not np.all(d2 > 0):
            return None
        if lipschitz:
            dz = add_up(delta, z)
            num, dd = add_up(mul_up(vup, 2.0 * dz), mul_up(dz, dz)), down(d2 * d2)
        else:
            num, dd = add_up(mul_up(vup, add_up(2.0 * delta, z)), mul_up(delta, add_up(delta, z))), down(den * d2)
        return add_up(div_up(add_up(tau, div_up(mul_up(aS0, num), vlo2)), dd), e_kap)

    # the least fixed point of z -> R z + T (phi + omega(z) z), in plain floating point
    z = T @ phi
    it, reason = 0, "the iteration did not settle"
    for it in range(1, max_iter + 1):
        w = omega(z)
        if w is None:
            reason = "the tube reaches zero voltage"
            break
        new = R @ z + T @ (phi + w * z)
        if not np.all(np.isfinite(new)):
            reason = "the iteration diverged"
            break
        settled = np.abs(new - z).max() <= 1e-12 * max(new.max(), 1e-300)
        z = new
        if settled:
            reason = ""
            break
    if reason:
        return Tube(False, False, nan, it, reason, np.nan)

    # the verification proper: the inflated z is mapped into itself, with every operation rounded up
    for f in inflate:
        zz = up(z * (1.0 + f))
        w = omega(zz)
        if w is None:
            continue
        Rz = matmul_up(R, zz, m)
        if np.all(add_up(Rz, matmul_up(T, add_up(phi, mul_up(w, zz)), m)) <= zz) and np.all(Rz < zz):
            lam = omega(zz, True)
            unique = lam is not None and bool(np.all(add_up(Rz, matmul_up(T, mul_up(lam, zz), m)) < zz))
            share = float(np.max(matmul_up(T, phi_round, m) / zz))
            return Tube(True, unique, zz, it, "", share, dict(C=C, R=R, T=T, phi=phi, kap=kap) if keep else None)
    return Tube(False, False, nan, it, "the inflated fixed point is not mapped into itself", np.nan)


# ----------------------------------------------------------------------------- quadratic form and magnitudes
def raak_range(H, max_k: int, node_budget: int = 0):
    """Rigorous (lo, hi, blocks) of eps' H eps on [-1, 1]^K for the floating-point matrix ``H``.

    The blocks come from the (non-rigorous) search on H / max|H|; with s = max|H| and R = H - s sum_b c_b u_b u_b',
    eps' H eps = s sum_b c_b (u_b' eps)^2 + eps' R eps exactly, whatever the blocks. R is computed in floating
    point and its error is bounded entrywise.
    """
    k = len(H)
    s = float(np.abs(H).max()) if k else 0.0
    if s == 0.0:
        return 0.0, 0.0, 0
    c, Ub, _, _ = knomial_rewrite(H / s, None, max_k, node_budget=node_budget)
    nb = len(c)
    aU = np.abs(Ub)
    R = H - s * ((Ub.T * c) @ Ub)
    e_R = up(err_const(nb + 2) * add_up(np.abs(H), up(s * matmul_up(aU.T, aU, max(nb, 1)))))
    l1 = sum_up(aU, axis=1)
    blk = up(s * up(l1 * l1))                                        # s |u_b|_1^2
    d = np.diag(R)
    off = add_up(sum_up(np.abs(R) * (1.0 - np.eye(k))), sum_up(e_R))
    hi = add_up(sum_up(blk[c > 0]), sum_up(np.maximum(d, 0.0)), off)
    lo = -add_up(sum_up(blk[c < 0]), sum_up(np.maximum(-d, 0.0)), off)
    return float(lo), float(hi), nb


def _lin_range(C, half_disks):
    """Rigorous (lo, hi) of C_j' eps over the box, or over the product of half disks, for every row of C."""
    if not half_disks:
        r = sum_up(np.abs(C), axis=1)
        return -r, r
    n = C.shape[1] // 2
    c1, c2 = C[:, :n], C[:, n:]

    def hyp(a):
        return up(np.sqrt(up(up(a * a) + up(c2 * c2))))

    return ssum_down(down(-c1 - hyp(np.maximum(-2.0 * c1, 0.0))), axis=1), ssum_up(up(hyp(np.maximum(2.0 * c1, 0.0)) - c1), axis=1)


def _best_point(C, half_disks):
    """A maximiser of C_j' eps over the set, in plain floating point (it only chooses the tangents)."""
    if not half_disks:
        return np.where(C >= 0, 1.0, -1.0)
    n = C.shape[1] // 2
    a, b = np.maximum(2.0 * C[:, :n], 0.0), C[:, n:]
    nrm = np.hypot(a, b)
    u = np.divide(a, nrm, out=np.zeros_like(a), where=nrm > 0)
    w = np.divide(b, nrm, out=np.zeros_like(b), where=nrm > 0)
    return np.hstack([2.0 * u - 1.0, w])


def joint_range(C, H, half_disks: bool = False, refinements: int = 4):
    """Rigorous (lo, hi) of C_j' eps + eps' H_j eps over the set, linear and quadratic terms together (the
    "joint" bound of the module docstring). C (m, K) and H (m, K, K), symmetric, are floating-point data."""
    m, K = C.shape
    lam, W = np.linalg.eigh(H)
    eye = np.eye(K)
    pk = 1.0 + 4.0 * (K + 1) * U
    out = {}
    for sign in (1.0, -1.0):                                          # +1: upper bound; -1: lower bound
        Uf = W * np.sqrt(np.maximum(-sign * lam, 0.0))[:, None, :]    # columns: the tangent directions u_b
        aU = np.abs(Uf)
        rest = H + sign * np.einsum("jkb,jlb->jkl", Uf, Uf)           # H = rest - sign * sum_b u_b u_b'
        e_rest = up(err_const(K + 1) * add_up(np.abs(H), up(np.einsum("jkb,jlb->jkl", aU, aU) * pk)))
        ct = C
        for _ in range(refinements + 1):
            a = np.einsum("jkb,jk->jb", Uf, _best_point(sign * ct, half_disks))
            ct = C - sign * 2.0 * np.einsum("jkb,jb->jk", Uf, a)
        e_ct = up(err_const(K + 1) * add_up(np.abs(C), up(2.0 * up(np.einsum("jkb,jb->jk", aU, np.abs(a)) * pk))))
        lo_l, hi_l = _lin_range(ct, half_disks)
        d = np.einsum("jkk->jk", rest)
        const = add_up(sum_up(up(a * a), axis=1), sum_up(e_ct, axis=1), sum_up(np.abs(rest) * (1.0 - eye), axis=(1, 2)), sum_up(e_rest, axis=(1, 2)))
        if sign > 0:
            out["hi"] = add_up(hi_l, const, sum_up(np.maximum(d, 0.0), axis=1))
        else:
            out["lo"] = -add_up(-lo_l, const, sum_up(np.maximum(-d, 0.0), axis=1))
    return out["lo"], out["hi"]


def magnitude_bounds(Vn, E, B, z, half_disks: bool = False, max_k: int | None = None, node_budget: int = 0, joint: bool = True,
                     details: bool = False):
    """(lo, hi) of |V_j| for V_j = Vn_j + E_j eps + eps' B_j eps + s_j, |s_j| <= z_j, eps in the box (or, with
    ``half_disks``, in the product of half disks whose bounding boxes are the box: P symbols, then Q symbols).
    ``joint``: also bound the linear and quadratic terms together and keep the better bound on each side.
    ``details``: return in addition the share of buses at which the joint bound is the better one (lo, hi)."""
    m, K = E.shape
    vr, vi = Vn.real, Vn.imag
    vup, vlo = abs_up(Vn), abs_down(Vn)
    va = np.sqrt(vr * vr + vi * vi)
    alpha = sum_up(abs_up(E), axis=1)
    b_l1 = sum_up(abs_up(B), axis=(1, 2))
    # components along (cl, Hp) and across (cp, Hq) v; every entry within 16 u |E_jk| (|B_jkl|) of its exact value
    cl = (vr[:, None] * E.real + vi[:, None] * E.imag) / va[:, None]
    cp = (vr[:, None] * E.imag - vi[:, None] * E.real) / va[:, None]
    Hp = (vr[:, None, None] * B.real + vi[:, None, None] * B.imag) / va[:, None, None]
    Hq = (vr[:, None, None] * B.imag - vi[:, None, None] * B.real) / va[:, None, None]
    H = Hp + cp[:, :, None] * cp[:, None, :] / (2.0 * va)[:, None, None]
    H = 0.5 * (H + H.transpose(0, 2, 1))
    s_lin = up(16.0 * U * alpha)                                                   # sum_k |cl_jk - exact|
    s_quad = up(32.0 * U * add_up(b_l1, div_up(2.0 * mul_up(alpha, alpha), vlo)))   # sum_kl |H_jkl - exact|

    lin_lo, lin_hi = _lin_range(cl, half_disks)
    k = K if max_k is None else max_k
    q = np.array([raak_range(h, k, node_budget)[:2] for h in H]).reshape(m, 2)
    f_hi, f_lo = add_up(lin_hi, q[:, 1]), down(lin_lo + q[:, 0])              # separate bound of c' eps + eps' H eps
    share = (0.0, 0.0)
    if joint and K > 0:
        j_lo, j_hi = joint_range(cl, H, half_disks)
        share = (float(np.mean(j_lo > f_lo)), float(np.mean(j_hi < f_hi)))
        f_hi, f_lo = np.minimum(f_hi, j_hi), np.maximum(f_lo, j_lo)
    slack = add_up(s_lin, s_quad)                                             # rounding of the coefficients
    f_hi, f_lo = add_up(f_hi, slack), sub_down(f_lo, slack)

    Q1 = add_up(sum_up(np.abs(cp), axis=1), s_lin)
    P1 = add_up(sum_up(np.abs(cl), axis=1), s_lin)
    slack2 = up(16.0 * U * b_l1)
    P2, Q2 = add_up(_form_sup(Hp), slack2), add_up(_form_sup(Hq), slack2)
    Pm = add_up(P1, P2)
    den = sub_down(vlo, Pm)
    ok = den > 0
    den = np.where(ok, den, 1.0)
    Qs = add_up(Q1, Q2)
    Qs2 = mul_up(Qs, Qs)
    e3 = add_up(div_up(mul_up(mul_up(Q1, Q1), Pm), down(2.0 * down(vlo * den))),
                div_up(add_up(mul_up(2.0 * Q1, Q2), mul_up(Q2, Q2)), 2.0 * den),
                div_up(mul_up(Qs2, Qs2), down(8.0 * down(down(den * den) * den))))
    hi = add_up(vup, f_hi, e3, z)
    lo = sub_down(sub_down(down(vlo + f_lo), e3), z)
    lo, hi = np.where(ok, lo, np.nan), np.where(ok, hi, np.nan)
    return (lo, hi, share) if details else (lo, hi)


@dataclass
class VerifiedEnclosure:
    """Bounds of |V_i| for every PQ bus over the whole uncertainty set, rigorous in floating point."""

    lo: np.ndarray
    hi: np.ndarray
    verified: bool              # a solution exists inside the bounds for every point of the set
    unique: bool                # and it is the only one within ``s_abs`` of the second-order model
    s_abs: np.ndarray           # radius of the tube around the second-order model, per PQ bus [p.u.]
    iterations: int
    reason: str
    rounding_share: float
    joint_share: tuple = (0.0, 0.0)     # share of buses whose lower / upper bound comes from the joint bound


def verified_enclosure(so: SecondOrder, half_disks: bool = False, max_k: int | None = None, node_budget: int = 0,
                       Hs: np.ndarray | None = None, s_err=None, joint: bool = True) -> VerifiedEnclosure:
    """Verified enclosure of the voltage magnitudes of a network without PV buses over the box of ``so``
    (over the half disks inscribed in it with ``half_disks``; see cdhim.halfdisk). ``s_err``: see ``verify_tube``."""
    g, st = so.st.grid, so.st
    if len(g.pv):
        raise ValueError("the verified enclosure needs all non-slack buses to be PQ")
    ns, m = g.pq, so.m
    Hs = state_forms(so) if Hs is None else Hs
    E, B = so.E[ns], Hs[:m] + 1j * Hs[m:]
    tube = verify_tube(g.Ybus, ns, st.V_star, st.sbus[ns], injection_changes(so)[ns], E, B, s_err)
    nan = np.full(m, np.nan)
    if not tube.exists:
        return VerifiedEnclosure(nan, nan, False, False, nan, tube.iterations, tube.reason, np.nan)
    lo, hi, share = magnitude_bounds(st.V_star[ns], E, B, tube.z, half_disks, max_k, node_budget, joint, details=True)
    ok = bool(np.all(np.isfinite(lo)))
    return VerifiedEnclosure(lo, hi, ok, ok and tube.unique, tube.z, tube.iterations,
                             "" if ok else "the bound of the modulus breaks down", tube.rounding_share, share)


def verified_point(grid, sbus, V, s_err=None, cache: dict | None = None):
    """Rigorous (lo, hi, unique) of |V_i| at every PQ bus for the power-flow solution next to the computed ``V``
    with the specified injections ``sbus`` (known within ``s_err`` per PQ bus): the tube with no uncertainty.
    NaN if the verification fails. ``cache``: see ``verify_tube``; if the shared inverse is not good enough for
    this point, the verification is repeated with its own."""
    ns = grid.pq
    m = len(ns)
    Z = np.zeros((m, 0), complex)
    args = (grid.Ybus, ns, V, np.asarray(sbus, dtype=complex)[ns], Z, Z, np.zeros((m, 0, 0), complex), s_err)
    tube = verify_tube(*args, cache=cache)
    if cache is not None and not (tube.exists and tube.unique):
        tube = verify_tube(*args)
    if not tube.exists:
        return np.full(m, np.nan), np.full(m, np.nan), False
    return sub_down(abs_down(V[ns]), tube.z), add_up(abs_up(V[ns]), tube.z), tube.unique
