"""Verified tube around the second-order model for networks with PV buses: voltage and current remainders.

The equations of the network are written with the injected currents I = Y V + c as they appear in every bus:

    PQ bus j     N_j = conj(V_j) I_j - conj(S_j(eps))                       (complex)
    PV bus j     N_j = ( Re(conj(V_j) I_j) - P_j(eps),  |V_j|^2 - v_j^2 )   (two real equations)

so the non-linearity is a product of the voltage and the current of the SAME bus. With the model m(eps) of the
voltages and s the remainder, t = Y s is the remainder of the currents and, exactly,

    N(m + s) = N(m) + Jc[s] + rem(s),      Jc[s]_j = conj(s_j) I*_j + conj(V*_j) t_j   (2 Re(conj(V*_j) s_j) in the
                                                                                         second PV equation),
    PQ and first PV equation   |rem_j| <= |s_j| |t_j| + |a_j + b_j| |t_j| + |s_j| |(Y(a + b))_j|,
    second PV equation         |rem_j| <= |s_j|^2 + 2 |a_j + b_j| |s_j|,

a = E eps and b = eps' B eps being the first- and second-order parts of the model. No entry of Y multiplies a
bound of s: the current remainder is bounded on its own. On the set

    D(z, zeta) = { s : |s_j| <= z_j,  |(Y s)_j| <= zeta_j },

convex and compact, the operator Psi(s) = s - C N(m + s) (C any real 2m x 2m matrix, meant to be Jc^-1) gives

    |Psi(s)_j|     <= (R^V z)_j + (T^V rho)_j,          R^V, T^V: 2x2 block norms of I - C Jc and of C,
    |(Y Psi(s))_j| <= (R^I z)_j + (T^I rho)_j,          R^I, T^I: of Y (I - C Jc) and of Y C,
    rho_j = phi_j + z_j zeta_j + dl_j zeta_j + Ibar_j z_j  (+ z_j^2 + 2 dl_j z_j at PV buses),

phi_j >= |N_j(m(eps))| (defect of the model), dl_j >= |a_j + b_j|, Ibar_j >= |(Y(a + b))_j|. If both right-hand
sides are below z and zeta, Psi maps D into itself and Brouwer's theorem gives a zero of C N in it; R^V z < z
makes C non-singular. The same inequalities with the Lipschitz bound of rem (2 z zeta instead of z zeta, 2 z^2
instead of z^2) make Psi a contraction in max(|ds_j|/z_j, |(Y ds)_j|/zeta_j), hence uniqueness in D.

The defect is the polynomial of cdhim.verified (degrees 0 to 2 evaluated with error bounds, degrees 3 and 4
bounded with the norms of Y E and Y B); at PV buses only its real part is an equation, and the second equation
has the defect |m_j|^2 - v_j^2, of degrees 0 to 4 as well.

Rounding: as in cdhim.verified. Jc is never formed with rounded products: C Jc is evaluated as (C Dv) Yr + C Kc
with Dv and Kc made of stored numbers, and the error of every product is bounded.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .enclosure import injection_changes, state_forms
from .second_order import SecondOrder
from .verified import (VerifiedEnclosure, _cform_sup, _sparse_abs, abs_down, abs_up, add_up, csum_up, down, err_const, magnitude_bounds,
                       matmul_up, mul_up, sum_up, up)


def _blocks(A, m):
    return A[:m, :m], A[:m, m:], A[m:, :m], A[m:, m:]


def block_frobenius_up(A, m):
    """Upper bound of the spectral norm of every 2x2 block of a matrix known through entrywise bounds A >= 0."""
    a, b, c, d = _blocks(A, m)
    return up(np.sqrt(add_up(mul_up(a, a), mul_up(b, b), mul_up(c, c), mul_up(d, d))))


def block_spectral_up(A, m):
    """Upper bound of the spectral norm of every 2x2 block of the floating-point matrix A."""
    a, b, c, d = _blocks(A, m)
    s = add_up(up(a * a), up(b * b), up(c * c), up(d * d))
    e_det = up(err_const(2) * add_up(up(np.abs(a * d)), up(np.abs(b * c))))
    det_lo = np.maximum(down(np.abs(a * d - b * c) - e_det), 0.0)
    disc = np.maximum(up(up(s * s) - 4.0 * np.maximum(down(det_lo * det_lo), 0.0)), 0.0)
    return up(np.sqrt(up(0.5 * up(s + up(np.sqrt(disc))))))


@dataclass
class TubeVI:
    exists: bool
    unique: bool
    z: np.ndarray               # radius of the voltage remainder per non-slack bus [p.u.]
    zeta: np.ndarray            # radius of the current remainder per non-slack bus [p.u.]
    iterations: int
    reason: str
    rounding_share: float
    aux: dict | None = None


def verify_tube_vi(Ybus, ns, is_pv, v_set, V, S0, M, E, B, s_err=None, inflate=(1e-9, 1e-6, 1e-3), max_iter: int = 5000,
                   keep: bool = False) -> TubeVI:
    """Radii (z, zeta) of the set around m(eps) = V[ns] + E eps + eps' B eps, verified with directed bounds.

    ``ns``: the non-slack buses; ``is_pv`` (m,) marks the voltage-controlled ones and ``v_set`` (m,) holds their
    set-points; of ``S0`` and ``M`` only the real part is used at PV buses. Other arguments as in
    cdhim.verified.verify_tube.
    """
    Y = Ybus.tocsr()
    m, K = E.shape
    is_pv = np.asarray(is_pv, dtype=bool)
    pqf = (~is_pv).astype(float)
    Yn = Y[ns].tocsr()
    YL = Yn[:, ns].tocsr()
    p = max(int(np.diff(Yn.indptr).max()), 1)
    aYn, aYL = _sparse_abs(Yn), _sparse_abs(YL)
    Vn, Vc = V[ns], np.conj(V[ns])
    vr, vi = Vn.real, Vn.imag
    vup = abs_up(Vn)
    nan = np.full(m, np.nan)
    S0 = np.where(is_pv, np.asarray(S0).real, S0)
    M = np.where(is_pv[:, None], np.asarray(M).real, M)

    def part(g, mask):
        """Modulus at PQ buses, |real part| at PV buses."""
        return np.where(mask, np.abs(g.real), abs_up(g))

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

    # first equation of every bus: conj(m)(Y m + c) - conj(S), real part only at PV buses
    g0 = Vc * cur - np.conj(S0)
    G0 = add_up(part(g0, is_pv), mul_up(vup, e_cur), up(err_const(2) * add_up(mul_up(vup, a_cur), aS0)))
    if s_err is not None:
        G0 = add_up(G0, s_err)
    g1 = Vc[:, None] * YE + np.conj(E) * cur[:, None] - np.conj(M)
    e1 = add_up(mul_up(vup[:, None], e_YE), mul_up(aE, e_cur[:, None]),
                up(err_const(3) * add_up(mul_up(vup[:, None], aYE), mul_up(aE, a_cur[:, None]), aM)))
    G1 = sum_up(add_up(part(g1, is_pv[:, None]), e1), axis=1)
    cross = np.conj(E)[:, :, None] * YE[:, None, :]
    g2 = Vc[:, None, None] * YB + np.conj(B) * cur[:, None, None] + 0.5 * (cross + cross.transpose(0, 2, 1))
    a_cross, e_cross = mul_up(aE[:, :, None], aYE[:, None, :]), mul_up(aE[:, :, None], e_YE[:, None, :])
    e2 = add_up(mul_up(vup[:, None, None], e_YB), mul_up(aB, e_cur[:, None, None]), 0.5 * add_up(e_cross, e_cross.transpose(0, 2, 1)),
                up(err_const(4) * add_up(mul_up(vup[:, None, None], aYB), mul_up(aB, a_cur[:, None, None]),
                                         0.5 * add_up(a_cross, a_cross.transpose(0, 2, 1)))))
    G2 = sum_up(add_up(part(g2, is_pv[:, None, None]), e2), axis=(1, 2))
    rot = Vc / np.sqrt(vr * vr + vi * vi)                              # axes along and across V*, as in cdhim.verified
    alpha = csum_up(E, 1, rot)
    beta = _cform_sup(B, rot)
    YEn = add_up(csum_up(YE, 1, rot), sum_up(e_YE, axis=1))
    YBn = add_up(_cform_sup(YB, rot), sum_up(e_YB, axis=(1, 2)))
    G3 = add_up(mul_up(alpha, YBn), mul_up(beta, YEn), mul_up(beta, YBn))
    low, high = add_up(G0, G1, G2), G3

    # second equation of the PV buses: |m|^2 - v_set^2
    if is_pv.any():
        h0 = vr * vr + vi * vi - v_set * v_set
        H0 = add_up(np.abs(h0), up(err_const(3) * add_up(up(vr * vr), up(vi * vi), up(v_set * v_set))))
        h1 = 2.0 * (vr[:, None] * E.real + vi[:, None] * E.imag)
        H1 = sum_up(add_up(np.abs(h1), up(err_const(2) * 2.0 * mul_up(vup[:, None], aE))), axis=1)
        ee = E.real[:, :, None] * E.real[:, None, :] + E.imag[:, :, None] * E.imag[:, None, :]
        h2 = 2.0 * (vr[:, None, None] * B.real + vi[:, None, None] * B.imag) + ee
        H2 = sum_up(add_up(np.abs(h2), up(err_const(4) * add_up(2.0 * mul_up(vup[:, None, None], aB), mul_up(aE[:, :, None], aE[:, None, :])))), axis=(1, 2))
        H3 = add_up(2.0 * mul_up(alpha, beta), mul_up(beta, beta))
        low = np.where(is_pv, add_up(low, H0, H1, H2), low)
        high = np.where(is_pv, add_up(high, H3), high)
    phi = add_up(low, high)                                            # >= 2-norm of the defect of the bus
    delta = add_up(alpha, beta)                                        # >= |a + b|
    Ibar = add_up(YEn, YBn, e_cur)                                     # >= |Y (a + b)|, plus the error of I* in Kc

    # Jc = Dv Yr + Kc in real form (rows: first equations, second equations; columns: Re s, Im s) and C ~ Jc^-1
    Gs, Bs = YL.real, YL.imag
    Yr = sp.bmat([[Gs, -Bs], [Bs, Gs]], format="csr")
    d11, d12, d21, d22 = vr, vi, -vi * pqf, vr * pqf
    ir, ii = cur.real, cur.imag
    k11, k12, k21, k22 = ir, ii, np.where(is_pv, 2.0 * vr, ii), np.where(is_pv, 2.0 * vi, -ir)
    Dv = sp.bmat([[sp.diags(d11), sp.diags(d12)], [sp.diags(d21), sp.diags(d22)]], format="csr")
    Kc = sp.bmat([[sp.diags(k11), sp.diags(k12)], [sp.diags(k21), sp.diags(k22)]], format="csr")
    C = np.linalg.inv((Dv @ Yr + Kc).toarray())
    aC, aYr = np.abs(C), abs(Yr)
    pc = max(int(np.diff(Yr.tocsc().indptr).max()), 1)
    pr = max(int(np.diff(Yr.indptr).max()), 1)

    def two_terms(A, c1, c2, c3, c4):
        """A [[diag c1, diag c2], [diag c3, diag c4]] and the bound of its rounding error."""
        a1, a2 = A[:, :m], A[:, m:]
        val = np.hstack([a1 * c1 + a2 * c3, a1 * c2 + a2 * c4])
        b1, b2 = np.abs(a1), np.abs(a2)
        err = up(err_const(2) * np.hstack([add_up(mul_up(b1, np.abs(c1)), mul_up(b2, np.abs(c3))), add_up(mul_up(b1, np.abs(c2)), mul_up(b2, np.abs(c4)))]))
        return val, err

    CD, e_CD = two_terms(C, d11, d12, d21, d22)
    P1 = (Yr.T.tocsr() @ CD.T).T
    e_P1 = add_up(up(err_const(pc) * matmul_up(np.abs(CD), aYr, pc)), matmul_up(e_CD, aYr, pc))
    P2, e_P2 = two_terms(C, k11, k12, k21, k22)
    eye = np.eye(2 * m)
    res = add_up(np.abs(eye - P1 - P2), e_P1, e_P2, up(err_const(3) * add_up(eye, np.abs(P1), np.abs(P2))))
    RV = block_frobenius_up(res, m)
    TV = block_spectral_up(C, m)
    YC = Yr @ C
    TI = add_up(block_spectral_up(YC, m), block_frobenius_up(up(err_const(pr) * matmul_up(aYr, aC, pr)), m))
    RI = block_frobenius_up(matmul_up(aYr, res, pr), m)

    def rho(z, zeta, lipschitz=False):
        f = 2.0 if lipschitz else 1.0
        r = add_up(mul_up(f * z, zeta), mul_up(delta, zeta), mul_up(Ibar, z))
        return np.where(is_pv, add_up(r, mul_up(f * z, z), mul_up(2.0 * delta, z)), r)

    z, zeta = TV @ phi, TI @ phi
    it, reason = 0, "the iteration did not settle"
    for it in range(1, max_iter + 1):
        r = phi + rho(z, zeta)
        nz, nt = RV @ z + TV @ r, RI @ z + TI @ r
        if not (np.all(np.isfinite(nz)) and np.all(np.isfinite(nt))) or nz.max() > 1e6:
            reason = "the iteration diverged"
            break
        settled = max(np.abs(nz - z).max() / max(nz.max(), 1e-300), np.abs(nt - zeta).max() / max(nt.max(), 1e-300)) <= 1e-12
        z, zeta = nz, nt
        if settled:
            reason = ""
            break
    if reason:
        return TubeVI(False, False, nan, nan, it, reason, np.nan)

    for f in inflate:
        zz, ze = up(z * (1.0 + f)), up(zeta * (1.0 + f))
        Rz, Ri = matmul_up(RV, zz, m), matmul_up(RI, zz, m)
        r = add_up(phi, rho(zz, ze))
        if (np.all(add_up(Rz, matmul_up(TV, r, m)) <= zz) and np.all(add_up(Ri, matmul_up(TI, r, m)) <= ze) and np.all(Rz < zz)):
            lam = rho(zz, ze, True)
            unique = bool(np.all(add_up(Rz, matmul_up(TV, lam, m)) < zz) and np.all(add_up(Ri, matmul_up(TI, lam, m)) < ze))
            share = float(np.max(matmul_up(TV, low, m) / zz))
            return TubeVI(True, unique, zz, ze, it, "", share, dict(C=C, TV=TV, TI=TI, RV=RV, RI=RI, phi=phi) if keep else None)
    return TubeVI(False, False, nan, nan, it, "the inflated fixed point is not mapped into itself", np.nan)


def verified_enclosure_vi(so: SecondOrder, half_disks: bool = False, max_k: int | None = None, node_budget: int = 0,
                          Hs: np.ndarray | None = None, s_err=None) -> VerifiedEnclosure:
    """Verified enclosure of |V_i| at the PQ buses of a network that may have PV buses (voltage-current form)."""
    g, st = so.st.grid, so.st
    ns, m, npv = g.pvpq, so.m, len(g.pv)
    is_pv = np.r_[np.ones(npv, bool), np.zeros(m - npv, bool)]
    Hs = state_forms(so) if Hs is None else Hs
    E, B = so.E[ns], Hs[:m] + 1j * Hs[m:]
    tube = verify_tube_vi(g.Ybus, ns, is_pv, np.asarray(g.v_set, dtype=float)[ns], st.V_star, st.sbus[ns], injection_changes(so)[ns], E, B, s_err)
    nan = np.full(m - npv, np.nan)
    if not tube.exists:
        return VerifiedEnclosure(nan, nan, False, False, nan, tube.iterations, tube.reason, np.nan)
    lo, hi = magnitude_bounds(st.V_star[g.pq], E[npv:], B[npv:], tube.z[npv:], half_disks, max_k, node_budget)
    ok = bool(np.all(np.isfinite(lo)))
    return VerifiedEnclosure(lo, hi, ok, ok and tube.unique, tube.z[npv:], tube.iterations,
                             "" if ok else "the bound of the modulus breaks down", tube.rounding_share)
