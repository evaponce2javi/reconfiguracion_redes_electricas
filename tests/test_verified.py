"""The verified enclosure: the directed bounds against exact rational arithmetic, the pieces of the derivation
against exact and sampled evaluations, and the result against AC solutions."""
import warnings
from fractions import Fraction as Fr

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import injection_changes, state_forms, voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.powerflow import newton_pf
from cdhim.verified import (U, abs_down, abs_up, err_const, magnitude_bounds, matmul_up, raak_range, ssum_down, ssum_up, sum_up,
                            verified_enclosure, verified_point, verify_tube)


def _close_ties(net):
    net.line["in_service"] = True


def _add_shunts(net):
    import pandapower as pp

    net.line["c_nf_per_km"] = 300.0
    for bus in (17, 32):
        pp.create_shunt(net, bus, q_mvar=-0.3)


def _grid(case, prepare=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case, prepare=prepare)


CASES = [("case33bw", None, "F", 0.5), ("case33bw", _close_ties, "R", 0.8), ("case33bw", _add_shunts, "F", 0.5)]


@pytest.fixture(scope="module", params=CASES, ids=["radial", "meshed", "shunts"])
def model(request):
    case, prepare, scen, rho_U = request.param
    g = _grid(case, prepare)
    sc = SCENARIOS[scen]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], rho_U)
    st = half_disk_study(g, s0, cap)
    so = SecondOrder(st, 1.0)
    Hs = state_forms(so)
    m = so.m
    data = dict(ns=g.pq, V=st.V_star, S0=st.sbus[g.pq], M=injection_changes(so)[g.pq], E=so.E[g.pq], B=Hs[:m] + 1j * Hs[m:])
    return g, s0, cap, st, so, data


# ----------------------------------------------------------------------------- directed bounds, exactly
def test_sums_bound_the_exact_rational_sum():
    rng = np.random.default_rng(0)
    for n in (1, 7, 1000):
        x = np.abs(rng.standard_normal(n)) * 10.0 ** rng.uniform(-12, 3, n)
        exact = sum(Fr(v) for v in x)
        assert Fr(float(sum_up(x))) >= exact
        y = x * rng.choice([-1.0, 1.0], n)
        exact = sum(Fr(v) for v in y)
        assert Fr(float(ssum_down(y))) <= exact <= Fr(float(ssum_up(y)))
    # heavy cancellation: the bracket must still hold
    y = np.array([1e16, 1.0, -1e16, 1e-8, 3.0, -4.0])
    assert Fr(float(ssum_down(y))) <= sum(Fr(v) for v in y) <= Fr(float(ssum_up(y)))


def test_modulus_bounds_are_exact_brackets():
    rng = np.random.default_rng(1)
    z = (rng.standard_normal(500) + 1j * rng.standard_normal(500)) * 10.0 ** rng.uniform(-9, 6, 500)
    lo, hi = abs_down(z), abs_up(z)
    for v, a, b in zip(z, lo, hi):
        sq = Fr(v.real) ** 2 + Fr(v.imag) ** 2
        assert Fr(float(a)) ** 2 <= sq <= Fr(float(b)) ** 2


def test_products_and_their_error_bound_exactly():
    rng = np.random.default_rng(2)
    n, p, q = 4, 9, 3
    A, B = np.abs(rng.standard_normal((n, p))), np.abs(rng.standard_normal((p, q))) * 10.0 ** rng.uniform(-6, 6, (p, q))
    P = matmul_up(A, B, p)
    for i in range(n):
        for j in range(q):
            assert Fr(float(P[i, j])) >= sum(Fr(A[i, k]) * Fr(B[k, j]) for k in range(p))
    # complex inner products with cancellation: |fl - exact| <= err_const(p) * sum |x||y|
    x = (rng.standard_normal((50, p)) + 1j * rng.standard_normal((50, p))) * 10.0 ** rng.uniform(-3, 6, (50, p))
    y = (rng.standard_normal((p, 50)) + 1j * rng.standard_normal((p, 50)))
    y[-1] = -(x[:, :-1] * y[:-1].T).sum(axis=1) / x[:, -1]              # row r of x against column r of y nearly cancels
    for r in range(50):
        fl = x[r] @ y[:, r]
        re = sum(Fr(a.real) * Fr(b.real) - Fr(a.imag) * Fr(b.imag) for a, b in zip(x[r], y[:, r]))
        im = sum(Fr(a.real) * Fr(b.imag) + Fr(a.imag) * Fr(b.real) for a, b in zip(x[r], y[:, r]))
        bound = Fr(err_const(p)) * sum(Fr(float(abs_up(a))) * Fr(float(abs_up(b))) for a, b in zip(x[r], y[:, r]))
        assert (Fr(fl.real) - re) ** 2 + (Fr(fl.imag) - im) ** 2 <= bound ** 2


def test_quadratic_range_holds_exactly_at_the_vertices():
    rng = np.random.default_rng(3)
    for k, scale in ((5, 1.0), (6, 1e-7), (4, 1e5)):
        A = rng.standard_normal((k, k))
        H = scale * 0.5 * (A + A.T)
        lo, hi, _ = raak_range(H, k)
        for v in range(2 ** k):
            e = [1 if (v >> j) & 1 else -1 for j in range(k)]
            val = sum(Fr(H[a, b]) * e[a] * e[b] for a in range(k) for b in range(k))
            assert Fr(lo) <= val <= Fr(hi)
        # and it is no wider than the term-by-term bound by more than rounding
        d = np.diag(H)
        off = np.abs(H).sum() - np.abs(d).sum()
        assert hi - lo <= (np.abs(d).sum() + 2 * off) * (1 + 1e-9)


# ----------------------------------------------------------------------------- the derivation
class Q2:
    """Complex numbers with rational parts: exact +, -, *, /, conj."""

    def __init__(self, re, im=0):
        self.re, self.im = Fr(re), Fr(im)

    @classmethod
    def of(cls, z):
        z = complex(z)
        return cls(z.real, z.imag)

    def __add__(self, o):
        return Q2(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return Q2(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        return Q2(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def conj(self):
        return Q2(self.re, -self.im)

    def __truediv__(self, o):
        n = o.re ** 2 + o.im ** 2
        t = self * o.conj()
        return Q2(t.re / n, t.im / n)

    def abs2(self):
        return self.re ** 2 + self.im ** 2


def test_defect_bound_holds_in_exact_arithmetic(model):
    """phi_j >= |F_j(m(eps), eps)| at vertices of the box, with F evaluated in rational arithmetic."""
    g, s0, cap, st, so, d = model
    tube = verify_tube(g.Ybus, keep=True, **d)
    assert tube.exists
    ns, m, K = d["ns"], *d["E"].shape
    Y = g.Ybus.tocsr()
    pos = np.full(g.n_bus, -1)
    pos[ns] = np.arange(m)
    Vq = [Q2.of(v) for v in d["V"]]
    rng = np.random.default_rng(4)
    for _ in range(3):
        eps = [int(e) for e in rng.choice([-1, 1], K)]
        mod = []
        for j in range(m):
            v = Q2.of(d["V"][ns[j]])
            for k in range(K):
                v = v + Q2.of(d["E"][j, k] * eps[k])
                for l in range(K):
                    v = v + Q2.of(d["B"][j, k, l]) * Q2(eps[k] * eps[l])
            mod.append(v)
        for j in rng.choice(m, 6, replace=False):
            cur = Q2(0)
            for idx in range(Y.indptr[ns[j]], Y.indptr[ns[j] + 1]):
                c = Y.indices[idx]
                cur = cur + Q2.of(Y.data[idx]) * (mod[pos[c]] if pos[c] >= 0 else Vq[c])
            S = Q2.of(d["S0"][j])
            for k in range(K):
                S = S + Q2.of(d["M"][j, k]) * Q2(eps[k])
            F = cur - S.conj() / mod[j].conj()
            assert F.abs2() <= Fr(float(tube.aux["phi"][j])) ** 2


def test_operator_maps_the_tube_into_itself(model):
    """Phi(s) = s - C F(m + s) evaluated at sampled eps and s on the boundary of the tube stays inside it."""
    g, s0, cap, st, so, d = model
    tube = verify_tube(g.Ybus, keep=True, **d)
    assert tube.exists and tube.unique
    ns, m, K = d["ns"], *d["E"].shape
    Y, C, z = g.Ybus.tocsr(), tube.aux["C"], tube.z
    rng = np.random.default_rng(5)
    worst = 0.0
    for trial in range(200):
        eps = rng.choice([-1.0, 1.0], K) if trial % 2 else rng.uniform(-1, 1, K)
        mod = d["V"][ns] + d["E"] @ eps + np.einsum("jkl,k,l->j", d["B"], eps, eps)
        s = z * np.exp(1j * rng.uniform(0, 2 * np.pi, m))
        V = d["V"].copy()
        V[ns] = mod + s
        F = (Y @ V)[ns] - np.conj(d["S0"] + d["M"] @ eps) / np.conj(V[ns])
        CF = C @ np.r_[F.real, F.imag]
        phi_s = s - (CF[:m] + 1j * CF[m:])
        worst = max(worst, float((np.abs(phi_s) / z).max()))
    assert worst <= 1.0


def test_tube_contains_the_ac_solution(model):
    g, s0, cap, st, so, d = model
    tube = verify_tube(g.Ybus, **d)
    K = d["E"].shape[1]
    eps = np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** K)] for j in range(K)])
    V, conv, _ = st.solver.solve(st.injections(eps, 1.0))
    assert conv.all()
    mod = d["V"][g.pq][:, None] + d["E"] @ eps + np.einsum("jkl,kv,lv->jv", d["B"], eps, eps)
    assert np.all(np.abs(V[g.pq] - mod) <= tube.z[:, None] + 1e-9)          # 1e-9: tolerance of the AC solver
    assert tube.rounding_share < 1e-3


# ----------------------------------------------------------------------------- the result
def test_enclosure_contains_ac_solutions_over_the_half_disks(model):
    g, s0, cap, st, so, d = model
    enc = verified_enclosure(so, half_disks=True)
    assert enc.verified and enc.unique
    dg = np.flatnonzero(cap > 0)
    pts = sample_half_disks(np.random.default_rng(6), len(dg), 400)
    S = np.repeat(np.asarray(s0, dtype=complex)[:, None], pts.shape[1], axis=1)
    S[dg] += cap[dg, None] * pts
    V, conv, _ = st.solver.solve(S)
    vm = np.abs(V[g.pq][:, conv])
    assert conv.sum() > 300
    assert np.all(vm >= enc.lo[:, None] - 1e-9) and np.all(vm <= enc.hi[:, None] + 1e-9)
    for sign in (+1, -1):
        att = attained_extreme(g, s0, cap, sign)
        assert np.all(att >= enc.lo - 1e-9) and np.all(att <= enc.hi + 1e-9)
    # close to the enclosure without directed rounding (a different bound of the remainder, same model)
    ref = voltage_enclosure(so, first_order=first_order_range(so.c_lin))
    assert ref.converged
    ratio = (enc.hi - enc.lo) / (ref.hi - ref.lo)
    assert 0.8 < ratio.min() and ratio.max() < 1.2
    box = verified_enclosure(so)
    assert np.all(box.lo <= enc.lo) and np.all(box.hi >= enc.hi)


def test_magnitude_bounds_hold_for_points_of_the_tube(model):
    g, s0, cap, st, so, d = model
    z = np.full(len(g.pq), 1e-4)
    lo, hi = magnitude_bounds(d["V"][g.pq], d["E"], d["B"], z)
    rng = np.random.default_rng(7)
    K = d["E"].shape[1]
    for _ in range(300):
        eps = rng.choice([-1.0, 1.0], K)
        v = d["V"][g.pq] + d["E"] @ eps + np.einsum("jkl,k,l->j", d["B"], eps, eps) + z * np.exp(1j * rng.uniform(0, 2 * np.pi, len(z)))
        assert np.all(np.abs(v) >= lo) and np.all(np.abs(v) <= hi)


def test_point_verification(model):
    g, s0, cap, st, so, d = model
    lo, hi, unique = verified_point(g, st.sbus, st.V_star)
    vm = np.abs(st.V_star[g.pq])
    assert unique and np.all(lo <= vm) and np.all(vm <= hi) and (hi - lo).max() < 1e-7
    # a point that is not a solution is not certified as one within a small radius
    V, ok, _, _ = newton_pf(g.Ybus, st.sbus, g.V0, g.pv, g.pq, tol=1e-9, max_iter=40)
    assert ok
    bad = V.copy()
    bad[g.pq] *= 1.01
    lo2, hi2, _ = verified_point(g, st.sbus, bad)
    assert not np.all(np.isfinite(lo2)) or (hi2 - lo2).max() > 1e-3


def test_pv_buses_are_rejected():
    g = _grid("case30")
    with pytest.raises(ValueError):
        from cdhim import make_study

        verified_enclosure(SecondOrder(make_study(g), 0.01))


# ----------------------------------------------------------------------------- linear and quadratic terms together
def _random_forms(rng, k, scale, definite=0):
    A = rng.standard_normal((k, k))
    H = 0.5 * (A + A.T) if definite == 0 else definite * (A @ A.T)
    return scale * rng.standard_normal(k), scale * scale * H


def test_joint_range_holds_exactly_at_the_vertices_of_the_box():
    from cdhim.verified import joint_range

    rng = np.random.default_rng(8)
    for k, scale, definite in ((5, 1.0, 0), (6, 1e-3, -1), (4, 30.0, 1), (5, 1e-2, 0)):
        c, H = _random_forms(rng, k, scale, definite)
        lo, hi = joint_range(c[None, :], H[None, :, :])
        for v in range(2 ** k):
            e = [1 if (v >> j) & 1 else -1 for j in range(k)]
            val = sum(Fr(c[a]) * e[a] for a in range(k)) + sum(Fr(H[a, b]) * e[a] * e[b] for a in range(k) for b in range(k))
            assert Fr(float(lo[0])) <= val <= Fr(float(hi[0]))


def test_joint_range_holds_inside_the_box_and_over_half_disks():
    from cdhim.verified import joint_range

    rng = np.random.default_rng(9)
    for half in (False, True):
        for definite in (0, -1, 1):
            c, H = _random_forms(rng, 6, 0.3, definite)
            lo, hi = joint_range(c[None, :], H[None, :, :], half)
            if half:
                pts = sample_half_disks(rng, 3, 20000)
                eps = np.vstack([2 * pts.real - 1, pts.imag])
            else:
                eps = rng.uniform(-1, 1, (6, 20000))
            val = c @ eps + np.einsum("kv,kl,lv->v", eps, H, eps)
            assert lo[0] <= val.min() and val.max() <= hi[0]


def test_joint_range_is_exact_for_a_concave_form_with_a_dominant_linear_term():
    from cdhim.verified import joint_range

    rng = np.random.default_rng(10)
    c, H = _random_forms(rng, 6, 1.0, -1)
    H *= 0.01                                                      # gradient keeps its sign over the box
    lo, hi = joint_range(c[None, :], H[None, :, :])
    e = np.sign(c)
    best = c @ e + e @ H @ e
    assert 0.0 <= hi[0] - best <= 1e-12 * abs(best)
    # the separate bound cannot see that the form is negative at the maximiser of the linear term
    assert np.abs(c).sum() + raak_range(H, 6)[1] > hi[0] + 0.5 * abs(e @ H @ e)


def test_joint_bound_never_loosens_the_magnitude_bounds(model):
    g, s0, cap, st, so, d = model
    z = np.full(len(g.pq), 1e-6)
    for half in (False, True):
        a = magnitude_bounds(d["V"][g.pq], d["E"], d["B"], z, half, joint=False)
        b, share = magnitude_bounds(d["V"][g.pq], d["E"], d["B"], z, half, joint=True, details=True)[:2], None
        assert np.all(b[0] >= a[0]) and np.all(b[1] <= a[1])
        assert np.median((b[1] - b[0]) / (a[1] - a[0])) < 1.0
