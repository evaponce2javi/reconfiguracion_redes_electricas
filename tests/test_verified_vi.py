"""The verified tube in voltage-current form (PV buses allowed): the defect in exact rational arithmetic, the set
against AC solutions and against the operator itself, and the enclosure against the current form on PQ networks."""
import sys
import warnings
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import pytest

from cdhim import SecondOrder, load_grid, make_study, spectral_zones
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import Enclosure, injection_changes, state_forms
from cdhim.halfdisk import half_disk_study
from cdhim.study import custom_study
from cdhim.verified import verified_enclosure
from cdhim.verified_vi import block_frobenius_up, block_spectral_up, verified_enclosure_vi, verify_tube_vi
from cdhim.zones import laplacian_embedding

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))


def _close_ties(net):
    net.line["in_service"] = True


def _grid(case, prepare=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid(case, prepare=prepare)


def _data(so):
    g, st = so.st.grid, so.st
    ns, m, npv = g.pvpq, so.m, len(g.pv)
    Hs = state_forms(so)
    return dict(ns=ns, is_pv=np.r_[np.ones(npv, bool), np.zeros(m - npv, bool)], v_set=np.asarray(g.v_set, dtype=float)[ns], V=st.V_star,
                S0=st.sbus[ns], M=injection_changes(so)[ns], E=so.E[ns], B=Hs[:m] + 1j * Hs[m:])


def _pv_feeder():
    """The 33-bus feeder with ties and its X8 sources as voltage-controlled generators (set-up of experiments/x14_pv.py)."""
    from x14_pv import pv_variant

    g = _grid("case33bw", _close_ties)
    sc = SCENARIOS["R"]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], 0.5)
    gpv, s_c, dg = pv_variant(g, s0, cap)
    st = custom_study(gpv, s_c, dg, np.zeros(len(dg), bool), cap[dg] / 2.0)
    return st, 1.0


def _transmission():
    g = _grid("case30")
    return make_study(g, zones=spectral_zones(g, 4, laplacian_embedding(g))), 0.05


def _pq_feeder():
    g = _grid("case33bw", _close_ties)
    sc = SCENARIOS["R"]
    s0, cap = net_injections(g, benchmark_ratings(g), sc["rho_D"], sc["rho_G"], sc["pf"], 0.5)
    return half_disk_study(g, s0, cap), 1.0


@pytest.fixture(scope="module", params=[_pv_feeder, _transmission, _pq_feeder], ids=["feeder with PV sources", "case30", "PQ feeder"])
def model(request):
    st, delta = request.param()
    so = SecondOrder(st, delta)
    return st, delta, so, _data(so)


def vertices(k):
    return np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** k)] for j in range(k)])


class Q2:
    """Complex numbers with rational parts."""

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


def test_block_norms_are_upper_bounds():
    rng = np.random.default_rng(0)
    m = 6
    A = rng.standard_normal((2 * m, 2 * m)) * 10.0 ** rng.uniform(-6, 3, (2 * m, 2 * m))
    sv, fr = block_spectral_up(A, m), block_frobenius_up(np.abs(A), m)
    for i in range(m):
        for j in range(m):
            blk = np.array([[A[i, j], A[i, m + j]], [A[m + i, j], A[m + i, m + j]]])
            n2 = np.linalg.norm(blk, 2)
            assert n2 <= sv[i, j] <= n2 * (1 + 1e-12) + 1e-300 and sv[i, j] <= fr[i, j] * (1 + 1e-12)


def test_defect_bound_holds_in_exact_arithmetic(model):
    """phi_j >= 2-norm of the defect of bus j at vertices of the box, evaluated in rational arithmetic: the complex
    equation of a PQ bus, and the pair (active power, squared magnitude) of a PV bus."""
    st, delta, so, d = model
    tube = verify_tube_vi(st.grid.Ybus, keep=True, **d)
    assert tube.exists
    g = st.grid
    ns, (m, K) = d["ns"], d["E"].shape
    Y = g.Ybus.tocsr()
    pos = np.full(g.n_bus, -1)
    pos[ns] = np.arange(m)
    Vq = [Q2.of(v) for v in d["V"]]
    rng = np.random.default_rng(1)
    buses = np.unique(np.r_[np.flatnonzero(d["is_pv"]), rng.choice(m, 6, replace=False)])
    for _ in range(2):
        eps = [int(e) for e in rng.choice([-1, 1], K)]
        mod = []
        for j in range(m):
            v = Q2.of(d["V"][ns[j]])
            for k in range(K):
                v = v + Q2.of(d["E"][j, k]) * Q2(eps[k])
                for l in range(K):
                    v = v + Q2.of(d["B"][j, k, l]) * Q2(eps[k] * eps[l])
            mod.append(v)
        for j in buses:
            cur = Q2(0)
            for idx in range(Y.indptr[ns[j]], Y.indptr[ns[j] + 1]):
                c = Y.indices[idx]
                cur = cur + Q2.of(Y.data[idx]) * (mod[pos[c]] if pos[c] >= 0 else Vq[c])
            S = Q2.of(d["S0"][j])
            for k in range(K):
                S = S + Q2.of(d["M"][j, k]) * Q2(eps[k])
            w = mod[j].conj() * cur
            if d["is_pv"][j]:
                n1 = w.re - S.re
                n2 = mod[j].re ** 2 + mod[j].im ** 2 - Fr(float(d["v_set"][j])) ** 2
            else:
                n1, n2 = w.re - S.re, w.im + S.im                      # conj(m) I - conj(S)
            assert n1 ** 2 + n2 ** 2 <= Fr(float(tube.aux["phi"][j])) ** 2


def test_set_contains_the_ac_solution(model):
    st, delta, so, d = model
    g = st.grid
    tube = verify_tube_vi(g.Ybus, **d)
    assert tube.exists and tube.unique
    ns, K = d["ns"], d["E"].shape[1]
    eps = vertices(K)
    V, conv, _ = st.solver.solve(st.injections(eps, delta))
    assert conv.all()
    s = V[ns] - (d["V"][ns][:, None] + d["E"] @ eps + np.einsum("jkl,kv,lv->jv", d["B"], eps, eps))
    YL = g.Ybus.tocsr()[ns][:, ns]
    tol = 1e-9 * (1.0 + abs(YL).sum(axis=1).max())                    # tolerance of the AC solver, seen through Y
    assert np.all(np.abs(s) <= tube.z[:, None] + 1e-9)
    assert np.all(np.abs(YL @ s) <= tube.zeta[:, None] + tol)


def test_operator_maps_the_set_into_itself(model):
    """Psi(s) = s - C N(m + s) at sampled eps and s of the set D(z, zeta) stays in the set."""
    st, delta, so, d = model
    g = st.grid
    tube = verify_tube_vi(g.Ybus, keep=True, **d)
    ns, (m, K), C = d["ns"], d["E"].shape, tube.aux["C"]
    Y = g.Ybus.tocsr()
    YL = Y[ns][:, ns]
    lu = __import__("scipy.sparse.linalg", fromlist=["splu"]).splu(YL.tocsc())
    S0 = np.where(d["is_pv"], d["S0"].real, d["S0"])
    M = np.where(d["is_pv"][:, None], d["M"].real, d["M"])
    rng = np.random.default_rng(2)
    worst_v = worst_i = 0.0
    for trial in range(100):
        eps = rng.choice([-1.0, 1.0], K) if trial % 2 else rng.uniform(-1, 1, K)
        s = lu.solve(tube.zeta * np.exp(1j * rng.uniform(0, 2 * np.pi, m)))          # |Y s| = zeta
        s = s * min(1.0, float((tube.z / np.maximum(np.abs(s), 1e-300)).min()))       # ... pulled back until |s| <= z too
        Vf = d["V"].copy()
        Vf[ns] = d["V"][ns] + d["E"] @ eps + np.einsum("jkl,k,l->j", d["B"], eps, eps) + s
        w = np.conj(Vf[ns]) * (Y @ Vf)[ns] - np.conj(S0 + M @ eps)
        n2 = np.where(d["is_pv"], np.abs(Vf[ns]) ** 2 - d["v_set"] ** 2, w.imag)
        CN = C @ np.r_[w.real, n2]
        psi = s - (CN[:m] + 1j * CN[m:])
        worst_v = max(worst_v, float((np.abs(psi) / tube.z).max()))
        worst_i = max(worst_i, float((np.abs(YL @ psi) / tube.zeta).max()))
    assert worst_v <= 1.0 + 1e-6 and worst_i <= 1.0 + 1e-6


def test_enclosure_contains_ac_solutions(model):
    st, delta, so, d = model
    g = st.grid
    enc = verified_enclosure_vi(so)
    assert enc.verified and enc.unique
    K = d["E"].shape[1]
    eps = np.hstack([vertices(K), np.random.default_rng(3).uniform(-1, 1, (K, 300))])
    V, conv, _ = st.solver.solve(st.injections(eps, delta))
    vm = np.abs(V[g.pq][:, conv])
    assert conv.sum() >= 2 ** K
    assert np.all(vm >= enc.lo[:, None] - 1e-9) and np.all(vm <= enc.hi[:, None] + 1e-9)
    if len(g.pv):                                                      # the controlled magnitudes stay at their set-points
        assert np.abs(np.abs(V[g.pv][:, conv]) - np.asarray(g.v_set)[g.pv][:, None]).max() < 1e-8


def test_agrees_with_the_current_form_on_a_pq_network():
    st, delta = _pq_feeder()
    so = SecondOrder(st, delta)
    a, b = verified_enclosure(so, half_disks=True), verified_enclosure_vi(so, half_disks=True)
    assert a.verified and b.verified and a.unique and b.unique
    r = b.s_abs / a.s_abs
    assert 0.9 < r.min() and r.max() < 1.1
    assert np.abs((b.hi - b.lo) / (a.hi - a.lo) - 1.0).max() < 1e-3


def test_diagnostics_of_the_experiment_match_the_module(model):
    """experiments/x14_pv.py recomputes the matrices of the radii map and the linear radius of Enclosure."""
    from x14_pv import linear_radius, majorant_matrices

    st, delta, so, d = model
    tube = verify_tube_vi(st.grid.Ybus, keep=True, **d)
    q = majorant_matrices(so)
    for key in ("RV", "TV", "RI", "TI"):
        assert np.array_equal(q[key], tube.aux[key])
    assert linear_radius(so) == pytest.approx(Enclosure(so).spectral_radius, rel=1e-12)
