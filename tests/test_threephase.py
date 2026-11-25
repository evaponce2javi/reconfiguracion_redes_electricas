"""Three-phase networks in phase coordinates: the model against pandapower, then the enclosure on node-phases."""
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

from cdhim import SecondOrder, custom_study
from cdhim.discs import disc_bound
from cdhim.enclosure import InjectionEnclosure, injection_changes, voltage_enclosure
from cdhim.halfdisk import attained_extreme, first_order_range, half_disk_study, sample_half_disks
from cdhim.powerflow import newton_pf
from cdhim.threephase import (PHASE_SHIFT, Branch, build_grid, from_pandapower, merge_series, pandapower_phase_voltages,
                              sequence_to_phase)

# An untransposed line: unequal self and mutual impedances between the phases [p.u. per unit length]
Z_LINE = 0.15 * np.array([[0.3465 + 1.0179j, 0.1560 + 0.5017j, 0.1580 + 0.4236j],
                          [0.1560 + 0.5017j, 0.3375 + 1.0478j, 0.1535 + 0.3849j],
                          [0.1580 + 0.4236j, 0.1535 + 0.3849j, 0.3414 + 1.0348j]])
Y_LINE = 0.02j * np.array([[1.0, -0.2, -0.1], [-0.2, 0.9, -0.15], [-0.1, -0.15, 0.95]])     # shunt admittance of a section


def sub(M, phases):
    return M[np.ix_(phases, phases)]


def unbalanced_grid():
    """Source, a three-phase trunk with a loop, a two-phase lateral (a, c) and a single-phase lateral (b)."""
    br = [Branch(0, 1, (0, 1, 2), 1.0 * Z_LINE, 1.0 * Y_LINE), Branch(1, 2, (0, 1, 2), 0.8 * Z_LINE, 0.8 * Y_LINE),
          Branch(2, 5, (0, 1, 2), 0.7 * Z_LINE), Branch(5, 1, (0, 1, 2), 0.9 * Z_LINE),
          Branch(2, 4, (0, 2), 0.6 * sub(Z_LINE, [0, 2])), Branch(1, 3, (1,), 0.6 * sub(Z_LINE, [1]))]
    load = np.zeros((6, 3), complex)
    load[1] = [0.05 + 0.02j, 0.03 + 0.01j, 0.04 + 0.015j]
    load[2] = [0.06 + 0.02j, 0.0, 0.08 + 0.03j]
    load[3, 1] = 0.07 + 0.02j
    load[4] = [0.03 + 0.01j, 0.0, 0.05 + 0.02j]
    load[5] = [0.0, 0.09 + 0.03j, 0.02 + 0.01j]
    return br, build_grid("unbalanced", 6, br, 0, 1.02, load, tol=1e-12)


def half_disk_case(g, sources, cap_pu):
    """Nominal injections (the loads) and the half-disk capacity of single-phase sources at (bus, phase) pairs."""
    cap = np.zeros(g.n_bus)
    cap[[g.node[b, p] for b, p in sources]] = cap_pu
    return -g.s_load.copy(), cap


@pytest.fixture(scope="module")
def model():
    br, g = unbalanced_grid()
    s0, cap = half_disk_case(g, [(3, 1), (4, 2), (5, 0)], 0.06)
    st = half_disk_study(g, s0, cap, tol=1e-12)
    return br, g, s0, cap, st, SecondOrder(st, 1.0)


def vertices(k):
    return np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(2 ** k)] for j in range(k)])


def solve(st, eps):
    V, conv, _ = st.solver.solve(st.injections(eps, 1.0))
    assert conv.all()
    return V


# ------------------------------------------------------------------------------------------------ the model
def test_sequence_to_phase_has_the_sequence_values_as_eigenvalues():
    z0, z1 = 1.2 + 0.3j, 0.32 + 0.08j
    a = np.exp(2j * np.pi / 3)
    A = np.array([[1, 1, 1], [1, a * a, a], [1, a, a * a]])
    np.testing.assert_allclose(np.linalg.solve(A, sequence_to_phase(z0, z1) @ A), np.diag([z0, z1, z1]), atol=1e-14)


def test_node_phases_and_kirchhoff_laws_branch_by_branch(model):
    """The solved voltages satisfy the current balance computed from the branch impedances, not from Y."""
    br, g = model[0], model[1]
    assert g.n_bus == 3 * 4 + 2 + 1 and len(g.ref) == 3 and len(g.pv) == 0          # missing phases have no node
    assert g.node[3].tolist() == [-1, g.node[3, 1], -1] and g.node[4, 1] == -1
    np.testing.assert_allclose(g.V0[g.ref], 1.02 * PHASE_SHIFT, atol=0)
    out = np.zeros(g.n_bus, complex)                                                # current leaving every node into the lines
    for b in br:
        ph = list(b.phases)
        vf, vt = g.V0[g.node[b.f, ph]], g.V0[g.node[b.t, ph]]
        i = np.linalg.solve(np.atleast_2d(b.z), vf - vt)
        sh = np.zeros((len(ph), len(ph))) if b.y is None else 0.5 * np.atleast_2d(b.y)
        out[g.node[b.f, ph]] += i + sh @ vf
        out[g.node[b.t, ph]] += -i + sh @ vt
    drawn = np.conj(g.s_load / g.V0)                                                # wye constant-power load: conj(S / V)
    np.testing.assert_allclose(out[g.pq], -drawn[g.pq], atol=1e-11)
    vm = np.abs(g.V0)
    spread = max(np.ptp(vm[g.node[b][g.node[b] >= 0]]) for b in range(1, 6))
    assert spread > 0.01                                                            # the case is unbalanced


def _pp_net(kind="wye"):
    import pandapower as pp

    net = pp.create_empty_network(sn_mva=10.0)                     # 20 kV cables: the line charging is not negligible
    b = [pp.create_bus(net, 20.0) for _ in range(6)]
    pp.create_ext_grid(net, b[0], vm_pu=1.02, s_sc_max_mva=1e9, s_sc_min_mva=1e9, rx_max=0.1, rx_min=0.1, x0x_max=1.0, r0x0_max=0.1)
    for f, t, km in ((0, 1, 8.0), (1, 2, 6.0), (2, 3, 5.0), (1, 4, 7.0), (4, 5, 9.0), (3, 5, 6.0)):   # one loop
        pp.create_line_from_parameters(net, b[f], b[t], km, r_ohm_per_km=0.20, x_ohm_per_km=0.12, c_nf_per_km=300.0, max_i_ka=0.4,
                                       r0_ohm_per_km=0.80, x0_ohm_per_km=0.45, c0_nf_per_km=180.0)
    for bus, mw in ((2, (1.2, 0, 0.4)), (3, (0, 1.8, 0)), (4, (0.6, 0.6, 1.6)), (5, (0, 0, 1.4))):
        pp.create_asymmetric_load(net, b[bus], p_a_mw=mw[0], q_a_mvar=0.3 * mw[0], p_b_mw=mw[1], q_b_mvar=0.3 * mw[1],
                                  p_c_mw=mw[2], q_c_mvar=0.3 * mw[2], type=kind)
    return net


def test_matches_the_three_phase_power_flow_of_pandapower():
    import pandapower as pp

    net = _pp_net()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp_3ph(net)
    assert net.converged
    g = from_pandapower(net, tol=1e-12)
    v_pp = pandapower_phase_voltages(net)[g.bus_of, g.phase_of]
    vm = np.abs(v_pp).reshape(-1, 3)
    assert np.ptp(vm, axis=1).max() > 0.01                                          # unbalanced: > 1 % between phases
    assert abs(np.asarray(g.Ybus.sum(axis=1))).max() > 0.01                         # line charging [p.u.]: it matters here
    # the two solutions agree to the accuracy at which pandapower stops (a power mismatch of 3e-8 MVA)
    assert np.abs(g.V0 - v_pp).max() < 1e-7
    residual = (v_pp * np.conj(g.Ybus @ v_pp) + g.s_load)[g.pq] * g.base_mva        # of pandapower's solution, in this model
    assert np.abs(residual).max() < 1e-6                                            # [MVA]


def test_delta_loads_and_transformers_are_rejected():
    import pandapower as pp

    with pytest.raises(ValueError, match="delta"):
        from_pandapower(_pp_net("delta"))
    net = _pp_net()
    hv = pp.create_bus(net, 110.0)
    pp.create_transformer_from_parameters(net, hv, net.bus.index[1], sn_mva=25.0, vn_hv_kv=110.0, vn_lv_kv=20.0, vk_percent=12.0,
                                          vkr_percent=0.4, pfe_kw=0.0, i0_percent=0.0)
    with pytest.raises(ValueError, match="trafo"):
        from_pandapower(net)


def test_merging_series_sections_does_not_change_the_voltages():
    z = sequence_to_phase(0.9 + 0.3j, 0.3 + 0.1j) * 0.02
    br = [Branch(i, i + 1, (0, 1, 2), (1.0 + 0.1 * i) * z) for i in range(8)] + [Branch(4, 9, (0, 1, 2), z), Branch(9, 10, (0, 1, 2), 2 * z)]
    load = np.zeros((11, 3), complex)
    load[8], load[10], load[4] = [0.3 + 0.1j, 0.0, 0.1], [0.0, 0.25 + 0.1j, 0.0], [0.0, 0.0, 0.2 + 0.05j]
    full = build_grid("chain", 11, br, 0, 1.0, load, tol=1e-12)
    merged_br, new = merge_series(br, 11, keep={0, 4, 8, 10})
    kept = np.flatnonzero(new >= 0)
    assert kept.tolist() == [0, 4, 8, 10] and len(merged_br) == 3
    small = build_grid("chain", len(kept), merged_br, 0, 1.0, load[kept], tol=1e-12)
    np.testing.assert_allclose(small.V0, full.V0[full.node[kept].ravel()], atol=1e-11)


# ------------------------------------------------------------------------------------------- the enclosure
def test_exact_identity_of_the_current_injection_form_on_node_phases(model):
    """J D = dS + V*^2 (-tau d + sig d^2) / (1 + d): the equations have the single-phase form on every node-phase."""
    br, g, s0, cap, st, so = model
    ns, m = g.pvpq, so.m
    eps = np.random.default_rng(3).uniform(-1, 1, (st.src.n, 20))
    Vs = st.V_star[ns][:, None]
    D = solve(st, eps)[ns] - Vs
    d = D / Vs
    dS = injection_changes(so)[ns] @ eps
    rhs = dS + Vs ** 2 * (-(dS / Vs ** 2) * d + (st.sbus[ns][:, None] / Vs ** 2) * d * d) / (1 + d)
    x = so.lu.solve(np.vstack([rhs.real, rhs.imag]))
    np.testing.assert_allclose(x[:m] + 1j * x[m:], D, rtol=0, atol=1e-10)


def test_partial_bounds_and_remainder_hold_at_every_vertex(model):
    br, g, s0, cap, st, so = model
    enc = InjectionEnclosure(so)
    fp = enc.solve(1.0)
    assert fp.converged and fp.rel_radius < 1.0
    ns, m, k = g.pvpq, so.m, st.src.n
    assert k == 6
    eps = vertices(k)
    Vs = st.V_star[ns][:, None]
    D = solve(st, eps)[ns] - Vs
    D1 = so.E[ns] @ eps
    r2 = enc.Hs.reshape(2 * m, -1) @ (eps[:, None, :] * eps[None, :, :]).reshape(k * k, -1)
    D2 = r2[:m] + 1j * r2[m:]
    tau = injection_changes(so)[ns] @ eps / Vs ** 2
    tol = 1e-12
    assert np.all(np.abs(D1 / Vs).max(axis=1) <= enc.a1 + tol)
    assert np.all(np.abs(D2 / Vs).max(axis=1) <= enc.a2 + tol)
    assert np.all(np.abs(tau).max(axis=1) <= enc.t + tol)
    s_true = np.abs(D - D1 - D2).max(axis=1)
    assert np.all(s_true <= fp.s_abs[ns] + 1e-11)
    assert s_true.max() > 0.05 * fp.s_abs[ns].max()                                # the bound is not vacuous


def test_enclosure_contains_the_ac_voltages_of_every_node_phase(model):
    br, g, s0, cap, st, so = model
    enc = voltage_enclosure(so, first_order_range(so.c_lin))
    box = voltage_enclosure(so)
    assert enc.converged and box.converged
    dg = np.flatnonzero(cap > 0)
    tol = 1e-10
    vm = np.abs(solve(st, vertices(st.src.n))[g.pq])                               # the box, corners included
    assert np.all(vm.min(axis=1) >= box.lo - tol) and np.all(vm.max(axis=1) <= box.hi + tol)
    pts = sample_half_disks(np.random.default_rng(6), len(dg), 3000)              # the half disks
    S = np.repeat(s0[:, None], 3000, axis=1)
    S[dg] += cap[dg][:, None] * pts
    V, conv, _ = st.solver.solve(S)
    assert conv.all()
    vm = np.abs(V[g.pq])
    assert np.all(vm.min(axis=1) >= enc.lo - tol) and np.all(vm.max(axis=1) <= enc.hi + tol)
    ext = {sign: attained_extreme(g, s0, cap, sign, tol=1e-12) for sign in (1, -1)}
    assert np.all(enc.hi >= ext[1] - tol) and np.all(enc.lo <= ext[-1] + tol)      # and the extremes found by local search
    width = ext[1] - ext[-1]
    assert width.min() > 1e-4                                                      # every node-phase moves: the phases are coupled
    assert np.median((enc.hi - enc.lo) / width) < 1.1                              # within 10 % of the attained range
    # a source on one phase moves the other phases the other way: the minimum is not at (dp, dq) = (0, -C) everywhere
    corner, ok, _, _ = newton_pf(g.Ybus, s0 - 1j * cap, g.V0, g.pv, g.pq, tol=1e-12)
    assert ok and (np.abs(corner[g.pq]) > ext[-1] + 1e-6).any()


def test_disc_bound_of_wang_et_al_holds_and_is_wider(model):
    br, g, s0, cap, st, so = model
    v_hat, ok, _, _ = newton_pf(g.Ybus, s0, g.V0, g.pv, g.pq, tol=1e-12)
    db = disc_bound(g, s0, v_hat, cap)
    assert ok and db.ok
    dg = np.flatnonzero(cap > 0)                                                   # the discs cover |ds| <= C: the half disks,
    S = np.repeat(s0[:, None], 3000, axis=1)                                       # not the corners of their bounding boxes
    S[dg] += cap[dg][:, None] * sample_half_disks(np.random.default_rng(7), len(dg), 3000)
    V, conv, _ = st.solver.solve(S)
    assert conv.all()
    V = V[g.pq]
    assert np.all(np.abs(V - v_hat[g.pq][:, None]).max(axis=1) <= db.radius_bus + 1e-10)
    enc = voltage_enclosure(so, first_order_range(so.c_lin))
    assert np.median(2 * db.radius_bus / (enc.hi - enc.lo)) > 1.2


# ------------------------------------------------------------------------- balanced cases: the single-phase limit
def _balanced(z0, z1):
    """A balanced three-phase feeder with a loop and the single-phase network of its positive sequence."""
    topo = ((0, 1, 1.0), (1, 2, 0.8), (2, 3, 0.6), (1, 3, 0.9))
    per_phase = {1: 0.10 + 0.04j, 2: 0.15 + 0.05j, 3: 0.12 + 0.03j}
    load3, load1 = np.zeros((4, 3), complex), np.zeros((4, 3), complex)
    for b, s in per_phase.items():
        load3[b], load1[b, 0] = s, s
    g3 = build_grid("balanced", 4, [Branch(f, t, (0, 1, 2), km * sequence_to_phase(z0, z1)) for f, t, km in topo], 0, 1.03, load3, tol=1e-12)
    g1 = build_grid("positive sequence", 4, [Branch(f, t, (0,), np.array([[km * z1]])) for f, t, km in topo], 0, 1.03, load1, tol=1e-12)
    return g3, g1


def test_balanced_power_flow_is_the_positive_sequence_network():
    g3, g1 = _balanced(0.20 + 0.50j, 0.05 + 0.08j)
    for p in range(3):
        np.testing.assert_allclose(g3.V0[g3.node[:, p]], g1.V0 * PHASE_SHIFT[p], atol=1e-11)


def test_without_coupling_between_phases_the_enclosure_is_the_single_phase_one():
    z1 = 0.05 + 0.08j
    g3, g1 = _balanced(z1, z1)                                                     # z0 = z1: no mutual impedance
    s3, cap3 = half_disk_case(g3, [(2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)], 0.08)
    s1, cap1 = half_disk_case(g1, [(2, 0), (3, 0)], 0.08)
    so3, so1 = SecondOrder(half_disk_study(g3, s3, cap3, tol=1e-12), 1.0), SecondOrder(half_disk_study(g1, s1, cap1, tol=1e-12), 1.0)
    e1 = voltage_enclosure(so1, first_order_range(so1.c_lin))
    e3 = voltage_enclosure(so3, first_order_range(so3.c_lin))
    assert e1.converged and e3.converged
    pos = np.full(g3.n_bus, -1)
    pos[g3.pq] = np.arange(len(g3.pq))
    for p in range(3):
        o = pos[g3.node[1:, p]]
        np.testing.assert_allclose(e3.lo2[o], e1.lo2, rtol=0, atol=1e-11)          # first and second order: identical
        np.testing.assert_allclose(e3.hi2[o], e1.hi2, rtol=0, atol=1e-11)
        # The bound of the remainder is identical on phase a only. InjectionEnclosure bounds the second-order state
        # term through its real and imaginary parts, which depends on the angle of the phasor: on phases b and c
        # (-120 and +120 degrees) the same term gets a slightly different, equally valid, bound.
        np.testing.assert_allclose(e3.allowance[o], e1.allowance, rtol=1e-11 if p == 0 else 1e-2, atol=0)
    np.testing.assert_allclose(e3.lo[pos[g3.node[1:, 0]]], e1.lo, rtol=0, atol=1e-11)
    np.testing.assert_allclose(e3.hi[pos[g3.node[1:, 0]]], e1.hi, rtol=0, atol=1e-11)


def test_balanced_source_on_coupled_phases_second_order_equal_remainder_bound_not_smaller():
    """With mutual impedances and a balanced three-phase source (one pair of symbols for its three phases) the
    voltages stay balanced, so the first- and second-order bounds are those of the positive-sequence network.
    The bound of the remainder adds the three phases with absolute values and cannot be smaller."""
    g3, g1 = _balanced(0.20 + 0.50j, 0.05 + 0.08j)
    cap = 0.08
    nodes = g3.node[2]
    centre = -g3.s_load.copy()
    centre[nodes] += cap / 2.0
    Z = sp.csr_matrix((np.ones(6), (np.arange(6), np.repeat([0, 1], 3))), shape=(6, 2))
    st3 = custom_study(g3, centre, np.r_[nodes, nodes], np.r_[np.zeros(3, bool), np.ones(3, bool)], np.r_[np.full(3, cap / 2.0), np.full(3, cap)],
                       Z=Z, tol=1e-12)
    s1, cap1 = half_disk_case(g1, [(2, 0)], cap)
    so3, so1 = SecondOrder(st3, 1.0), SecondOrder(half_disk_study(g1, s1, cap1, tol=1e-12), 1.0)
    e1 = voltage_enclosure(so1, first_order_range(so1.c_lin))
    e3 = voltage_enclosure(so3, first_order_range(so3.c_lin))
    assert e1.converged and e3.converged
    pos = np.full(g3.n_bus, -1)
    pos[g3.pq] = np.arange(len(g3.pq))
    for p in range(3):
        o = pos[g3.node[1:, p]]
        np.testing.assert_allclose(e3.lo2[o], e1.lo2, rtol=0, atol=1e-11)
        np.testing.assert_allclose(e3.hi2[o], e1.hi2, rtol=0, atol=1e-11)
        assert np.all(e3.allowance[o] >= e1.allowance * (1.0 - 1e-9))
    assert np.all(e3.allowance[pos[g3.node[1:, 0]]] > 1.5 * e1.allowance)           # and here it is clearly larger
    pts = sample_half_disks(np.random.default_rng(8), 1, 2000)                     # (u, w) of the half disk -> the two symbols
    V = np.abs(solve(st3, np.vstack([2.0 * pts.real - 1.0, pts.imag]))[g3.pq])
    assert np.all(V.min(axis=1) >= e3.lo - 1e-10) and np.all(V.max(axis=1) <= e3.hi + 1e-10)
