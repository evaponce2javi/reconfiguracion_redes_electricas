"""Three-phase unbalanced networks in phase coordinates, as a :class:`Grid` whose "buses" are node-phases.

Model. Every bus has up to three phase conductors; the neutral is taken at ground potential at every bus
(multi-grounded neutral, Kron reduction), so the state is one phase-to-ground voltage per node-phase. A branch
with k phases is a k x k series impedance matrix (mutual coupling included) plus an optional k x k shunt
admittance, half at each end. The source is an ideal three-phase voltage source: the node-phases of the slack
bus are reference nodes with given phasors (1 at 0, -120 and +120 degrees times the source magnitude).

Loads are constant power and connected between a phase and the neutral (wye). The current drawn at node-phase
n is then conj(S_n / V_n), a function of the voltage of the same node-phase only, and the network equations are

    V_n conj((Y V)_n) = S_n      for every non-reference node-phase n,

with Y the admittance matrix over node-phases: the same equations as a single-phase network with PQ buses
only. Everything in this package that takes a :class:`Grid` with PQ buses therefore applies unchanged
(``half_disk_study``, ``SecondOrder``, ``voltage_enclosure``, ``disc_bound``, ``attained_extreme``).

Not modelled: loads between two phases (delta), whose current depends on two voltages; transformers; voltage
regulators; an explicit neutral conductor with its own voltage; the short-circuit impedance of the source.

Per unit: phase-to-neutral voltage vn / sqrt(3) and the per-phase power sn / 3, hence the impedance base
vn^2 / sn of the three-phase system. ``Grid.base_mva`` holds the per-phase base.
"""
from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
import scipy.sparse as sp

from .network import Grid
from .powerflow import newton_pf

PHASE_SHIFT = np.exp(-2j * np.pi * np.arange(3) / 3.0)        # phases a, b, c of a positive-sequence source


def sequence_to_phase(x0: complex, x1: complex) -> np.ndarray:
    """3 x 3 phase matrix of a transposed (balanced) element from its zero- and positive-sequence value.

    Works for impedances and for admittances: A diag(x0, x1, x1) A^-1 = m 11' + (s - m) I with
    s = (x0 + 2 x1) / 3 and m = (x0 - x1) / 3."""
    s, m = (x0 + 2.0 * x1) / 3.0, (x0 - x1) / 3.0
    return m * np.ones((3, 3), complex) + (s - m) * np.eye(3)


@dataclass
class Branch:
    """A line section between buses ``f`` and ``t`` on the phases ``phases`` (0 = a, 1 = b, 2 = c)."""

    f: int
    t: int
    phases: tuple
    z: np.ndarray                      # (k, k) series impedance [p.u.]
    y: np.ndarray | None = None        # (k, k) total shunt admittance [p.u.], half at each end; None = none


def merge_series(branches: list[Branch], n_bus: int, keep) -> tuple[list[Branch], np.ndarray]:
    """Remove the buses that only join two sections of the same phases without shunt admittance.

    Such a bus carries no injection (it must not be in ``keep``), so the two sections are one branch whose
    impedance matrix is the sum of theirs; the voltages of the remaining buses do not change. Returns the new
    branches and, for every old bus, its new index (-1 if removed)."""
    keep = set(int(b) for b in keep)
    live = {i: br for i, br in enumerate(branches)}
    at = [set() for _ in range(n_bus)]
    for i, br in live.items():
        at[br.f].add(i)
        at[br.t].add(i)
    nxt = len(branches)
    stack = [b for b in range(n_bus) if b not in keep and len(at[b]) == 2]
    gone = np.zeros(n_bus, bool)
    while stack:
        b = stack.pop()
        if gone[b] or b in keep or len(at[b]) != 2:
            continue
        i, j = sorted(at[b])
        p, q = live[i], live[j]
        u = p.t if p.f == b else p.f
        v = q.t if q.f == b else q.f
        if u == v or u == b or v == b or tuple(p.phases) != tuple(q.phases) or p.y is not None or q.y is not None:
            continue
        for e, end in ((i, u), (j, v)):
            del live[e]
            at[end].discard(e)
        at[b].clear()
        gone[b] = True
        live[nxt] = Branch(u, v, tuple(p.phases), p.z + q.z)
        at[u].add(nxt)
        at[v].add(nxt)
        nxt += 1
        stack += [w for w in (u, v) if w not in keep and len(at[w]) == 2]
    new = np.full(n_bus, -1)
    new[~gone] = np.arange(int((~gone).sum()))
    return [Branch(int(new[br.f]), int(new[br.t]), br.phases, br.z, br.y) for _, br in sorted(live.items())], new


def phase_ybus(n_bus: int, branches: list[Branch]):
    """Admittance matrix over node-phases. Returns (Y, node) with ``node[bus, phase]`` the index of that
    node-phase or -1 if no branch brings the phase to the bus."""
    has = np.zeros((n_bus, 3), bool)
    for br in branches:
        has[br.f, list(br.phases)] = True
        has[br.t, list(br.phases)] = True
    node = np.full((n_bus, 3), -1)
    node[has] = np.arange(int(has.sum()))
    rows, cols, vals = [], [], []
    for br in branches:
        ph = list(br.phases)
        z = np.atleast_2d(np.asarray(br.z, complex))
        if z.shape != (len(ph), len(ph)):
            raise ValueError("the impedance matrix of a branch does not match its phases")
        ys = np.linalg.inv(z)
        sh = np.zeros_like(ys) if br.y is None else 0.5 * np.atleast_2d(np.asarray(br.y, complex))
        nf, nt = node[br.f, ph], node[br.t, ph]
        for a, b, blk in ((nf, nf, ys + sh), (nt, nt, ys + sh), (nf, nt, -ys), (nt, nf, -ys)):
            rows.append(np.repeat(a, len(ph)))
            cols.append(np.tile(b, len(ph)))
            vals.append(blk.ravel())
    n = int(has.sum())
    Y = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    return Y, node


@dataclass
class ThreePhaseGrid(Grid):
    """A :class:`Grid` over node-phases. ``ref`` holds the node-phases of the source bus; there are no PV nodes."""

    node: np.ndarray | None = None        # (n_bus_3ph, 3): node-phase index of every bus and phase, -1 if absent
    bus_of: np.ndarray | None = None      # bus of every node-phase
    phase_of: np.ndarray | None = None    # phase of every node-phase (0 = a, 1 = b, 2 = c)
    v_ref: np.ndarray | None = None       # no-load profile used as flat start: the source phasor of each phase

    def flat_voltage(self) -> np.ndarray:
        return self.v_ref.copy()


def build_grid(name: str, n_bus: int, branches: list[Branch], slack_bus: int, v_slack, s_load: np.ndarray,
               bus_label=None, base_mva: float = 1.0, base_kv=None, source: str = "", tol: float = 1e-9) -> ThreePhaseGrid:
    """Grid over node-phases with an ideal three-phase source at ``slack_bus``.

    ``v_slack``: phasor of phase a of the source (b and c lag by 120 and 240 degrees), or the three phasors.
    ``s_load``: (n_bus, 3) complex demand per bus and phase [p.u. of the per-phase base], wye-connected.
    The base case is solved with the Newton solver of this package from the no-load profile."""
    Y, node = phase_ybus(n_bus, branches)
    v_src = np.asarray(v_slack, complex)
    v_src = v_src * PHASE_SHIFT if v_src.ndim == 0 else v_src
    s_load = np.asarray(s_load, complex)
    if np.any((np.abs(s_load) > 0) & (node < 0)):
        raise ValueError("a load is connected to a phase that does not reach its bus")
    n = Y.shape[0]
    bus_of, phase_of = np.nonzero(node >= 0)                  # row-major: the order in which nodes were numbered
    ref = node[slack_bus][node[slack_bus] >= 0]
    pq = np.setdiff1d(np.arange(n), ref)
    v_ref = v_src[phase_of]
    load = s_load[bus_of, phase_of]
    if np.any(np.abs(load[ref]) > 0):
        raise ValueError("a load at the source bus has no effect in this model; remove it")
    t0 = perf_counter()
    V0, ok, _, _ = newton_pf(Y, -load, v_ref, np.array([], dtype=int), pq, tol=tol, max_iter=50)
    t_pf = perf_counter() - t0
    if not ok:
        raise RuntimeError("the three-phase base case did not converge")
    f, t = [], []
    for br in branches:
        ph = list(br.phases)
        f.append(node[br.f, ph])
        t.append(node[br.t, ph])
    label = np.arange(n_bus) if bus_label is None else np.asarray(bus_label)
    kv = None if base_kv is None else np.asarray(base_kv, float)[bus_of]
    return ThreePhaseGrid(name, base_mva, Y, np.asarray(ref, dtype=int), np.array([], dtype=int), pq, V0, np.zeros(n, complex), load,
                          np.ones(n), np.concatenate(f), np.concatenate(t), label[bus_of], t_pf, source, kv,
                          node=node, bus_of=bus_of, phase_of=phase_of, v_ref=v_ref)


# ----------------------------------------------------------------------------------------------- pandapower
_UNSUPPORTED = ("trafo", "trafo3w", "impedance", "shunt", "gen", "ward", "xward", "dcline", "storage", "motor",
                "asymmetric_sgen", "sgen", "load", "switch")


def _value(row, name: str) -> float:
    """A numeric field of a table row; 0 if the column is missing or empty."""
    x = getattr(row, name, 0.0)
    return 0.0 if x is None or x != x else float(x)


def lv_side(net, vm_pu: float | None = None, s_sc_mva: float = 1e9):
    """Copy of a pandapower network with one transformer, cut at its low-voltage bus.

    The transformer, the external grid and every bus upstream of the transformer are removed and an external
    grid with (practically) no short-circuit impedance is placed at the low-voltage bus, with the voltage
    magnitude of the original source unless ``vm_pu`` is given. What is left has lines only."""
    import copy

    import pandapower as pp
    from pandapower.toolbox import drop_buses

    if len(net.trafo) != 1 or len(net.ext_grid) != 1:
        raise ValueError("expected exactly one transformer and one external grid")
    out = copy.deepcopy(net)
    lv, hv = int(out.trafo.lv_bus.iloc[0]), int(out.trafo.hv_bus.iloc[0])
    vm = float(out.ext_grid.vm_pu.iloc[0]) if vm_pu is None else vm_pu
    out.trafo.drop(out.trafo.index, inplace=True)
    out.ext_grid.drop(out.ext_grid.index, inplace=True)
    drop_buses(out, [hv])
    pp.create_ext_grid(out, lv, vm_pu=vm, va_degree=0.0, s_sc_max_mva=s_sc_mva, s_sc_min_mva=s_sc_mva,
                       rx_max=0.1, rx_min=0.1, x0x_max=1.0, r0x0_max=0.1)
    return out


def from_pandapower(net, merge: bool = False, name: str = "", tol: float = 1e-9) -> ThreePhaseGrid:
    """Phase-coordinate model of a pandapower network made of lines, wye-connected asymmetric loads and one
    external grid, i.e. what ``pandapower.runpp_3ph`` solves for such a network with an ideal source.

    Lines are given by their positive- and zero-sequence parameters (transposed lines): the 3 x 3 phase blocks
    follow from ``sequence_to_phase``. ``merge``: remove the buses that only join two line sections
    (``merge_series``; exact). Raises for anything this model does not cover, delta loads in particular."""
    for el in _UNSUPPORTED:
        if el in net and len(net[el]):
            raise ValueError(f"element '{el}' is not covered by the phase-coordinate model")
    if len(net.ext_grid) != 1:
        raise ValueError("expected exactly one external grid")
    al = net.asymmetric_load[net.asymmetric_load.in_service]
    if (al.type.str.lower() != "wye").any():
        raise ValueError("delta-connected loads are outside this model (their current depends on two phase voltages)")
    idx = net.bus.index.values
    pos = {int(b): i for i, b in enumerate(idx)}
    sn = float(net.sn_mva)
    omega = 2.0 * np.pi * float(net.f_hz)
    branches = []
    for ln in net.line[net.line.in_service].itertuples():
        f, t = pos[int(ln.from_bus)], pos[int(ln.to_bus)]
        zb = float(net.bus.vn_kv.at[ln.from_bus]) ** 2 / sn
        par = float(ln.parallel)
        z1 = complex(ln.r_ohm_per_km, ln.x_ohm_per_km) * ln.length_km / par / zb
        z0 = complex(ln.r0_ohm_per_km, ln.x0_ohm_per_km) * ln.length_km / par / zb
        y1 = complex(_value(ln, "g_us_per_km") * 1e-6, omega * _value(ln, "c_nf_per_km") * 1e-9) * ln.length_km * par * zb
        y0 = complex(_value(ln, "g0_us_per_km") * 1e-6, omega * _value(ln, "c0_nf_per_km") * 1e-9) * ln.length_km * par * zb
        y = sequence_to_phase(y0, y1) if (y0 != 0 or y1 != 0) else None
        branches.append(Branch(f, t, (0, 1, 2), sequence_to_phase(z0, z1), y))
    s_load = np.zeros((len(idx), 3), complex)
    for ld in al.itertuples():
        s_load[pos[int(ld.bus)]] += ld.scaling * np.array([complex(ld.p_a_mw, ld.q_a_mvar), complex(ld.p_b_mw, ld.q_b_mvar),
                                                            complex(ld.p_c_mw, ld.q_c_mvar)]) / (sn / 3.0)
    eg = net.ext_grid.iloc[0]
    slack = pos[int(eg.bus)]
    label, n_bus = idx, len(idx)
    if merge:
        keep = set(np.flatnonzero(np.abs(s_load).sum(axis=1) > 0)) | {slack}
        branches, new = merge_series(branches, n_bus, keep)
        kept = np.flatnonzero(new >= 0)
        s_load, label, slack, n_bus = s_load[kept], idx[kept], int(new[slack]), len(kept)
    v_a = float(eg.vm_pu) * np.exp(1j * np.deg2rad(float(eg.va_degree)))
    kv = net.bus.vn_kv.loc[label].values / np.sqrt(3.0)
    return build_grid(name or "pandapower_3ph", n_bus, branches, slack, v_a, s_load, bus_label=label, base_mva=sn / 3.0,
                      base_kv=kv, source="phase-coordinate model built from the line data of a pandapower network"
                      + (" (series sections merged)" if merge else ""), tol=tol)


def pandapower_phase_voltages(net) -> np.ndarray:
    """(n_bus, 3) complex phase voltages [p.u.] of the last ``pandapower.runpp_3ph`` run, in the order of net.bus."""
    r = net.res_bus_3ph.loc[net.bus.index]
    return np.stack([r[f"vm_{p}_pu"].values * np.exp(1j * np.deg2rad(r[f"va_{p}_degree"].values)) for p in "abc"], axis=1)
