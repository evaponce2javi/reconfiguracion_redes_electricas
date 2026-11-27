"""Test systems: load a case with pandapower once, keep only sparse arrays afterwards."""
from __future__ import annotations

import hashlib
import re
import urllib.request
import warnings
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy.sparse as sp

from .powerflow import mismatch, newton_pf, stack_mismatch

PGLIB_URL = "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/master/{name}.m"
MATPOWER_URL = "https://raw.githubusercontent.com/MATPOWER/matpower/7.1/data/{name}.m"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@dataclass
class Grid:
    """Bus-branch model in per unit, in the internal (consecutive) bus numbering of pandapower."""

    name: str
    base_mva: float
    Ybus: sp.csr_matrix
    ref: np.ndarray
    pv: np.ndarray
    pq: np.ndarray
    V0: np.ndarray            # converged complex bus voltages of the base case (load factor 1)
    s_gen: np.ndarray         # complex generation per bus [p.u.]; only Re at PV buses is a set-point
    s_load: np.ndarray        # complex net demand per bus [p.u.] (PD + jQD; negative PD = injection)
    v_set: np.ndarray         # |V| set-point at PV/slack buses, 1.0 elsewhere (flat-start profile)
    br_f: np.ndarray          # from/to bus of every in-service branch
    br_t: np.ndarray
    bus_label: np.ndarray     # pandapower bus index of each internal bus
    base_pf_seconds: float    # wall time of one (warm) pandapower base power flow
    source: str = ""          # where the data came from (and a hash for downloaded files)
    base_kv: np.ndarray | None = None     # nominal voltage of every bus [kV]
    Yf: sp.csr_matrix | None = None       # branch current at the from / to end of every in-service branch: I = Yf V
    Yt: sp.csr_matrix | None = None

    @property
    def n_bus(self) -> int:
        return self.Ybus.shape[0]

    @property
    def pvpq(self) -> np.ndarray:
        return np.concatenate([self.pv, self.pq])

    def sbus(self, load_factor: float = 1.0) -> np.ndarray:
        """Specified injections when demand and scheduled generation are scaled by ``load_factor``."""
        return load_factor * self.s_gen.real + 1j * self.s_gen.imag - load_factor * self.s_load

    def flat_voltage(self) -> np.ndarray:
        return self.v_set.astype(complex)

    def pq_pq_branches(self):
        """Positions (within ``pq``) of both ends of every branch whose two terminals are PQ buses."""
        pos = np.full(self.n_bus, -1)
        pos[self.pq] = np.arange(len(self.pq))
        keep = (pos[self.br_f] >= 0) & (pos[self.br_t] >= 0) & (self.br_f != self.br_t)
        return pos[self.br_f[keep]], pos[self.br_t[keep]], np.flatnonzero(keep)


def fetch_pglib(name: str, data_dir: Path = DATA_DIR) -> Path:
    """Download a PGLib-OPF case file on first use (same source as the original script)."""
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / f"{name}.m"
    if not path.exists():
        urllib.request.urlretrieve(PGLIB_URL.format(name=name), path)
    return path


def fetch_matpower(name: str, data_dir: Path = DATA_DIR) -> Path:
    """Download a case file of MATPOWER 7.1 (e.g. the radial feeders case33bw, case69, case85, case141)."""
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / f"{name}.m"
    if not path.exists():
        urllib.request.urlretrieve(MATPOWER_URL.format(name=name), path)
    return path


# Executable statements at the end of the MATPOWER distribution feeders (case33bw, case69, case85, case141, ...).
# The data matrices of those files hold ohms and kW; these statements turn them into p.u. and MW. A plain
# parser of the matrices skips them, so they are applied here, in file order.
_M_OHM = "mpc.branch(:,[BR_R BR_X])=mpc.branch(:,[BR_R BR_X])/(Vbase^2/Sbase);"
_M_KW = "mpc.bus(:,[PD,QD])=mpc.bus(:,[PD,QD])/1e3;"
_M_QPF = "mpc.bus(:,QD)=mpc.bus(:,PD)*sin(acos(pf));"
_M_PPF = "mpc.bus(:,PD)=mpc.bus(:,PD)*pf;"
_M_VBASE = "Vbase=mpc.bus(1,BASE_KV)*1e3;"
_M_SBASE = "Sbase=mpc.baseMVA*1e6;"


def _squeeze(line: str) -> str:
    """A MATLAB statement without its comment and with single spaces only inside brackets."""
    code = line.split("%")[0].strip()
    code = re.sub(r"\s+", " ", code)
    return re.sub(r"\s*([=/*^(),;:])\s*", r"\1", code)


def matpower_unit_statements(text: str) -> list[str]:
    """Statements of a MATPOWER case file that modify ``mpc`` after the data matrices; raises on unknown ones."""
    lines = [_squeeze(ln) for ln in text.splitlines()]
    found = [ln for ln in lines if re.match(r"mpc\.\w+\(", ln)]
    unknown = [ln for ln in found if ln not in (_M_OHM, _M_KW, _M_QPF, _M_PPF)]
    if unknown:
        raise ValueError(f"the case file modifies its data with statements this loader does not know: {unknown}")
    if _M_OHM in found and not (_M_VBASE in lines and _M_SBASE in lines):
        raise ValueError("ohm-to-p.u. conversion without the expected definitions of Vbase and Sbase")
    return found


def apply_matpower_units(ppc: dict, text: str) -> list[str]:
    """Apply the unit conversions written as MATLAB code in a MATPOWER case file to ``ppc`` (in place)."""
    from pandapower.pypower.idx_brch import BR_R, BR_X
    from pandapower.pypower.idx_bus import BASE_KV, PD, QD

    done = matpower_unit_statements(text)
    pf = None
    m = re.search(r"^\s*pf\s*=\s*([0-9.]+)\s*;", text, flags=re.M)
    if m:
        pf = float(m.group(1))
    bus, branch = ppc["bus"], ppc["branch"]
    for stmt in done:
        if stmt == _M_OHM:
            z_base = (bus[0, BASE_KV] * 1e3) ** 2 / (ppc["baseMVA"] * 1e6)
            branch[:, [BR_R, BR_X]] = branch[:, [BR_R, BR_X]] / z_base
        elif stmt == _M_KW:
            bus[:, [PD, QD]] = bus[:, [PD, QD]] / 1e3
        else:
            if pf is None:
                raise ValueError("the case file scales the loads with a power factor that could not be read")
            if stmt == _M_QPF:
                bus[:, QD] = bus[:, PD] * np.sin(np.arccos(pf))
            else:
                bus[:, PD] = bus[:, PD] * pf
    return done


def _build_net(case: str, data_dir: Path):
    import pandapower.networks as pn

    if case.startswith("pglib_") or case.startswith("matpower_"):
        from pandapower.converter.matpower.from_mpc import _m2ppc
        from pandapower.converter.pypower.from_ppc import from_ppc

        path = fetch_pglib(case, data_dir) if case.startswith("pglib_") else fetch_matpower(case[len("matpower_"):], data_dir)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
        origin = "PGLib-OPF" if case.startswith("pglib_") else "MATPOWER"
        ppc = _m2ppc(str(path))
        done = apply_matpower_units(ppc, path.read_text(encoding="utf-8", errors="replace"))
        note = f" ({len(done)} unit statements of the file applied)" if done else ""
        return from_ppc(ppc, f_hz=50), f"{origin} {path.name} sha256:{digest}{note}"
    if not hasattr(pn, case):
        raise ValueError(f"unknown case {case!r}")
    return getattr(pn, case)(), f"pandapower.networks.{case}"


def load_grid(case: str, data_dir: Path = DATA_DIR, tol: float = 1e-9, prepare=None) -> Grid:
    """Build a :class:`Grid` from a pandapower built-in case, a ``pglib_opf_*`` file or ``matpower_<case>``.

    ``prepare``: optional function applied to the pandapower net before it is solved (e.g. to close tie lines)."""
    import pandapower as pp
    from pandapower.pypower.idx_brch import BR_STATUS, F_BUS, T_BUS
    from pandapower.pypower.idx_bus import BASE_KV, PD, QD
    from pandapower.pypower.idx_gen import GEN_BUS, GEN_STATUS, VG
    from pandapower.pypower.makeSbus import makeSbus

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net, source = _build_net(case, data_dir)
        if prepare is not None:
            prepare(net)
            source += " (modified)"
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, max_iteration=50)
        t0 = perf_counter()                       # second, warm call: the first one pays import and JIT costs
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, max_iteration=50)
        t_pf = perf_counter() - t0
    ppci = net._ppc["internal"]
    bus, gen, branch = ppci["bus"], ppci["gen"], ppci["branch"]
    base = float(ppci["baseMVA"])
    Ybus = sp.csr_matrix(ppci["Ybus"])
    ref, pv, pq = (np.asarray(ppci[k], dtype=int) for k in ("ref", "pv", "pq"))
    V0 = np.asarray(ppci["V"], dtype=complex)

    s_load = (bus[:, PD] + 1j * bus[:, QD]) / base
    s_gen = np.asarray(makeSbus(base, bus, gen)).ravel() + s_load

    v_set = np.ones(len(bus))
    on = gen[:, GEN_STATUS] > 0
    v_set[gen[on, GEN_BUS].astype(int)] = gen[on, VG]
    v_set[pq] = 1.0

    live = branch[:, BR_STATUS].real > 0
    br_f = branch[live, F_BUS].real.astype(int)
    br_t = branch[live, T_BUS].real.astype(int)

    lookup = np.asarray(net._pd2ppc_lookups["bus"])
    label = np.full(len(bus), -1)
    pp_idx = net.bus.index.values
    valid = lookup[pp_idx] < len(bus)
    label[lookup[pp_idx][valid][::-1]] = pp_idx[valid][::-1]     # first pandapower bus wins on fused buses
    vm_pp = net.res_bus.vm_pu.values[np.searchsorted(pp_idx, label)]
    if not np.allclose(vm_pp, np.abs(V0), atol=1e-9, equal_nan=False):
        raise RuntimeError("internal bus numbering does not match the pandapower result table")

    pvpq = np.concatenate([pv, pq])
    grid = Grid(case, base, Ybus, ref, pv, pq, V0, s_gen, s_load, v_set, br_f, br_t, label, t_pf, source,
                np.asarray(bus[:, BASE_KV].real, dtype=float),
                Yf=sp.csr_matrix(ppci["Yf"])[np.flatnonzero(live)] if "Yf" in ppci else None,
                Yt=sp.csr_matrix(ppci["Yt"])[np.flatnonzero(live)] if "Yt" in ppci else None)
    res = np.max(np.abs(stack_mismatch(mismatch(Ybus, V0, grid.sbus()), pvpq, pq)))
    if res > 1e-5:
        raise RuntimeError(f"extracted model is inconsistent with the pandapower solution (mismatch {res:.2e})")
    V_pol, ok, _, _ = newton_pf(Ybus, grid.sbus(), V0, pv, pq, tol=tol)   # polish to the solver tolerance
    if not ok:
        raise RuntimeError("base case did not converge to the requested tolerance")
    grid.V0 = V_pol
    return grid
