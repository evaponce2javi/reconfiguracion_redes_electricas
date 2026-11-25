"""X13 -- How far the verified enclosure goes: network size, number of symbols, cost, and what limits it.

Networks (all non-slack buses PQ), larger than those of X8:
  MATPOWER 7.1 distribution cases with their tie lines closed (case118zh, case136ma), case141, the 533-bus network of
  MATPOWER 8.0 (case533mt_hi, radial and with ties) and two pandapower networks with every switch closed
  (mv_oberrhein, lv_schutterwald).

  case533mt_hi is not in MATPOWER 7.1 and writes arithmetic expressions (50/3, 135/sqrt(3)) inside its matrices,
  which the parser of pandapower does not evaluate. ``numeric_case`` downloads the file of the 8.0 tag and writes
  a copy with those expressions replaced by their values (data/<case>_numeric.m); nothing else is changed.
  lv_schutterwald comes with an unsorted bus table, which cdhim.network.load_grid does not support; the table is
  sorted by index before loading (same buses, same numbers).

Uncertainty: half disks at the sources, as in X8 (rule of Sou et al., Sec. VI-I: equal ratings, total rating equal
to |sum of the loads|), scenarios F and R, two shares of uncertain capacity. Two ways of placing the sources:
  every tenth bus   the rule of X8: one source every ten non-reference buses, so K = 2 x sources grows with the network;
  8 sources         eight sources evenly spaced over the non-reference buses in label order: K = 16 in every network;
  32 sources        the same with thirty-two sources: K = 64 in every network.

Per case: the stages are timed separately (second-order model; verification of the tube; bounds of the magnitudes
with three bounds of the quadratic forms: term by term, RAA-K blocks of up to 3 symbols, RAA-K with all symbols) and
the peak memory of the process is recorded. Every case runs in its own process, so a case that runs out of memory
or of time is reported as such and does not stop the others. A case whose second-order tensors (m x K x K complex
numbers each) would not fit is not attempted and says so. After the timed stages the tube is verified once more
with no uncertainty at all: that radius is the floor set by rounding and by the residual of the centre.

Reference: only cheap attained values (nine boundary load flows and random samples), so width / attained here is
an upper estimate of the looseness, coarser than the one of X8 and X12.

Timings are wall-clock on a shared machine; rerun alone for the figures to quote.

Writes results/x13_scale.csv.

    python -X utf8 experiments/x13_scale.py              # every case (case136ma needs UTF-8 mode to be read)
    python -X utf8 experiments/x13_scale.py 533          # only the networks whose name contains the given text
"""
import json
import math
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from time import perf_counter

import numpy as np

from common import CFG, ROOT, seed_for
from cdhim import SecondOrder, load_grid, save_table
from cdhim.distflow import SCENARIOS, benchmark_ratings, net_injections
from cdhim.enclosure import injection_changes, state_forms
from cdhim.halfdisk import half_disk_study, sample_half_disks
from cdhim.network import DATA_DIR
from cdhim.powerflow import BatchACSolver
from cdhim.validation import VIOLATION_ATOL
from cdhim.verified import magnitude_bounds, verify_tube

TOL = CFG["powerflow"]["tol_pu"]
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]
RHO_U = (0.1, 0.5)                    # shares of uncertain capacity
FIXED = (8, 32)                       # sources of the designs with a fixed number of them
N_SAMPLES = 400                       # random injections of the cheap reference
MIN_WIDTH = 1e-9                      # buses with a narrower attained range are left out of the ratios [p.u.]
TIME_LIMIT_S = 3600                   # per case; what is finished by then is kept
TENSOR_LIMIT_GB = 1.0                 # a single (m, K, K) complex tensor above this is not attempted (see TENSOR_COPIES)
TENSOR_COPIES = 12                    # tensors of that size alive at once, measured on the 533-bus case (column peak_rss_gb)
DENSE_COPIES = 15                     # dense 2m x 2m matrices alive at once in the tube, measured on Schutterwald (idem)
MATPOWER8_URL = "https://raw.githubusercontent.com/MATPOWER/matpower/8.0/data/{name}.m"
QUAD_BOUNDS = (("termwise", 1, 0), ("k3", 3, NODE_BUDGET), ("full", None, NODE_BUDGET))


def evaluate_expressions(text: str) -> str:
    """A MATPOWER case file with the expressions a/b and a/sqrt(b) of its code replaced by their values."""
    num = r"\d+(?:\.\d+)?"
    out = []
    for line in text.splitlines():
        code, sep, comment = line.partition("%")
        code = re.sub(rf"(-?{num})/sqrt\(({num})\)", lambda m: repr(float(m.group(1)) / math.sqrt(float(m.group(2)))), code)
        code = re.sub(rf"(?<![\w.])(-?{num})/({num})(?![\w.(])", lambda m: repr(float(m.group(1)) / float(m.group(2))), code)
        out.append(code + sep + comment)
    return "\n".join(out) + "\n"


def numeric_case(name: str, data_dir: Path = DATA_DIR) -> str:
    """Download ``name`` from the MATPOWER 8.0 tag on first use and write data/<name>_numeric.m; returns that name."""
    data_dir.mkdir(parents=True, exist_ok=True)
    raw, out = data_dir / f"{name}.m", data_dir / f"{name}_numeric.m"
    if not raw.exists():
        urllib.request.urlretrieve(MATPOWER8_URL.format(name=name), raw)
    if not out.exists():
        text = evaluate_expressions(raw.read_text(encoding="utf-8"))
        out.write_text(text.replace(f"function mpc = {name}", f"function mpc = {name}_numeric"), encoding="utf-8")
    return f"{name}_numeric"


def close_ties(net):
    net.line["in_service"] = True


def close_switches(net):
    net.switch["closed"] = True


def sort_and_close_switches(net):
    net.bus.sort_index(inplace=True)
    net.switch["closed"] = True


PREPARE = {"": None, "ties": close_ties, "switches": close_switches, "sorted, switches": sort_and_close_switches}

# (name in the tables, case, key of PREPARE, description)
NETWORKS = [("118 + ties", "matpower_case118zh", "ties", "MATPOWER 7.1 case118zh with its 15 tie lines closed"),
            ("136 + ties", "matpower_case136ma", "ties", "MATPOWER 7.1 case136ma with its 21 tie lines closed"),
            ("141", "matpower_case141", "", "MATPOWER 7.1 case141, radial"),
            ("Oberrhein", "mv_oberrhein", "switches", "pandapower mv_oberrhein, 20 kV, all switches closed"),
            ("533", "numeric:case533mt_hi", "", "MATPOWER 8.0 case533mt_hi, radial"),
            ("533 + ties", "numeric:case533mt_hi", "ties", "MATPOWER 8.0 case533mt_hi with its 45 tie lines closed"),
            ("Schutterwald", "lv_schutterwald", "sorted, switches", "pandapower lv_schutterwald, 0.4 kV with 14 MV/LV transformers, all switches closed")]
DESIGNS = ("every tenth bus",) + tuple(f"{n} sources" for n in FIXED)


def load(case, prepare):
    if case.startswith("numeric:"):
        case = "matpower_" + numeric_case(case[len("numeric:"):])
    return load_grid(case, tol=TOL, prepare=PREPARE[prepare])


def fixed_ratings(grid, n: int) -> np.ndarray:
    """``n`` sources evenly spaced over the non-reference buses in label order, equal ratings, total = |sum of loads|."""
    order = np.argsort(np.asarray(grid.bus_label))
    order = order[~np.isin(order, grid.ref)]
    dg = order[((np.arange(n) + 0.5) * len(order) / n).astype(int)]
    rating = np.zeros(grid.n_bus)
    rating[dg] = abs(grid.s_load.sum()) / len(dg)
    return rating


def peak_gb() -> float:
    import psutil

    info = psutil.Process().memory_info()
    return float(getattr(info, "peak_wset", info.rss)) / 2 ** 30


def run_case(spec: dict, out: Path):
    """One case, in this process; the row is written to ``out`` after every stage."""
    row = dict(spec)

    def save(**kw):
        row.update(kw)
        out.write_text(json.dumps(row, default=float), encoding="utf-8")

    save(status="loading the network")
    name, case, prepare = spec["network"], spec["case"], spec["prepare"]
    g = load(case, prepare)
    edges = {(min(a, b), max(a, b)) for a, b in zip(g.br_f, g.br_t) if a != b}
    rating = benchmark_ratings(g) if spec["design"] == DESIGNS[0] else fixed_ratings(g, int(spec["design"].split()[0]))
    sc = SCENARIOS[spec["scenario"]]
    s0, cap = net_injections(g, rating, sc["rho_D"], sc["rho_G"], sc["pf"], spec["rho_U"])
    dg = np.flatnonzero(cap > 0)
    m, K = len(g.pq), 2 * len(dg)
    tensor_gb = m * K * K * 16 / 2 ** 30
    save(n_bus=g.n_bus, slack_buses=len(g.ref), loops=len(edges) - g.n_bus + len(g.ref), largest_admittance_pu=float(abs(g.Ybus).max()),
         sources=len(dg), symbols=K, tensor_gb=tensor_gb, dense_matrix_gb=(2 * m) ** 2 * 8 / 2 ** 30)
    if len(g.pv):
        return save(status="not attempted: the network has PV buses")
    if tensor_gb > TENSOR_LIMIT_GB:
        return save(status=f"not attempted: one second-order tensor of {m} x {K} x {K} complex numbers is {tensor_gb:.1f} GB "
                           f"and about {TENSOR_COPIES} of them are alive at once")

    save(status="second-order model")
    t0 = perf_counter()
    try:
        st = half_disk_study(g, s0, cap, tol=TOL)
    except RuntimeError as e:
        return save(status=str(e))
    so = SecondOrder(st, 1.0)
    Hs = state_forms(so)
    ns = g.pq
    E, B = so.E[ns], Hs[:m] + 1j * Hs[m:]
    del Hs
    mis = st.V_star * np.conj(g.Ybus @ st.V_star) - st.sbus
    save(model_s=perf_counter() - t0, peak_rss_after_model_gb=peak_gb(), centre_mismatch_max_pu=float(np.abs(mis[ns]).max()),
         status="verification of the tube")

    t0 = perf_counter()
    tube = verify_tube(g.Ybus, ns, st.V_star, st.sbus[ns], injection_changes(so)[ns], E, B)
    save(tube_s=perf_counter() - t0, peak_rss_after_tube_gb=peak_gb(), verified=bool(tube.exists), unique=bool(tube.unique), reason=tube.reason,
         tube_iterations=tube.iterations, status="reference")
    if not tube.exists:
        return save(status="ok", peak_rss_gb=peak_gb())
    save(tube_radius_max_pu=float(tube.z.max()), rounding_share_of_tube=tube.rounding_share)

    # cheap reference: nine boundary load flows and random samples
    t0 = perf_counter()
    k = np.arange(9)
    pts = np.hstack([np.repeat((np.sin(k * np.pi / 8) + 1j * np.cos(k * np.pi / 8))[None, :], len(dg), axis=0),
                     sample_half_disks(np.random.default_rng(spec["seed"]), len(dg), N_SAMPLES)])
    pts.real = np.maximum(pts.real, 0.0)
    S = np.repeat(np.asarray(s0, dtype=complex)[:, None], pts.shape[1], axis=1)
    S[dg] += cap[dg][:, None] * pts
    V, conv, _ = BatchACSolver(g.Ybus, g.pv, g.pq, st.V_star, tol=TOL).solve(S)
    vm = np.abs(V[g.pq][:, conv])
    vmax, vmin = vm.max(axis=1), vm.min(axis=1)
    width = vmax - vmin
    act = width > MIN_WIDTH
    save(reference_s=perf_counter() - t0, reference_points=int(conv.sum()), reference_ac_failures=int((~conv).sum()), buses=int(act.sum()),
         attained_width_max_pu=float(width.max()))

    widths = {}
    for tag, max_k, budget in QUAD_BOUNDS:
        save(status=f"bounds of the magnitudes, quadratic forms {tag}")
        t0 = perf_counter()
        lo, hi = magnitude_bounds(st.V_star[ns], E, B, tube.z, True, max_k, budget)
        dt = perf_counter() - t0
        widths[tag] = hi - lo
        ratio = (widths[tag] / width)[act]
        save(**{f"magnitudes_{tag}_s": dt, f"width_over_attained_{tag}_median": float(np.median(ratio)), f"width_over_attained_{tag}_max": float(ratio.max()),
                f"violations_{tag}": int(((hi < vmax - VIOLATION_ATOL) | (lo > vmin + VIOLATION_ATOL)).sum())})
        if tag != "termwise":
            save(**{f"width_{tag}_over_termwise_median": float(np.median((widths[tag] / widths["termwise"])[act]))})
    peak = peak_gb()
    # the same verification with no uncertainty: what rounding and the residual of the centre leave of the tube
    save(status="tube without uncertainty", peak_rss_gb=peak)
    t0 = perf_counter()
    none = np.zeros((m, 0), complex)
    floor = verify_tube(g.Ybus, ns, st.V_star, st.sbus[ns], none, none, np.zeros((m, 0, 0), complex))
    if floor.exists:
        save(tube_radius_without_uncertainty_max_pu=float(floor.z.max()),
             floor_share_at_largest_radius=float(floor.z[int(np.argmax(tube.z))] / tube.z.max()))
    save(status="ok", floor_s=perf_counter() - t0, floor_verified=bool(floor.exists))


def cases(only):
    for name, case, prepare, _ in NETWORKS:
        if only and not any(o.lower() in name.lower() for o in only):
            continue
        for design in DESIGNS:
            for scen in "FR":
                for rho_U in RHO_U:
                    yield {"network": name, "case": case, "prepare": prepare, "design": design, "scenario": scen, "rho_U": rho_U,
                           "seed": seed_for("x13", name, design, scen, rho_U)}


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--one":
        run_case(json.loads(sys.argv[2]), Path(sys.argv[3]))
        return
    rows = []
    tmp = Path(tempfile.mkdtemp(prefix="x13_"))
    for i, spec in enumerate(cases(sys.argv[1:])):
        out = tmp / f"case_{i}.json"
        t0 = perf_counter()
        note = ""
        try:
            done = subprocess.run([sys.executable, "-X", "utf8", str(Path(__file__).resolve()), "--one", json.dumps(spec), str(out)],
                                  cwd=str(ROOT), timeout=TIME_LIMIT_S, capture_output=True, text=True)
            if done.returncode != 0:
                tail = (done.stderr or "").strip().splitlines()[-1:] or ["no message"]
                note = f"the process ended with code {done.returncode}: {tail[0][:200]}"
        except subprocess.TimeoutExpired:
            note = f"stopped after {TIME_LIMIT_S} s"
        row = json.loads(out.read_text(encoding="utf-8")) if out.exists() else dict(spec, status="nothing was written")
        if note:
            row["status"] = f"{note}; last stage: {row.get('status', '')}"
        row["wall_s"] = perf_counter() - t0
        row.pop("case", None)
        row.pop("prepare", None)
        rows.append(row)
        print(f"[x13] {row['network']} | {row['design']} | {row['scenario']} rho_U={row['rho_U']:g}: n={row.get('n_bus')} K={row.get('symbols')} "
              f"status={row['status']} verified={row.get('verified')} unique={row.get('unique')} model {row.get('model_s', float('nan')):.2f}s "
              f"tube {row.get('tube_s', float('nan')):.2f}s magnitudes {row.get('magnitudes_termwise_s', float('nan')):.2f}/{row.get('magnitudes_k3_s', float('nan')):.2f}/"
              f"{row.get('magnitudes_full_s', float('nan')):.2f}s width/attained {row.get('width_over_attained_termwise_median', float('nan')):.4f}/"
              f"{row.get('width_over_attained_k3_median', float('nan')):.4f}/{row.get('width_over_attained_full_median', float('nan')):.4f} "
              f"peak {row.get('peak_rss_gb', float('nan')):.2f} GB", flush=True)
    meta = {"script": "experiments/x13_scale.py", "pf_tol": TOL, "rho_U": list(RHO_U), "fixed_sources": list(FIXED), "samples": N_SAMPLES,
            "raak_node_budget": NODE_BUDGET, "time_limit_s": TIME_LIMIT_S, "tensor_limit_gb": TENSOR_LIMIT_GB, "tensor_copies": TENSOR_COPIES,
            "dense_copies": DENSE_COPIES,
            "min_width_pu": MIN_WIDTH, "networks": {n: d for n, _, _, d in NETWORKS},
            "quadratic_bounds": {"termwise": "no blocks", "k3": "RAA-K blocks of up to 3 symbols", "full": "RAA-K blocks of any size"},
            "note": "one process per case; timings are wall-clock and depend on what else runs on the machine"}
    if sys.argv[1:]:
        meta["only"] = sys.argv[1:]
    print(save_table("x13_scale" if not sys.argv[1:] else "x13_scale_partial", rows, meta))


if __name__ == "__main__":
    main()
