"""C7 -- Scalability: run time by stage on systems of growing size, fill-in, and log-log exponents.

Stages (single thread, each repeated ``timing_repeats`` times):
  jacobian  sparse dS/dV at the converged point and assembly of J
  lu        SuperLU factorisation (scipy.sparse.linalg.splu)
  solves    N_u sparse triangular solves that give the |V| rows of J^-1 for the uncertain injections
  branch    covariances, 2x2 spectra, decision and radii for every PQ-PQ branch
The base power flow that provides the linearisation point is timed separately (pandapower and the
in-house Newton from the flat start), because the method needs it but it is not part of the method.

Run it alone: other jobs on the machine distort the timings.
Writes results/c7_timing.csv, results/c7_fill.csv and results/c7_fit.csv
"""
import os
from time import perf_counter

import numpy as np
from scipy import stats
from scipy.sparse.linalg import splu

from common import CFG, grid, method
from cdhim import save_table
from cdhim.core import run_cdhim
from cdhim.powerflow import newton_pf, polar_jacobian
from cdhim.study import build_sources, linearize, voltage_coefficients

REPEATS = CFG["experiments"]["timing_repeats"]
STAGES = ("jacobian", "lu", "solves", "branch", "total")


def time_case(case):
    g = grid(case)
    delta = CFG["transfer_delta"].get(case, CFG["transfer_delta"]["default"])
    src = build_sources(g)
    fi, ti, _ = g.pq_pq_branches()
    cfg = method()
    t = {s: [] for s in STAGES}
    for _ in range(REPEATS):
        lin = linearize(g, g.V0, "exact", CFG["linearization"]["permc_spec"])
        t0 = perf_counter()
        C1 = voltage_coefficients(lin, g, src, threads=1)
        t1 = perf_counter()
        res = run_cdhim(C1 * delta, fi, ti, cfg)
        t2 = perf_counter()
        for s, v in zip(STAGES, (lin.t_jacobian, lin.t_lu, t1 - t0, t2 - t1, lin.t_jacobian + lin.t_lu + t2 - t0)):
            t[s].append(v)
    t0 = perf_counter()
    C1 = voltage_coefficients(lin, g, src, threads=None)             # default thread pool, for information
    t_threads = perf_counter() - t0
    t0 = perf_counter()
    _, ok, iters, _ = newton_pf(g.Ybus, g.sbus(), g.flat_voltage(), g.pv, g.pq, tol=CFG["powerflow"]["tol_pu"], max_iter=50)
    t_newton = perf_counter() - t0
    row = {"case": case, "n_bus": g.n_bus, "n_pq": len(g.pq), "n_pv": len(g.pv), "n_branch": len(g.br_f),
           "jacobian_dim": lin.J.shape[0], "N_u": src.n, "n_pq_pq_branches": int(res.n_branches),
           "n_ellipse": res.n_ellipse, "repeats": REPEATS}
    for s in STAGES:
        row[f"{s}_mean_ms"] = 1e3 * float(np.mean(t[s]))
        row[f"{s}_std_ms"] = 1e3 * float(np.std(t[s], ddof=1))
    row.update({"solves_share_of_total_pct": 100 * np.mean(t["solves"]) / np.mean(t["total"]),
                "solves_default_threads_ms": 1e3 * t_threads, "solve_threads_default": max(1, (os.cpu_count() or 2) // 2),
                "base_pf_pandapower_ms": 1e3 * g.base_pf_seconds, "base_pf_newton_flat_ms": 1e3 * t_newton,
                "base_pf_newton_flat_converged": ok, "base_pf_newton_iterations": iters,
                "dense_block_entries": int(C1.size), "dense_block_MB": C1.nbytes / 2**20})
    fill = []
    for spec in ("COLAMD", "MMD_AT_PLUS_A", "MMD_ATA"):
        tt = []
        for _ in range(3):
            t0 = perf_counter()
            lu = splu(lin.J, permc_spec=spec)
            tt.append(perf_counter() - t0)
        fill.append({"case": case, "n_bus": g.n_bus, "permc_spec": spec, "nnz_J": lin.J.nnz, "nnz_L": lu.L.nnz, "nnz_U": lu.U.nnz,
                     "fill_ratio_(L+U)/J": (lu.L.nnz + lu.U.nnz) / lin.J.nnz, "lu_time_ms": 1e3 * float(np.median(tt))})
    print(f"[c7] {case}: total {row['total_mean_ms']:.1f} ms", flush=True)
    return row, fill


def loglog_fit(x, y):
    """OLS of log y on log x; exponent, 95 % CI (Student t), R^2."""
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    r = stats.linregress(lx, ly)
    half = stats.t.ppf(0.975, len(lx) - 2) * r.stderr
    return r.slope, r.slope - half, r.slope + half, r.rvalue**2


def main():
    rows, fill = [], []
    for case in CFG["scaling_cases"]:
        try:
            r, f = time_case(case)
        except MemoryError as exc:
            print(f"[c7] {case} skipped: {exc}", flush=True)
            continue
        rows.append(r)
        fill += f
    fits = []
    n = [r["n_bus"] for r in rows]
    series = {s: [r[f"{s}_mean_ms"] for r in rows] for s in STAGES}
    series["nnz(L+U), COLAMD"] = [f["nnz_L"] + f["nnz_U"] for f in fill if f["permc_spec"] == "COLAMD"]
    series["nnz(J)"] = [f["nnz_J"] for f in fill if f["permc_spec"] == "COLAMD"]
    series["dense |V| block entries (n_pq x N_u)"] = [r["dense_block_entries"] for r in rows]
    for subset, idx in (("all systems", range(len(n))), ("n_bus >= 118", [i for i, v in enumerate(n) if v >= 118])):
        for name, y in series.items():
            k, lo, hi, r2 = loglog_fit([n[i] for i in idx], [y[i] for i in idx])
            fits.append({"quantity": name, "systems": subset, "n_points": len(list(idx)), "exponent": k,
                         "ci95_low": lo, "ci95_high": hi, "r2": r2, "ci_excludes_1.2": not (lo <= 1.2 <= hi)})
    meta = {"repeats": REPEATS, "threads": 1, "permc_spec": CFG["linearization"]["permc_spec"],
            "note": "x = number of buses; with 5-6 systems the confidence intervals are wide by construction"}
    print(save_table("c7_timing", rows, meta))
    print(save_table("c7_fill", fill, meta))
    print(save_table("c7_fit", fits, meta))


if __name__ == "__main__":
    main()
