"""B4 -- Second-order correction of the CD-HIM bounds against nonlinear AC (design point and load x delta sweep).

B1 shows that, at the design point, every AC violation above the nominal 0.27 % comes from the linearisation
error: the same bounds tested on the linear model stay below it. This experiment adds the second-order term
of the voltage magnitude (cdhim.correction) and re-tests the bounds on the same AC samples.

Corrections (see cdhim/correction.py): none | shift (centre on E[q]) | moment (centre + variance of q) |
range (interval plus the sampled range of q) for intervals; minkowski (outer ellipse of E + range of q) for ellipses.
The deterministic affine box gets box2 (rigorous enclosure of q, loose) or box2v (q at the extreme vertices).
CD-HIM corrected = Ellipse buses with shift/moment/range, Box buses with box2v.

Every bound is also tested on the second-order surrogate (c'eps + q(eps)) to separate what the third and
higher orders still leave outside.

    python experiments/b4_second_order.py            # design point (same AC samples as B1) and the sweep
    python experiments/b4_second_order.py design     # design point only (also: sweep, timing)
Writes results/b4_design_bus.csv, results/b4_design_branch.csv, results/b4_sweep.csv, results/b4_timing.csv
"""
import sys
from time import perf_counter

import numpy as np

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import loadability_limit, sample_eps, save_table
from cdhim import geometry as geo
from cdhim.correction import second_order_correction
from cdhim.validation import VIOLATION_ATOL

EXP = CFG["experiments"]
N_SURROGATE = 4000          # samples of the quadratic term used for Var(q), Cov(q) and its range
ELLIPSE_KINDS = ("none", "shift", "moment", "minkowski")


def bounds_and_ellipses(res, res_c, cor):
    """name -> (lo, hi) offsets per PQ bus; name -> ellipse parameters per Ellipse branch."""
    b = {}
    for kind in ("none", "shift", "moment", "range"):
        b[f"cdhim_lmax|{kind}"] = cor.interval(res.r_out, res, kind) if kind == "none" else cor.cdhim(res, kind)
        b[f"cdhim_coord|{kind}"] = cor.interval(res_c.r_out, res_c, kind) if kind == "none" else cor.cdhim(res_c, kind)
        b[f"cltm_3sigma|{kind}"] = cor.interval(res.r_cltm, res, kind)
    b["affine_box|none"] = cor.interval(res.r_affine, res, "none")
    b["affine_box|box2"] = cor.interval(res.r_affine, res, "box2")
    b["affine_box|box2v"] = cor.interval(res.r_affine, res, "box2v")
    e = {kind: cor.ellipse(res, kind) for kind in ELLIPSE_KINDS}
    return b, e


class Counter:
    def __init__(self, bounds, ellipses, n_pq, n_br):
        self.n = 0
        self.bus = {k: np.zeros(n_pq, np.int64) for k in bounds}
        self.any = {k: 0 for k in bounds}
        self.br = {k: np.zeros(n_br, np.int64) for k in ellipses}
        self.any_br = {k: 0 for k in ellipses}

    def add(self, dv, bounds, ellipses, fi, ti, chi2):
        self.n += dv.shape[1]
        for k, (lo, hi) in bounds.items():
            out = (dv < lo[:, None] - VIOLATION_ATOL) | (dv > hi[:, None] + VIOLATION_ATOL)
            self.bus[k] += out.sum(axis=1)
            self.any[k] += int(out.any(axis=0).sum())
        for k, (cf, ct, q11, q22, q12) in ellipses.items():
            if not fi.size:
                continue
            out = geo.mahalanobis2(dv[fi] - cf[:, None], dv[ti] - ct[:, None], q11, q22, q12) > chi2
            self.br[k] += out.sum(axis=1)
            self.any_br[k] += int(out.any(axis=0).sum())

    @staticmethod
    def rates(counts, n):
        r = counts / max(n, 1)
        if r.size == 0:
            return dict(mean=np.nan, p99=np.nan, max=np.nan, n_over=0, n=0)
        return dict(mean=pct(r.mean()), p99=pct(np.quantile(r, 0.99)), max=pct(r.max()), n_over=int((r > 0.0027).sum()), n=int(r.size))


def evaluate(st, delta, n, seed, batch):
    """One Monte Carlo pass. Returns the result objects, the correction and the counters (AC, linear, surrogate)."""
    res = st.run(delta, method())
    res_c = st.run(delta, method(radius="coordinate"))
    t0 = perf_counter()
    cor = second_order_correction(st, delta, res, N_SURROGATE, seed + 1)
    t_cor = perf_counter() - t0
    bounds, ellipses = bounds_and_ellipses(res, res_c, cor)
    m = res.ellipse
    fi, ti = res.fi[m], res.ti[m]
    n_pq = st.C1.shape[0]
    cnt = {ref: Counter(bounds, ellipses, n_pq, int(m.sum())) for ref in ("ac", "lin", "so")}
    rng = np.random.default_rng(seed)
    C = st.C1 * delta
    n_drawn = 0
    for start in range(0, n, batch):
        b = min(batch, n - start)
        eps = sample_eps(rng, st.src.n, b)
        n_drawn += b
        lin = C @ eps
        so = lin + cor.so.quadratic(eps)
        cnt["lin"].add(lin, bounds, ellipses, fi, ti, res.chi2)
        cnt["so"].add(so, bounds, ellipses, fi, ti, res.chi2)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        if conv.any():
            dv = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
            cnt["ac"].add(dv, bounds, ellipses, fi, ti, res.chi2)
    return res, res_c, cor, bounds, ellipses, cnt, n_drawn, t_cor


def width_stats(bounds, base, mask=None):
    """Median and max of (corrected width / uncorrected width) over the buses with a nonzero width."""
    lo, hi = bounds
    blo, bhi = base
    w, w0 = hi - lo, bhi - blo
    ok = w0 > 1e-10 if mask is None else (w0 > 1e-10) & mask
    r = w[ok] / w0[ok]
    return float(np.median(r)), float(r.max())


# ----------------------------------------------------------------------------- design point
def run_design(key):
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], EXP["n_ac"]
    seed = seed_for("ac_design_point", key)                      # the AC samples of B1
    st = study(sysc["case"], threads=worker_threads())
    t0 = perf_counter()
    res, res_c, cor, bounds, ellipses, cnt, n_drawn, t_cor = evaluate(st, delta, n, seed, CFG["powerflow"]["batch"])
    base = {"system": key, "delta": delta, "n_samples": n_drawn, "seed": seed, "n_ac_converged": cnt["ac"].n,
            "n_surrogate": N_SURROGATE, "t_correction_s": round(t_cor, 3), "mean_check_pu": cor.mean_check}
    bus_rows, br_rows = [], []
    for name, bd in bounds.items():
        bound, kind = name.split("|")
        ref0 = bounds[f"{bound}|none"]
        med, mx = width_stats(bd, ref0)
        for ref in ("ac", "lin", "so"):
            s = Counter.rates(cnt[ref].bus[name], cnt[ref].n)
            bus_rows.append({**base, "bound": bound, "correction": kind, "reference": ref, "n_buses": s["n"],
                             "mean_pct": s["mean"], "p99_pct": s["p99"], "max_pct": s["max"], "n_buses_over_0.27pct": s["n_over"],
                             "system_level_pct": pct(cnt[ref].any[name] / max(cnt[ref].n, 1)),
                             "width_ratio_median": med, "width_ratio_max": mx})
    m = res.ellipse
    area0 = np.sqrt(np.maximum(res.q11[m] * res.q22[m] - res.q12[m] ** 2, 0))
    for kind, (cf, ct, q11, q22, q12) in ellipses.items():
        area = np.sqrt(np.maximum(q11 * q22 - q12 ** 2, 0))
        ok = area0 > 0
        for ref in ("ac", "lin", "so"):
            s = Counter.rates(cnt[ref].br[kind], cnt[ref].n)
            br_rows.append({**base, "correction": kind, "reference": ref, "n_ellipse_branches": s["n"], "mean_pct": s["mean"],
                            "p99_pct": s["p99"], "max_pct": s["max"], "n_branches_over_0.27pct": s["n_over"],
                            "system_level_pct": pct(cnt[ref].any_br[kind] / max(cnt[ref].n, 1)),
                            "area_ratio_median": float(np.median(area[ok] / area0[ok])) if ok.any() else np.nan,
                            "area_ratio_max": float(np.max(area[ok] / area0[ok])) if ok.any() else np.nan})
    print(f"[b4] design {key:11s} done in {perf_counter() - t0:.0f}s (correction {t_cor:.1f}s)", flush=True)
    return bus_rows, br_rows


# ----------------------------------------------------------------------------- load x delta sweep
SWEEP_BOUNDS = [("cdhim_lmax", "none"), ("cdhim_lmax", "shift"), ("cdhim_lmax", "moment"), ("cdhim_lmax", "range"),
                ("cdhim_coord", "none"), ("cdhim_coord", "moment"), ("cltm_3sigma", "none"), ("cltm_3sigma", "moment"),
                ("affine_box", "none"), ("affine_box", "box2v")]


def run_sweep(task):
    key, lf, tag, lam_max = task
    sysc = CFG["systems"][key]
    seed = seed_for("c1", key)                                   # the AC samples of C1
    rows = []
    try:
        st = study(sysc["case"], load_factor=lf, threads=worker_threads())
    except RuntimeError as exc:
        return [{"system": key, "load_factor": lf, "load_level": tag, "delta": d, "status": f"blocked: {exc}"} for d in EXP["deltas"]]
    for delta in EXP["deltas"]:
        t0 = perf_counter()
        res, _, cor, bounds, _, cnt, n_drawn, _ = evaluate(st, delta, EXP["n_sweep"], seed, CFG["powerflow"]["batch"])
        row = {"system": key, "load_factor": lf, "load_level": tag, "loadability_limit": lam_max, "delta": delta,
               "n_samples": n_drawn, "seed": seed, "status": "ok", "nonconvergence_pct": pct(1 - cnt["ac"].n / n_drawn)}
        for bound, kind in SWEEP_BOUNDS:
            name = f"{bound}|{kind}"
            for ref in ("ac", "lin", "so"):
                s = Counter.rates(cnt[ref].bus[name], cnt[ref].n)
                row[f"{ref}_{bound}_{kind}_max_pct"] = s["max"]
                row[f"{ref}_{bound}_{kind}_n_over"] = s["n_over"]
                row[f"{ref}_{bound}_{kind}_system_pct"] = pct(cnt[ref].any[name] / max(cnt[ref].n, 1))
            row[f"width_{bound}_{kind}_median"] = width_stats(bounds[name], bounds[f"{bound}|none"])[0]
        for kind in ELLIPSE_KINDS:
            for ref in ("ac", "lin", "so"):
                s = Counter.rates(cnt[ref].br[kind], cnt[ref].n)
                row[f"{ref}_ellipse_{kind}_max_pct"] = s["max"]
                row[f"{ref}_ellipse_{kind}_n_over"] = s["n_over"]
        rows.append(row)
        print(f"[b4] sweep {key:11s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s "
              f"nonconv={row['nonconvergence_pct']:.2f}%", flush=True)
    return rows


# ----------------------------------------------------------------------------- cost of the correction
SURROGATE_SIZES = (1000, 2000, 4000)


def run_timing(key, repeats=3):
    """Single-thread cost of each piece, as the method would run it (vertex values only at Box buses,
    no rigorous enclosure), next to the cost of the linear method itself."""
    from cdhim import SecondOrder
    from cdhim.correction import vertex_values
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    rows = []
    for rep in range(repeats):
        t0 = perf_counter()
        st = study(sysc["case"], threads=1)
        t_lin = st.timings["jacobian"] + st.timings["lu"] + st.timings["solves"]
        t1 = perf_counter()
        res = st.run(delta, method())
        t_branch = perf_counter() - t1
        t1 = perf_counter()
        so = SecondOrder(st, delta)
        t_build = perf_counter() - t1
        t1 = perf_counter()
        so.trace()
        t_trace = perf_counter() - t1
        t1 = perf_counter()
        vertex_values(so, ~res.bus_ellipse)
        t_vertex = perf_counter() - t1
        row = {"system": key, "repeat": rep, "linear_method_s": t_lin + t_branch, "second_order_build_s": t_build,
               "mean_trace_s": t_trace, "vertex_box_buses_s": t_vertex, "n_box_buses": int((~res.bus_ellipse).sum())}
        rng = np.random.default_rng(rep)
        for n in SURROGATE_SIZES:
            t1 = perf_counter()
            for s0 in range(0, n, 500):
                so.quadratic(sample_eps(rng, st.src.n, min(500, n - s0)))
            row[f"surrogate_{n}_s"] = perf_counter() - t1
        rows.append(row)
        print(f"[b4] timing {key} repeat {rep}: {row}", flush=True)
    return rows


def main():
    what = sys.argv[1:] or ["design", "sweep", "timing"]
    meta = {"n_surrogate": N_SURROGATE, "violation_atol_pu": VIOLATION_ATOL,
            "note": "offsets relative to the converged AC solution; the surrogate seed is the AC seed + 1; "
                    "reference so = linear + exact quadratic term (no power flow)"}
    if "design" in what:
        out = parallel(run_design, list(CFG["systems"]), workers=2)
        print(save_table("b4_design_bus", [r for b, _ in out for r in b], {**meta, "n_samples": EXP["n_ac"]}))
        print(save_table("b4_design_branch", [r for _, b in out for r in b], {**meta, "n_samples": EXP["n_ac"]}))
    if "sweep" in what:
        tasks = []
        for key, sysc in CFG["systems"].items():
            lam_max, _ = loadability_limit(grid(sysc["case"]), tol=CFG["powerflow"]["tol_pu"])
            levels = [(lf, "fixed") for lf in EXP["load_factors"]] + [(round(EXP["near_limit_fraction"] * lam_max, 4), "near_limit")]
            tasks += [(key, lf, tag, lam_max) for lf, tag in levels]
        rows = [r for out in parallel(run_sweep, tasks, workers=3) for r in out]
        print(save_table("b4_sweep", rows, {**meta, "n_samples_per_cell": EXP["n_sweep"]}))
    if "timing" in what:
        import os
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ[var] = "1"
        rows = [r for key in CFG["systems"] for r in run_timing(key)]
        print(save_table("b4_timing", rows, {"note": "single process; run alone. Thread variables are set after NumPy "
                                                     "is imported, so BLAS may still use several threads"}))


if __name__ == "__main__":
    main()
