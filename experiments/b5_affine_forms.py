"""B5 -- Outputs from the affine forms, without the ellipse: per-bus intervals and branch voltage drops.

Targets: every PQ bus |V_i| and the drop |V_i| - |V_j| of every PQ-PQ branch, each with its own affine form
c_t = W' C (cdhim.affine_forms). For every target:

  deterministic (worst case over the noise box)
    det_lin     +-sum|c_tk|                                    (first order only)
    det_v       det_lin plus the quadratic term at the first-order extreme vertices +-sign(c_t)
    det_s       the extremes of the second-order model found by a flip search from those vertices
  probabilistic (alpha = 1 - confidence, per target)
    bern        Bernstein radius of the linear part plus the sampled range of the quadratic term
    sig3        3 sigma plus the same range (reference: assumes normality)
  for the drops, also the naive subtraction of two bus intervals (what one gets without the affine forms):
    sub_det_v, sub_bern, sub_cdhim (CD-HIM output intervals corrected as in B4)

Checks: (1) AC Monte Carlo with the samples of B1; (2) the exact AC power flow at the vertices found by both
the first-order rule and the second-order search, against the deterministic bounds.

    python experiments/b5_affine_forms.py              # design point of both systems
Writes results/b5_mc.csv, results/b5_vertex.csv, results/b5_width.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, sample_eps, save_table
from cdhim.affine_forms import (TargetModel, affine_radius, bernstein_radius, bus_targets, drop_targets, sigma)
from cdhim.correction import second_order_correction
from cdhim.validation import VIOLATION_ATOL

EXP = CFG["experiments"]
N_SURROGATE = 4000
ALPHA = 1.0 - CFG["method"]["confidence"]


def sampled_range(so, W, n, seed, batch=500):
    """min and max over samples of the quadratic part of every target."""
    rng = np.random.default_rng(seed)
    Wd = W.toarray()
    lo, hi = np.full(W.shape[1], np.inf), np.full(W.shape[1], -np.inf)
    for s in range(0, n, batch):
        q = Wd.T @ so.quadratic(sample_eps(rng, so.c_lin.shape[1], min(batch, n - s)))
        np.minimum(lo, q.min(axis=1), out=lo)
        np.maximum(hi, q.max(axis=1), out=hi)
    return lo, hi


def build(st, delta, res, seed):
    """Second-order models, bounds and worst vertices for the bus and drop targets."""
    t = {}
    t0 = perf_counter()
    so = SecondOrder(st, delta)
    t["second_order_build_s"] = perf_counter() - t0
    n_pq = so.c_lin.shape[0]
    v = res.valid
    fi, ti = res.fi[v], res.ti[v]
    out = {}
    for name, W in (("bus", bus_targets(n_pq)), ("drop", drop_targets(n_pq, fi, ti))):
        t0 = perf_counter()
        tm = TargetModel(so, W)
        Ct = tm.c
        r_aff, r_bern, sg = affine_radius(Ct), bernstein_radius(Ct, ALPHA), sigma(Ct)
        t1 = perf_counter()
        qmin, qmax = sampled_range(so, W, N_SURROGATE, seed + (1 if name == "bus" else 2))
        t_range = perf_counter() - t1
        t1 = perf_counter()
        Xp, fp, fpl, nfp = tm.worst_vertices(+1)
        Xm, fm, fml, nfm = tm.worst_vertices(-1)
        t_search = perf_counter() - t1
        qv_hi, qv_lo = fpl - r_aff, fml + r_aff                       # quadratic part at the first-order vertices
        b = {"det_lin": (-r_aff, r_aff),
             "det_v": (-r_aff + np.minimum(0.0, qv_lo), r_aff + np.maximum(0.0, qv_hi)),
             "det_s": (np.minimum(-r_aff, fm), np.maximum(r_aff, fp)),
             "bern": (-r_bern + np.minimum(0.0, qmin), r_bern + np.maximum(0.0, qmax)),
             "sig3": (-3 * sg + np.minimum(0.0, qmin), 3 * sg + np.maximum(0.0, qmax))}
        Xlin_p = np.sign(Ct).T
        Xlin_p[Xlin_p == 0] = 1.0
        out[name] = dict(W=W, tm=tm, bounds=b, r_aff=r_aff, r_bern=r_bern, sigma=sg,
                         vertices={"lin+": Xlin_p, "lin-": -Xlin_p, "so+": Xp, "so-": Xm},
                         model_gain_max=fp - fpl, model_gain_min=fml - fm, flips=(nfp, nfm))
        t[f"{name}_targets_s"] = perf_counter() - t0
        t[f"{name}_range_s"] = t_range
        t[f"{name}_search_s"] = t_search
    # naive subtraction of bus intervals for the drops
    bb = out["bus"]["bounds"]
    cor = second_order_correction(st, delta, res, N_SURROGATE, seed + 3, enclosure=False, vertex_buses=~res.bus_ellipse)
    cd = cor.cdhim(res, "range")
    sub = {"sub_det_v": bb["det_v"], "sub_bern": bb["bern"], "sub_cdhim": cd}
    for k, (lo, hi) in sub.items():
        out["drop"]["bounds"][k] = (lo[fi] - hi[ti], hi[fi] - lo[ti])
    out["bus"]["bounds"]["cdhim"] = cd
    return so, out, t


def ac_values(st, X, W, delta, batch=500):
    """Target values of the exact AC solution at column t of X (its own vertex). NaN where not converged."""
    Wd = W.tocsc()
    vals = np.full(X.shape[1], np.nan)
    for s in range(0, X.shape[1], batch):
        cols = np.arange(s, min(s + batch, X.shape[1]))
        V, conv, _ = st.solver.solve(st.injections(np.ascontiguousarray(X[:, cols]), delta))
        dv = np.abs(V[st.grid.pq]) - st.vm_star[:, None]
        tv = np.asarray(Wd[:, cols].multiply(dv).sum(axis=0)).ravel()
        vals[cols] = np.where(conv, tv, np.nan)
    return vals


def run(key):
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], EXP["n_ac"]
    seed = seed_for("ac_design_point", key)
    st = study(sysc["case"], threads=worker_threads())
    res = st.run(delta, method())
    so, out, timing = build(st, delta, res, seed)
    base = {"system": key, "delta": delta, "seed": seed, "n_surrogate": N_SURROGATE, "alpha": ALPHA}

    # (1) AC Monte Carlo, same samples as B1
    t0 = perf_counter()
    rng = np.random.default_rng(seed)
    C = st.C1 * delta
    counts = {(tg, k, ref): np.zeros(o["W"].shape[1], np.int64) for tg, o in out.items() for k in o["bounds"] for ref in ("ac", "lin")}
    anyc = {key_: 0 for key_ in counts}
    n_ac = n_lin = 0
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        lin = C @ eps
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_lin += b
        n_ac += int(conv.sum())
        for tg, o in out.items():
            Wt = o["W"].T.tocsr()
            for ref, dv in (("ac", ac), ("lin", lin)):
                tv = Wt @ dv
                for k, (lo, hi) in o["bounds"].items():
                    bad = (tv < lo[:, None] - VIOLATION_ATOL) | (tv > hi[:, None] + VIOLATION_ATOL)
                    counts[(tg, k, ref)] += bad.sum(axis=1)
                    anyc[(tg, k, ref)] += int(bad.any(axis=0).sum())
    t_mc = perf_counter() - t0
    mc_rows = []
    for (tg, k, ref), c in counts.items():
        nn = n_ac if ref == "ac" else n_lin
        r = c / max(nn, 1)
        mc_rows.append({**base, "target": tg, "bound": k, "reference": ref, "n_samples": nn, "n_targets": r.size,
                        "mean_pct": pct(r.mean()), "p99_pct": pct(np.quantile(r, 0.99)), "max_pct": pct(r.max()),
                        "n_over_alpha": int((r > ALPHA).sum()), "system_level_pct": pct(anyc[(tg, k, ref)] / max(nn, 1))})

    # (2) exact AC at the vertices
    t0 = perf_counter()
    vx_rows = []
    for tg, o in out.items():
        lo_s, hi_s = o["bounds"]["det_s"]
        for vname, X in o["vertices"].items():
            v = ac_values(st, X, o["W"], delta)
            ok = np.isfinite(v)
            side = +1 if vname.endswith("+") else -1
            row = {**base, "target": tg, "vertex": vname, "n_targets": int(v.size), "n_not_converged": int((~ok).sum())}
            for k in ("det_lin", "det_v", "det_s"):
                lo, hi = o["bounds"][k]
                exc = (v - hi) if side > 0 else (lo - v)              # > 0: the AC value is outside
                exc = exc[ok]
                width = (hi - lo)[ok]
                row[f"{k}_n_outside"] = int((exc > VIOLATION_ATOL).sum())
                row[f"{k}_max_excess_pu"] = float(max(exc.max(), 0.0)) if exc.size else np.nan
                rel = exc / np.where(width > 0, width, np.inf)
                row[f"{k}_max_excess_over_width"] = float(max(rel.max(), 0.0)) if exc.size else np.nan
            vx_rows.append(row)
    t_vertex = perf_counter() - t0

    # widths
    w_rows = []
    for tg, o in out.items():
        w = {k: hi - lo for k, (lo, hi) in o["bounds"].items()}
        ref = w["det_lin"]
        ok = ref > 1e-10
        for k, wk in w.items():
            for denom in ("det_lin", "bern", "sig3"):
                d = w[denom]
                m = ok & (d > 1e-10)
                w_rows.append({**base, "target": tg, "bound": k, "relative_to": denom,
                               "ratio_p10": float(np.quantile(wk[m] / d[m], 0.1)), "ratio_median": float(np.median(wk[m] / d[m])),
                               "ratio_p90": float(np.quantile(wk[m] / d[m], 0.9)), "median_width_pu": float(np.median(wk[ok]))})
        flips = np.concatenate(o["flips"])
        w_rows.append({**base, "target": tg, "bound": "search", "relative_to": "",
                       "ratio_p10": float(np.mean(flips > 0)), "ratio_median": float(np.median(o["model_gain_max"] / np.where(o["r_aff"] > 0, o["r_aff"], np.inf))),
                       "ratio_p90": float(np.max(np.concatenate([o["model_gain_max"], o["model_gain_min"]]) / np.where(np.tile(o["r_aff"], 2) > 0, np.tile(o["r_aff"], 2), np.inf))),
                       "median_width_pu": np.nan})
    timing.update(mc_s=t_mc, vertex_s=t_vertex)
    for r in w_rows + vx_rows + mc_rows:
        r.update({f"t_{k}": round(v, 3) for k, v in timing.items()})
    print(f"[b5] {key}: {timing}", flush=True)
    return mc_rows, vx_rows, w_rows


def main():
    out = parallel(run, list(CFG["systems"]), workers=2)
    meta = {"n_samples": EXP["n_ac"], "n_surrogate": N_SURROGATE, "alpha": ALPHA,
            "note": "search row of b5_width: ratio_p10 = share of targets where the search flipped any symbol, "
                    "ratio_median = median model gain of the search over the affine radius, ratio_p90 = its maximum"}
    print(save_table("b5_mc", [r for o in out for r in o[0]], meta))
    print(save_table("b5_vertex", [r for o in out for r in o[1]], meta))
    print(save_table("b5_width", [r for o in out for r in o[2]], meta))


if __name__ == "__main__":
    main()
