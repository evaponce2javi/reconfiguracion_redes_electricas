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
    python experiments/b5_affine_forms.py sweep        # load x delta sweep of C1 (Monte Carlo only, no vertex search)
    python experiments/b5_affine_forms.py surrogate    # coverage and cost against the number of surrogate samples
    python experiments/b5_affine_forms.py radii        # Gaussian, Bernstein, Hoeffding, unimodal and Chebyshev radii
    python experiments/b5_affine_forms.py sweep250     # the sweep with 250 samples of the quadratic term
Writes results/b5_mc.csv, results/b5_vertex.csv, results/b5_width.csv (design) or results/b5_sweep.csv
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


def lin_vertex_values(tm, chunk=128):
    """f_t at the first-order vertices +sign(c_t) and -sign(c_t)."""
    nt = tm.c.shape[0]
    fp, fm = np.empty(nt), np.empty(nt)
    for s in range(0, nt, chunk):
        idx = np.arange(s, min(s + chunk, nt))
        X = np.sign(tm.c[idx]).T
        X[X == 0] = 1.0
        fp[idx], fm[idx] = tm.value(X, idx), tm.value(-X, idx)
    return fp, fm


def build(st, delta, res, seed, search=True):
    """Second-order models, bounds and worst vertices for the bus and drop targets.
    ``search=False``: no second-order vertex search (no det_s), for the sweep."""
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
        if search:
            Xp, fp, fpl, nfp = tm.worst_vertices(+1)
            Xm, fm, fml, nfm = tm.worst_vertices(-1)
        else:
            fpl, fml = lin_vertex_values(tm)
            fp, fm, Xp, Xm = fpl, fml, None, None
            nfp = nfm = np.zeros(fpl.size, int)
        t_search = perf_counter() - t1
        qv_hi, qv_lo = fpl - r_aff, fml + r_aff                       # quadratic part at the first-order vertices
        b = {"det_lin": (-r_aff, r_aff),
             "det_v": (-r_aff + np.minimum(0.0, qv_lo), r_aff + np.maximum(0.0, qv_hi)),
             **({"det_s": (np.minimum(-r_aff, fm), np.maximum(r_aff, fp))} if search else {}),
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


# ----------------------------------------------------------------------------- load x delta sweep (Monte Carlo only)
SWEEP = {"bus": ("det_lin", "det_v", "bern", "sig3", "cdhim"), "drop": ("det_v", "bern", "sig3", "sub_bern", "sub_cdhim")}


def run_sweep(task):
    key, lf, tag, lam_max = task
    sysc = CFG["systems"][key]
    seed = seed_for("c1", key)                                   # the AC samples of C1 and B4
    try:
        st = study(sysc["case"], load_factor=lf, threads=worker_threads())
    except RuntimeError as exc:
        return [{"system": key, "load_factor": lf, "load_level": tag, "delta": d, "status": f"blocked: {exc}"} for d in EXP["deltas"]]
    rows = []
    for delta in EXP["deltas"]:
        t0 = perf_counter()
        res = st.run(delta, method())
        so, out, timing = build(st, delta, res, seed, search=False)
        rng = np.random.default_rng(seed)
        C = st.C1 * delta
        counts = {(tg, k, ref): np.zeros(out[tg]["W"].shape[1], np.int64) for tg, ks in SWEEP.items() for k in ks for ref in ("ac", "lin")}
        anyc = {k_: 0 for k_ in counts}
        n_ac = n_lin = 0
        for s in range(0, EXP["n_sweep"], CFG["powerflow"]["batch"]):
            b = min(CFG["powerflow"]["batch"], EXP["n_sweep"] - s)
            eps = sample_eps(rng, st.src.n, b)
            lin = C @ eps
            V, conv, _ = st.solver.solve(st.injections(eps, delta))
            ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
            n_lin += b
            n_ac += int(conv.sum())
            for tg, ks in SWEEP.items():
                Wt = out[tg]["W"].T.tocsr()
                for ref, dv in (("ac", ac), ("lin", lin)):
                    tv = Wt @ dv
                    for k in ks:
                        lo, hi = out[tg]["bounds"][k]
                        bad = (tv < lo[:, None] - VIOLATION_ATOL) | (tv > hi[:, None] + VIOLATION_ATOL)
                        counts[(tg, k, ref)] += bad.sum(axis=1)
                        anyc[(tg, k, ref)] += int(bad.any(axis=0).sum())
        row = {"system": key, "load_factor": lf, "load_level": tag, "loadability_limit": lam_max, "delta": delta,
               "n_samples": n_lin, "seed": seed, "status": "ok", "nonconvergence_pct": pct(1 - n_ac / n_lin),
               "t_correction_s": round(timing["second_order_build_s"] + timing["bus_range_s"] + timing["drop_range_s"], 2)}
        for (tg, k, ref), c in counts.items():
            nn = n_ac if ref == "ac" else n_lin
            r = c / max(nn, 1)
            row[f"{ref}_{tg}_{k}_max_pct"] = pct(r.max())
            row[f"{ref}_{tg}_{k}_n_over"] = int((r > ALPHA).sum())
            row[f"{ref}_{tg}_{k}_system_pct"] = pct(anyc[(tg, k, ref)] / max(nn, 1))
        for tg, ks in SWEEP.items():
            w0 = out[tg]["bounds"]["det_lin"][1] - out[tg]["bounds"]["det_lin"][0]
            ok = w0 > 1e-10
            for k in ks:
                lo, hi = out[tg]["bounds"][k]
                row[f"width_{tg}_{k}_over_affine_median"] = float(np.median((hi - lo)[ok] / w0[ok]))
        rows.append(row)
        print(f"[b5] sweep {key:11s} lf={lf:.3f} delta={delta:.2f} {perf_counter() - t0:5.0f}s nonconv={row['nonconvergence_pct']:.2f}%", flush=True)
    return rows


# ----------------------------------------------------------------------------- how many surrogate samples are needed
SURROGATE_SIZES = (250, 500, 1000, 2000, 4000)
SURROGATE_SEEDS = 3


def run_surrogate(key):
    """Coverage of the Bernstein bounds (bus and drop) against the AC samples of B1 when the range of the
    quadratic term comes from fewer samples, and the single-process cost of the whole method."""
    from cdhim.affine_forms import quadratic_ranges
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], EXP["n_ac"]
    seed = seed_for("ac_design_point", key)
    t0 = perf_counter()
    st = study(sysc["case"], threads=worker_threads())
    t_study = perf_counter() - t0
    t0 = perf_counter()
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    t_so = perf_counter() - t0
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    rb = {k: bernstein_radius(np.asarray(W.T @ (st.C1 * delta)), ALPHA) for k, W in Ws.items()}
    bounds, times = {}, {}
    for m in SURROGATE_SIZES:
        for r in range(SURROGATE_SEEDS):
            t0 = perf_counter()
            rng_ = quadratic_ranges(so, Ws, m, seed + 100 + r)
            times[(m, r)] = perf_counter() - t0
            for k in Ws:
                lo, hi = rng_[k]
                bounds[(k, m, r)] = (-rb[k] + np.minimum(0.0, lo), rb[k] + np.maximum(0.0, hi))
    for k in Ws:                                                   # no correction, for reference
        bounds[(k, 0, 0)] = (-rb[k], rb[k])
    counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ ac for k in Ws}
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    rows = []
    for (k, m, r), c in counts.items():
        rate = c / max(n_ac, 1)
        lo, hi = bounds[(k, m, r)]
        w0 = 2 * rb[k]
        ok = w0 > 1e-10
        rows.append({"system": key, "delta": delta, "seed": seed, "n_ac": n_ac, "target": k, "n_surrogate": m, "surrogate_rep": r,
                     "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "system_level_pct": pct(anyc[(k, m, r)] / max(n_ac, 1)),
                     "width_over_uncorrected_median": float(np.median((hi - lo)[ok] / w0[ok])),
                     "t_range_s": round(times.get((m, r), 0.0), 3), "t_study_s": round(t_study, 3), "t_run_and_second_order_s": round(t_so, 3)})
    print(f"[b5] surrogate {key} done", flush=True)
    return rows


# ----------------------------------------------------------------------------- distribution-free radii compared
N_SUR_RADII = 1000


def radii(Ct, alpha):
    """Half-widths of the linear part at level alpha, all capped by the affine (worst-case) radius."""
    sg, r_aff = sigma(Ct), affine_radius(Ct)
    L = np.log(2.0 / alpha)
    raw = {
        "gauss_3sigma": 3.0 * sg,                                   # Gaussian, the manuscript's CLTM box (not distribution-free)
        "gauss_quantile": sg * float(np.sqrt(2.0) * __import__("scipy.special", fromlist=["erfinv"]).erfinv(1.0 - alpha)),
        "bernstein": bernstein_radius(Ct, alpha),                   # independence + support [-1, 1] + variance
        "hoeffding": np.sqrt(2.0 * L * (Ct ** 2).sum(axis=1)),      # independence + support only
        "unimodal_gauss_ineq": sg * np.sqrt(4.0 / (9.0 * alpha)),   # any unimodal distribution with this variance
        "chebyshev": sg / np.sqrt(alpha),                           # variance only
        "affine": r_aff,                                            # worst case of the linear model
    }
    return {k: np.minimum(v, r_aff) for k, v in raw.items()}


def run_radii(key):
    from cdhim.affine_forms import quadratic_ranges
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], EXP["n_ac"]
    seed = seed_for("ac_design_point", key)
    st = study(sysc["case"], threads=worker_threads())
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    rq = quadratic_ranges(so, Ws, N_SUR_RADII, seed + 200)
    bounds, sig = {}, {}
    for tg, W in Ws.items():
        Ct = np.asarray(W.T @ (st.C1 * delta))
        sig[tg] = sigma(Ct)
        lo_q, hi_q = rq[tg]
        for name, r in radii(Ct, ALPHA).items():
            bounds[(tg, name, "none")] = (-r, r)
            bounds[(tg, name, "second_order")] = (-r + np.minimum(0.0, lo_q), r + np.maximum(0.0, hi_q))
    counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
    anyc = {b: 0 for b in bounds}
    rng = np.random.default_rng(seed)
    n_ac = 0
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    for s in range(0, n, CFG["powerflow"]["batch"]):
        b = min(CFG["powerflow"]["batch"], n - s)
        eps = sample_eps(rng, st.src.n, b)
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
        n_ac += int(conv.sum())
        tv = {k: WT[k] @ ac for k in Ws}
        for bk, (lo, hi) in bounds.items():
            bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
            counts[bk] += bad.sum(axis=1)
            anyc[bk] += int(bad.any(axis=0).sum())
    rows = []
    for (tg, name, cor), c in counts.items():
        lo, hi = bounds[(tg, name, cor)]
        half = 0.5 * (hi - lo)
        ok = sig[tg] > 1e-12
        aff = bounds[(tg, "affine", "none")][1]
        g3 = bounds[(tg, "gauss_3sigma", "none")][1]
        rate = c / max(n_ac, 1)
        rows.append({"system": key, "delta": delta, "seed": seed, "n_ac": n_ac, "n_surrogate": N_SUR_RADII, "target": tg,
                     "radius": name, "correction": cor,
                     "halfwidth_over_sigma_median": float(np.median(half[ok] / sig[tg][ok])),
                     "width_over_3sigma_median": float(np.median(half[ok] / g3[ok])),
                     "width_over_affine_median": float(np.median(half[ok] / aff[ok])),
                     "capped_by_affine_pct": pct(np.mean(bounds[(tg, name, "none")][1][ok] >= aff[ok] * (1 - 1e-12))),
                     "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()), "n_over_alpha": int((rate > ALPHA).sum()),
                     "system_level_pct": pct(anyc[(tg, name, cor)] / max(n_ac, 1))})
    print(f"[b5] radii {key} done", flush=True)
    return rows


def main():
    import sys
    from common import grid
    from cdhim import loadability_limit
    if "radii" in sys.argv[1:]:
        rows = [r for o in parallel(run_radii, list(CFG["systems"]), workers=2) for r in o]
        print(save_table("b5_radii", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR_RADII,
                                            "note": "every radius capped by the affine radius; second_order adds the sampled range of q"}))
        return
    if "surrogate" in sys.argv[1:]:
        rows = [r for k in CFG["systems"] for r in run_surrogate(k)]
        print(save_table("b5_surrogate", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA,
                                                "note": "run alone: the t_* columns are the single-process cost of the method"}))
        return
    if "sweep250" in sys.argv[1:]:                                   # the sweep with 250 surrogate samples
        global N_SURROGATE
        N_SURROGATE = 250
        tasks = []
        for key, sysc in CFG["systems"].items():
            lam_max, _ = loadability_limit(grid(sysc["case"]), tol=CFG["powerflow"]["tol_pu"])
            levels = [(lf, "fixed") for lf in EXP["load_factors"]] + [(round(EXP["near_limit_fraction"] * lam_max, 4), "near_limit")]
            tasks += [(key, lf, tag, lam_max) for lf, tag in levels]
        rows = [r for o in parallel(run_sweep, tasks, workers=3) for r in o]
        print(save_table("b5_sweep_n250", rows, {"n_samples_per_cell": EXP["n_sweep"], "n_surrogate": 250, "alpha": ALPHA}))
        return
    if "sweep" in sys.argv[1:]:
        tasks = []
        for key, sysc in CFG["systems"].items():
            lam_max, _ = loadability_limit(grid(sysc["case"]), tol=CFG["powerflow"]["tol_pu"])
            levels = [(lf, "fixed") for lf in EXP["load_factors"]] + [(round(EXP["near_limit_fraction"] * lam_max, 4), "near_limit")]
            tasks += [(key, lf, tag, lam_max) for lf, tag in levels]
        rows = [r for o in parallel(run_sweep, tasks, workers=3) for r in o]
        print(save_table("b5_sweep", rows, {"n_samples_per_cell": EXP["n_sweep"], "n_surrogate": N_SURROGATE, "alpha": ALPHA}))
        return
    out = parallel(run, list(CFG["systems"]), workers=2)
    meta = {"n_samples": EXP["n_ac"], "n_surrogate": N_SURROGATE, "alpha": ALPHA,
            "note": "search row of b5_width: ratio_p10 = share of targets where the search flipped any symbol, "
                    "ratio_median = median model gain of the search over the affine radius, ratio_p90 = its maximum"}
    print(save_table("b5_mc", [r for o in out for r in o[0]], meta))
    print(save_table("b5_vertex", [r for o in out for r in o[1]], meta))
    print(save_table("b5_width", [r for o in out for r in o[2]], meta))


if __name__ == "__main__":
    main()
