"""B6 -- Correlated sources in the affine-forms method (bus intervals and branch drops, design point).

Correlation is written the affine-arithmetic way, with shared independent symbols (cdhim.affine_forms.FactorModel):
    eps_k = w zeta_{z(k)} + (1 - w) eta_k,   zeta, eta independent U[-1, 1],
with one global factor or one factor per spectral zone (8 zones). The affine forms are rewritten on the independent
symbols, so the Bernstein radius applies exactly, and the range of the quadratic term is sampled under the same law.

Also tested: the one-factor Gaussian copula of C2, which is NOT a sum of independent bounded symbols. There the
Bernstein radius has no justification; three ways out are measured: pretend independence, a global factor model
matched to the copula's Pearson correlation, and the worst case of the linear model.

Every model with the slack bus alone (manuscript) and with the generators sharing the imbalance (participation).
Bounds are always centred on the converged AC solution and carry the sampled range of q (1.000 samples).

    python experiments/b6_correlated.py
Writes results/b6_correlated.csv
"""
from time import perf_counter

import numpy as np

from common import CFG, grid, method, parallel, pct, seed_for, study, worker_threads
from cdhim import SecondOrder, participation_factors, sample_eps, save_table, spectral_zones
from cdhim.affine_forms import (FactorModel, affine_radius, bernstein_radius, bus_targets, drop_targets, quadratic_ranges)
from cdhim.uncertainty import copula_to_pearson_uniform
from cdhim.validation import VIOLATION_ATOL

EXP = CFG["experiments"]
ALPHA = 1.0 - CFG["method"]["confidence"]
N_SUR = 1000
RHOS = (0.155, 0.5, 0.845)            # w = 0.3, 0.5, 0.7
COPULA_R = (0.3, 0.6, 0.9)
N_ZONES = 8


def models(st, key):
    """name -> (kind, sampler, coefficient map or None, descriptive fields)."""
    n = st.src.n
    zone_bus = spectral_zones(grid(CFG["systems"][key]["case"]), N_ZONES)
    zone_src = np.asarray(zone_bus)[st.src.bus]
    out = {"independent": ("factor", FactorModel(np.zeros(n, int), 0.0), {"rho": 0.0, "factor": "none"})}
    for rho in RHOS:
        w = FactorModel.w_for_rho(rho)
        out[f"global_rho{rho}"] = ("factor", FactorModel(np.zeros(n, int), w), {"rho": rho, "factor": "global"})
        out[f"zones_rho{rho}"] = ("factor", FactorModel(zone_src, w), {"rho": rho, "factor": f"{N_ZONES} zones"})
    for r in COPULA_R:
        out[f"copula_r{r}"] = ("copula", r, {"rho": float(copula_to_pearson_uniform(r)), "factor": "Gaussian copula", "copula_r": r})
    return out


def run(task):
    key, balancing = task
    sysc = CFG["systems"][key]
    delta, n = sysc["delta"], EXP["n_ac"]
    g = grid(sysc["case"])
    st = study(sysc["case"], threads=worker_threads(), participation=participation_factors(g) if balancing == "participation" else None)
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    Ct = {k: np.asarray(W.T @ (st.C1 * delta)) for k, W in Ws.items()}
    r_aff = {k: affine_radius(c) for k, c in Ct.items()}
    r_ind = {k: bernstein_radius(c, ALPHA) for k, c in Ct.items()}        # Bernstein if the sources were independent
    WT = {k: W.T.tocsr() for k, W in Ws.items()}
    rows = []
    for name, (kind, spec, info) in models(st, key).items():
        t0 = perf_counter()
        seed = seed_for("b6", key, balancing, name)
        if kind == "factor":
            fm = spec
            sampler = fm.sample
            radii = {k: {"bern_model": bernstein_radius(fm.coefficients(c), ALPHA),
                         "affine_model": affine_radius(fm.coefficients(c)),
                         "sig3_true": 3.0 * np.sqrt((fm.coefficients(c) ** 2).sum(axis=1) / 3.0)} for k, c in Ct.items()}
        else:
            r = spec
            sampler = lambda rng, m, r=r: sample_eps(rng, st.src.n, m, copula_r=r)
            rho_u = info["rho"]
            fm = FactorModel(np.zeros(st.src.n, int), FactorModel.w_for_rho(rho_u))
            radii = {}
            for k, c in Ct.items():
                s1 = c.sum(axis=1)
                var_true = (1.0 / 3.0) * ((1 - rho_u) * (c ** 2).sum(axis=1) + rho_u * s1 ** 2)
                radii[k] = {"bern_matched_factor": bernstein_radius(fm.coefficients(c), ALPHA),
                            "sig3_true": 3.0 * np.sqrt(var_true)}
        for k in Ct:
            radii[k]["bern_independent"] = r_ind[k]
            radii[k]["affine"] = r_aff[k]
        rq = quadratic_ranges(so, Ws, N_SUR, seed + 1, sampler=sampler)
        bounds = {(k, b): (-rr + np.minimum(0.0, rq[k][0]), rr + np.maximum(0.0, rq[k][1]))
                  for k in Ct for b, rr in radii[k].items()}
        counts = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
        counts_lin = {b: np.zeros(Ws[b[0]].shape[1], np.int64) for b in bounds}
        anyc = {b: 0 for b in bounds}
        rng = np.random.default_rng(seed)
        n_ac = 0
        C = st.C1 * delta
        for s in range(0, n, CFG["powerflow"]["batch"]):
            b = min(CFG["powerflow"]["batch"], n - s)
            eps = sampler(rng, b)
            V, conv, _ = st.solver.solve(st.injections(eps, delta))
            ac = np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None]
            lin = C @ eps
            n_ac += int(conv.sum())
            tv = {k: WT[k] @ ac for k in Ws}
            tl = {k: WT[k] @ lin for k in Ws}
            for bk, (lo, hi) in bounds.items():
                bad = (tv[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tv[bk[0]] > hi[:, None] + VIOLATION_ATOL)
                counts[bk] += bad.sum(axis=1)
                anyc[bk] += int(bad.any(axis=0).sum())
                counts_lin[bk] += ((tl[bk[0]] < lo[:, None] - VIOLATION_ATOL) | (tl[bk[0]] > hi[:, None] + VIOLATION_ATOL)).sum(axis=1)
        for (k, bname), c in counts.items():
            lo, hi = bounds[(k, bname)]
            rate = c / max(n_ac, 1)
            rate_lin = counts_lin[(k, bname)] / n
            ok = r_ind[k] > 1e-12
            rows.append({"system": key, "balancing": balancing, "model": name, **{"copula_r": "", **info}, "delta": delta,
                         "n_samples": n, "n_ac_not_converged": n - n_ac, "seed": seed, "n_surrogate": N_SUR, "target": k,
                         "bound": bname, "max_pct": pct(rate.max()), "mean_pct": pct(rate.mean()),
                         "n_over_alpha": int((rate > ALPHA).sum()), "system_level_pct": pct(anyc[(k, bname)] / max(n_ac, 1)),
                         "lin_max_pct": pct(rate_lin.max()), "lin_n_over_alpha": int((rate_lin > ALPHA).sum()),
                         "width_over_bern_independent_median": float(np.median((0.5 * (hi - lo))[ok] / r_ind[k][ok])),
                         "width_over_affine_median": float(np.median((0.5 * (hi - lo))[ok] / r_aff[k][ok]))})
        print(f"[b6] {key} {balancing} {name}: {perf_counter() - t0:.0f}s, not converged {n - n_ac}", flush=True)
    return rows


def main():
    tasks = [(k, b) for k in CFG["systems"] for b in ("slack", "participation")]
    rows = [r for o in parallel(run, tasks, workers=3) for r in o]
    print(save_table("b6_correlated", rows, {"n_samples": EXP["n_ac"], "alpha": ALPHA, "n_surrogate": N_SUR,
                                             "zones": f"{N_ZONES} spectral zones (cdhim.spectral_zones)",
                                             "note": "factor models: eps = w zeta_zone + (1-w) eta, support [-1,1]; the per-source variance "
                                                     "falls as w grows (1/6 at rho = 0.5). Copula: C2's one-factor Gaussian copula with U[-1,1] marginals"}))


if __name__ == "__main__":
    main()
