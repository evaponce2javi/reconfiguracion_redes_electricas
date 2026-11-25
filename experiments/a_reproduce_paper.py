"""Before/after check: do the manuscript figures come out of the refactored code?

``legacy`` columns run the new package in legacy mode (flat-start Jacobian with the H sign error,
absolute sx*sy filter, 1e-9 jitter in the Mahalanobis distance, original seeds and sample sizes) and
must equal what the unmodified ``finalFinalCode.py`` prints. ``corrected`` columns use the exact
Jacobian at the converged AC solution and no jitter.

Writes results/a_reproduction.csv
"""
import numpy as np
from scipy import stats

from common import CFG, method, pct, study
from cdhim import save_table

PAPER = {
    "ieee33": {"branches": 31, "ellipse": 16, "ellipse_pct": 51.6, "cri_pct": 86.92, "fig4_rho": 0.9990, "fig4_gain_pct": 96.45,
               "table1": {0.05: 0.00, 0.10: 0.06, 0.20: 0.08, 0.30: 0.10}, "beta_max_pct": 0.40},
    "pegase2869": {"branches": 3218, "ellipse": None, "ellipse_pct": 92.7, "cri_pct": 58.18,
                   "table1": {0.05: 0.06, 0.10: 0.08, 0.20: 0.12, 0.30: 0.12}, "beta_max_pct": 0.19, "pp_at_0.50": 0.73},
}
CHI2 = stats.chi2.ppf(0.9973, 2)


def legacy_layout(st, delta):
    """Coefficients laid out exactly like the original noise vector: [P at every non-slack bus ; Q at
    every PQ bus] in ascending bus order, with positive half-widths (sens * |P_d| * delta)."""
    g = st.grid
    ns, pq = np.sort(g.pvpq), np.sort(g.pq)
    cols = np.where(st.src.is_q, len(ns) + np.searchsorted(pq, st.src.bus), np.searchsorted(ns, st.src.bus))
    full = np.zeros((st.C1.shape[0], len(ns) + len(pq)))
    full[:, cols] = st.C1 * np.sign(st.src.unit) * delta
    return full


def d2_jitter(Q, points, jitter):
    return np.sum(points * np.linalg.solve(Q + jitter * np.eye(2), points), axis=0)


def top_branches(res, rho_cut, k=1):
    ok = np.flatnonzero(res.valid & (np.abs(res.rho) < rho_cut))
    return ok[np.argsort(-res.gain[ok], kind="stable")[:k]]


def experiments(st, res, delta, jitter):
    """Table I, Beta and P-P exactly as coded in finalFinalCode.py (seeds 99, 42, 101)."""
    out = {}
    b = top_branches(res, 0.999)[0]
    rows = [res.fi[b], res.ti[b]]
    for tol in (0.05, 0.10, 0.20, 0.30):
        C = legacy_layout(st, tol)[rows]
        eps = np.random.default_rng(99).uniform(-1, 1, size=(C.shape[1], 5000))
        out[f"table1_{tol:.2f}_pct"] = pct(np.mean(d2_jitter(C @ C.T / 3.0, C @ eps, jitter) > CHI2))
    Cd = legacy_layout(st, delta)
    rng = np.random.default_rng(42)
    beta = []
    for b3 in top_branches(res, 0.9999, 3):
        C = Cd[[res.fi[b3], res.ti[b3]]]
        eps = (rng.beta(2, 5, size=(C.shape[1], 10000)) - 0.285) * 4.0
        beta.append(pct(np.mean(d2_jitter(C @ C.T / 3.0, C @ eps, jitter) > CHI2)))
    out["beta_case1_pct"], out["beta_case2_pct"], out["beta_case3_pct"] = beta
    out["beta_max_pct"] = max(beta)
    out["beta_cases_in_box_mode"] = int(sum(not res.ellipse[b3] for b3 in top_branches(res, 0.9999, 3)))
    C = Cd[rows]
    eps = np.random.default_rng(101).uniform(-1, 1, size=(C.shape[1], 20000))
    d2 = d2_jitter(C @ C.T / 3.0, C @ eps, jitter)
    for p in (0.50, 0.90, 0.9973):
        out[f"pp_at_{p:.4g}"] = float(np.mean(d2 <= stats.chi2.ppf(p, 2)))
    return out


def summary(key, st, res, delta, jitter):
    g = st.grid
    b = top_branches(res, 0.999)[0]
    out = {"branches": res.n_branches, "ellipse": res.n_ellipse, "ellipse_pct": pct(res.n_ellipse / res.n_branches),
           "cri_pct": pct(res.cri_area()), "mean_gain_all_pct": pct(res.cri_area("all")),
           "min_gain_pct": pct(res.gain[res.valid].min()), "n_rho_ge_0.999": int((np.abs(res.rho[res.valid]) >= 0.999).sum()),
           "fig4_branch": f"{g.bus_label[g.pq[res.fi[b]]]}-{g.bus_label[g.pq[res.ti[b]]]}",
           "fig4_rho": float(res.rho[b]), "fig4_gain_pct": pct(res.gain[b])}
    out.update(experiments(st, res, delta, jitter))
    return out


def main():
    rows = []
    for key, sysc in CFG["systems"].items():
        delta = sysc["delta"]
        cols = {}
        st = study(sysc["case"], legacy=True)
        cols["legacy_mode"] = summary(key, st, st.run(delta, method(legacy=True)), delta, jitter=1e-9)
        st = study(sysc["case"], point="flat", jacobian="exact")
        cols["flat_exact"] = summary(key, st, st.run(delta, method()), delta, jitter=0.0)
        st = study(sysc["case"])
        cols["corrected"] = summary(key, st, st.run(delta, method()), delta, jitter=0.0)
        paper = dict(PAPER[key])
        for tol, v in paper.pop("table1").items():
            paper[f"table1_{tol:.2f}_pct"] = v
        for q in cols["legacy_mode"]:
            rows.append({"system": key, "quantity": q, "paper": paper.get(q, ""), **{c: v[q] for c, v in cols.items()}})
    path = save_table("a_reproduction", rows, {
        "note": "legacy_mode = new code emulating finalFinalCode.py (incl. 1e-9 jitter, seeds 99/42/101); "
                "flat_exact = sign error fixed, still flat start; corrected = exact Jacobian at the AC solution, no jitter",
        "seeds": {"table1": 99, "beta": 42, "pp": 101}, "n": {"table1": 5000, "beta": 10000, "pp": 20000}})
    print(path)


if __name__ == "__main__":
    main()
