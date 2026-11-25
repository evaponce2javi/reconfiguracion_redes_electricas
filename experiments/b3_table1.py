"""B3 -- Table I regenerated: ellipse miscoverage inside the LINEAR model, large N, fixed seed.

For every delta the same noise draws are used. Two columns per row:
  exact   Mahalanobis distance without regularisation (invariant to delta by construction)
  jitter  with the absolute 1e-9 diagonal regularisation of finalFinalCode.py:472 (not invariant)

Writes results/b3_table1.csv (the manuscript's branch) and results/b3_all_branches.csv
"""
import numpy as np

from common import CFG, method, pct, seed_for, study
from cdhim import geometry as geo
from cdhim import sample_eps, save_table

DELTAS = (0.05, 0.10, 0.20, 0.30)
BATCH = 2000


def paper_branch(res):
    ok = np.flatnonzero(res.valid & (np.abs(res.rho) < 0.999))
    return ok[np.argmax(res.gain[ok])]


def run(key, model, legacy):
    sysc = CFG["systems"][key]
    n = CFG["experiments"]["n_linear"]
    seed = seed_for("b3", key)
    design = sysc["delta"]
    st = study(sysc["case"], legacy=legacy)
    res = st.run(design, method(legacy=legacy))                     # branch selection at the design point
    b = paper_branch(res)
    m = res.ellipse
    fi, ti = res.fi[m], res.ti[m]
    chi2 = res.chi2
    pair = st.C1[[res.fi[b], res.ti[b]]]
    Q1 = pair @ pair.T / 3.0                                        # covariance of the paper branch at delta = 1
    cnt_exact = {d: 0 for d in DELTAS}
    cnt_jit = {d: 0 for d in DELTAS}
    cnt_all = np.zeros(int(m.sum()), np.int64)
    max_dev = 0.0
    rng = np.random.default_rng(seed)
    for s in range(0, n, BATCH):
        nb = min(BATCH, n - s)
        eps = sample_eps(rng, st.src.n, nb)
        x = st.C1 @ eps                                             # all PQ buses at delta = 1
        d2_all = geo.mahalanobis2(design * x[fi], design * x[ti], res.q11[m], res.q22[m], res.q12[m])
        cnt_all += (d2_all > chi2).sum(axis=1)
        p = x[[res.fi[b], res.ti[b]]]
        ref = None
        for d in DELTAS:
            pd_, Qd = p * d, Q1 * d * d
            d2 = np.sum(pd_ * np.linalg.solve(Qd, pd_), axis=0)
            ref = d2 if ref is None else ref
            max_dev = max(max_dev, float(np.max(np.abs(d2 - ref) / np.maximum(ref, 1e-300))))
            cnt_exact[d] += int((d2 > chi2).sum())
            cnt_jit[d] += int((np.sum(pd_ * np.linalg.solve(Qd + 1e-9 * np.eye(2), pd_), axis=0) > chi2).sum())
    g = st.grid
    name = f"{g.bus_label[g.pq[res.fi[b]]]}-{g.bus_label[g.pq[res.ti[b]]]}"
    rows = [{"system": key, "model": model, "branch": name, "rho": float(res.rho[b]), "delta": d, "n_samples": n, "seed": seed,
             "violation_exact_pct": pct(cnt_exact[d] / n), "violation_with_1e-9_jitter_pct": pct(cnt_jit[d] / n),
             "std_error_pct": pct(np.sqrt(max(cnt_exact[d], 1)) / n),
             "max_rel_change_of_d2_across_delta": max_dev} for d in DELTAS]
    r = cnt_all / n
    allrow = {"system": key, "model": model, "n_samples": n, "seed": seed, "n_ellipse_branches": int(m.sum()),
              "mean_pct": pct(r.mean()), "median_pct": pct(np.median(r)), "p99_pct": pct(np.quantile(r, 0.99)),
              "max_pct": pct(r.max()), "n_branches_over_0.27pct": int((r > 0.0027).sum()),
              "n_branches_over_0.27pct_by_3_std_errors": int((r > 0.0027 + 3 * np.sqrt(0.0027 / n)).sum())}
    return rows, allrow


def main():
    t1, allb = [], []
    for key in CFG["systems"]:
        for model, legacy in (("legacy", True), ("corrected", False)):
            a, b = run(key, model, legacy)
            t1 += a
            allb.append(b)
            print(f"[b3] {key} {model} done", flush=True)
    meta = {"n_samples": CFG["experiments"]["n_linear"], "distribution": "eps ~ U[-1,1] independent",
            "note": "samples propagated through the linear model only (dX = C eps); same draws for every delta"}
    print(save_table("b3_table1", t1, meta))
    print(save_table("b3_all_branches", allb, meta))


if __name__ == "__main__":
    main()
