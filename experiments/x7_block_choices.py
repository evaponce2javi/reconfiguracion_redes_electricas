"""X7 -- Contribution 3 under scrutiny: RAA-K against the standard relaxation and a hybrid.

Bounds of eps' H_i eps on [-1, 1]^K for a sample of buses of every case with K <= 16:

  af2       term by term (what interval evaluation of the second-order expansion gives)
  raak      the group's k-nomial rewriting, cliques of up to K variables
  spectral  all eigenvectors as blocks, same theorem
  hybrid_r  the r eigenvectors of largest |eigenvalue| as blocks, RAA-K on what is left (r = 1, 2, 4)
  sdp       Shor's semidefinite relaxation with the box constraints X_ii <= 1, |X_ij| <= 1
            (the standard tighter bound in the literature on quadratic forms over a box)

Reference: the range of the form over all 2^K vertices (results/x2_bus.csv). For negative semidefinite
forms that range is exact, [min over vertices, 0]; for the others it is an inner estimate.

Needs cvxpy with an SDP solver (not in requirements.txt: only this script uses it).
Writes results/x7_block_choices.csv
"""
from time import perf_counter

import cvxpy as cp
import numpy as np
import pandas as pd

from common import CFG, ROOT, grid, method, seed_for, study
from cdhim import SecondOrder, participation_factors, save_table, spectral_zones
from cdhim.knomial import block_bound, knomial_bound, knomial_rewrite
from cdhim.second_order import bound_monomial
from cdhim.zones import laplacian_embedding

N_SAMPLE = 300
NODE_BUDGET = CFG["experiments"]["raak_node_budget"]


def hybrid_bound(H, lam, vec, r):
    """r leading spectral blocks, then RAA-K on the remainder."""
    top = np.argsort(-np.abs(lam))[:r]
    c1, U1 = np.sign(lam[top]), (vec[:, top] * np.sqrt(np.abs(lam[top]))).T
    c2, U2, _, _ = knomial_rewrite(H - (U1.T * c1) @ U1, None, len(H), node_budget=NODE_BUDGET)
    return block_bound(H, np.concatenate([c1, c2]), np.vstack([U1, U2]))


class Shor:
    """max / min of tr(H X) over [[1, x'], [x, X]] >= 0, X_ii <= 1, |X_ij| <= 1; H is a parameter."""

    def __init__(self, k):
        self.H = cp.Parameter((k, k), symmetric=True)
        Z = cp.Variable((k + 1, k + 1), PSD=True)
        X = Z[1:, 1:]
        cons = [Z[0, 0] == 1, cp.diag(X) <= 1, X <= 1, X >= -1]
        self.hi = cp.Problem(cp.Maximize(cp.trace(self.H @ X)), cons)
        self.lo = cp.Problem(cp.Minimize(cp.trace(self.H @ X)), cons)

    def bound(self, H):
        scale = np.abs(H).max()
        self.H.value = H / scale
        out = []
        for prob in (self.lo, self.hi):
            prob.solve(solver=cp.CLARABEL)
            assert prob.status in ("optimal", "optimal_inaccurate"), prob.status
            out.append(prob.value * scale)
        return out[0], out[1]


def run_case(key, balancing, k, embedding, ref):
    sysc = CFG["systems"][key]
    g = grid(sysc["case"])
    st = study(sysc["case"], zones=spectral_zones(g, k, embedding), participation=participation_factors(g) if balancing == "participation" else None)
    delta = sysc["delta"]
    res = st.run(delta, method())
    H = SecondOrder(st, delta).forms()
    ref = ref.set_index("bus")
    label = np.array([int(g.bus_label[b]) for b in g.pq])
    active = np.flatnonzero(res.r_affine > 1e-6)
    seed = seed_for("x7", key, balancing, k)
    pick = np.sort(np.random.default_rng(seed).choice(active, min(N_SAMPLE, len(active)), replace=False))
    shor = Shor(st.src.n)
    rows = []
    names = ("af2", "raak", "spectral", "hybrid_1", "hybrid_2", "hybrid_4", "sdp")
    acc = {n: [] for n in names}
    t = {n: 0.0 for n in names}
    q_lo, q_hi, nsd = [], [], []
    for i in pick:
        lam, vec = np.linalg.eigh(H[i])
        fns = {"af2": lambda: bound_monomial(H[i]), "raak": lambda: knomial_bound(H[i], st.src.n, NODE_BUDGET)[:2],
               "spectral": lambda: block_bound(H[i], np.sign(lam), (vec * np.sqrt(np.abs(lam))).T),
               "hybrid_1": lambda: hybrid_bound(H[i], lam, vec, 1), "hybrid_2": lambda: hybrid_bound(H[i], lam, vec, 2),
               "hybrid_4": lambda: hybrid_bound(H[i], lam, vec, 4), "sdp": lambda: shor.bound(H[i])}
        r = ref.loc[label[i]]
        q_lo.append(min(r.vertex_q_min, 0.0))
        q_hi.append(max(r.vertex_q_max, 0.0))
        nsd.append(bool(r.h_negative_semidefinite))
        for n, fn in fns.items():
            t0 = perf_counter()
            lo, hi = fn()
            t[n] += perf_counter() - t0
            slack = 1e-6 * (q_hi[-1] - q_lo[-1])                 # the SDP solver is accurate to about 1e-8 relative
            assert lo <= q_lo[-1] + slack and hi >= q_hi[-1] - slack, f"{n} does not enclose an attained value"
            acc[n].append((lo, hi))
    q_lo, q_hi, nsd = np.array(q_lo), np.array(q_hi), np.array(nsd)
    for n in names:
        b = np.array(acc[n])
        w = (b[:, 1] - b[:, 0]) / (q_hi - q_lo)
        row = {"system": key, "balancing": balancing, "K": st.src.n, "bound": n, "n_forms": len(pick), "n_negative_semidefinite": int(nsd.sum()),
               "seed": seed, "width_over_vertex_range_median": float(np.median(w)), "width_over_vertex_range_max": float(w.max()),
               "ms_per_form": 1e3 * t[n] / len(pick)}
        if nsd.any():
            depth = -q_lo[nsd]
            row.update({"exact_width_ratio_median": float(np.median((b[nsd, 1] - b[nsd, 0]) / depth)),
                        "exact_width_ratio_max": float(((b[nsd, 1] - b[nsd, 0]) / depth).max()),
                        "lower_side_over_depth_median": float(np.median(-b[nsd, 0] / depth)),
                        "upper_side_over_depth_median": float(np.median(b[nsd, 1] / depth))})
        if (~nsd).any():
            ind = ~nsd
            row.update({"indefinite_width_ratio_median": float(np.median(w[ind])),
                        "indefinite_lower_side_ratio_median": float(np.median(b[ind, 0] / np.minimum(q_lo[ind], -1e-300))),
                        "indefinite_upper_side_ratio_median": float(np.median(b[ind, 1] / np.maximum(q_hi[ind], 1e-300)))})
        rows.append(row)
    print(f"[x7] {key} {balancing} K={st.src.n}: " + "  ".join(f"{r['bound']} {r['width_over_vertex_range_median']:.3f}" for r in rows), flush=True)
    return rows


def main():
    bus = pd.read_csv(ROOT / "results" / "x2_bus.csv")
    rows = []
    for key, counts in CFG["experiments"]["zone_counts"].items():
        emb = laplacian_embedding(grid(CFG["systems"][key]["case"]))
        delta = CFG["systems"][key]["delta"]
        for balancing in CFG["experiments"]["balancing"][key]:
            for k in [c for c in counts if c in (8, 16)]:
                ref = bus[(bus.system == key) & (bus.balancing == balancing) & (bus.K == k) & np.isclose(bus.delta, delta)]
                rows += run_case(key, balancing, k, emb, ref)
    print(save_table("x7_block_choices", rows, {"n_sample": N_SAMPLE, "raak_node_budget": NODE_BUDGET, "cvxpy": cp.__version__,
                                                "note": "reference = range of the form over all 2^K vertices (results/x2_bus.csv)"}))


if __name__ == "__main__":
    main()
