"""X3 -- Does the second-order model find the worst scenario of every bus?

For each bus the candidate worst vertices of the noise box are taken from the first-order model
(-sign(c_i) for the lowest voltage, +sign(c_i) for the highest) and from the second-order model
(the vertex that minimises / maximises c_i'x + x'H_i x, by 1-flip descent from the first-order one).
The exact AC power flow is solved at those vertices only, and the voltage reached is compared with
the true extreme over ALL 2^K vertices, which x2_few_sources.py computed by enumeration.

Needs results/x2_bus.csv (run x2 first). No random numbers are involved.

Writes results/x3_worst_vertex.csv
"""
import numpy as np
import pandas as pd

from common import CFG, ROOT, grid, method, pct, study
from cdhim import participation_factors, save_table, spectral_zones
from cdhim.zones import laplacian_embedding
from x1_aa_second_order import SecondOrder
from x2_few_sources import BATCH, K_ENUM, R_FLOOR, model_minimiser_vertices, quadratic_forms

TOL = 1e-9          # p.u.; a candidate "attains" the extreme if it is within the accuracy of the AC solver


def ac_at_own_vertex(st, X, rows, delta):
    """d|V_i| of the exact AC solution at the vertex X[i] proposed for bus i, for the buses in ``rows``."""
    uniq, inv = np.unique(X[rows], axis=0, return_inverse=True)
    dv = np.empty((X.shape[0], len(uniq)))
    for s in range(0, len(uniq), BATCH):
        V, conv, _ = st.solver.solve(st.injections(np.ascontiguousarray(uniq[s:s + BATCH].T), delta))
        assert conv.all(), "an AC solve did not converge"
        dv[:, s:s + BATCH] = np.abs(V[st.grid.pq]) - st.vm_star[:, None]
    return dv[rows, np.ravel(inv)], len(uniq)


def run_case(key, balancing, k, embedding, ref):
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    g = grid(sysc["case"])
    alpha = participation_factors(g) if balancing == "participation" else None
    st = study(sysc["case"], zones=spectral_zones(g, k, embedding), participation=alpha)
    res = st.run(delta, method())
    C, r = st.coefficients(delta), res.r_affine
    H = quadratic_forms(SecondOrder(st, delta))
    active = np.flatnonzero(r > R_FLOOR)
    first_lo = -np.sign(C)
    first_lo[first_lo == 0] = 1.0
    cand = {"lo_first": first_lo, "hi_first": -first_lo,
            "lo_second": model_minimiser_vertices(C, H), "hi_second": model_minimiser_vertices(-C, -H)}
    ref = ref.set_index("bus").loc[[int(g.bus_label[g.pq[i]]) for i in active]]
    true = {"lo": ref.vertex_dv_min.to_numpy(), "hi": ref.vertex_dv_max.to_numpy()}
    row = {"system": key, "balancing": balancing, "K": st.src.n, "delta": delta, "n_buses": len(active), "n_vertices_enumerated": 2 ** st.src.n}
    solves = 0
    for name, X in cand.items():
        side, order = name.split("_")
        dv, n_unique = ac_at_own_vertex(st, X, active, delta)
        solves = max(solves, n_unique)
        gap = dv - true[side] if side == "lo" else true[side] - dv            # >= 0: how far from the true extreme
        assert gap.min() > -1e-8, "a candidate beats the enumeration: the reference is not the true extreme"
        row[f"{name}_buses_attaining_pct"] = pct((gap <= TOL).mean())
        row[f"{name}_buses_missing"] = int((gap > TOL).sum())
        row[f"{name}_gap_max_pu"] = float(gap.max())
        row[f"{name}_gap_over_r_max_pct"] = pct((gap / r[active]).max())
        row[f"{name}_ac_solves"] = n_unique
        row[f"{name}_changed_vs_first_pct"] = pct((X[active] != cand[f"{side}_first"][active]).any(axis=1).mean())
    return row


def main():
    bus = pd.read_csv(ROOT / "results" / "x2_bus.csv")
    rows = []
    for key, counts in CFG["experiments"]["zone_counts"].items():
        emb = laplacian_embedding(grid(CFG["systems"][key]["case"]))
        delta = CFG["systems"][key]["delta"]
        for balancing in CFG["experiments"]["balancing"][key]:
            for k in [c for c in counts if c <= K_ENUM]:
                ref = bus[(bus.system == key) & (bus.balancing == balancing) & (bus.K == k) & np.isclose(bus.delta, delta)]
                rows.append(run_case(key, balancing, k, emb, ref))
                print("[x3] " + "  ".join(f"{n}={rows[-1][n]}" for n in ("system", "balancing", "K", "lo_first_buses_missing", "lo_second_buses_missing",
                                                                           "hi_first_buses_missing", "hi_second_buses_missing")), flush=True)
    print(save_table("x3_worst_vertex", rows, {"reference": "exact AC at all 2^K vertices, results/x2_bus.csv", "tolerance_pu": TOL,
                                               "note": "deterministic; no random numbers"}))


if __name__ == "__main__":
    main()
