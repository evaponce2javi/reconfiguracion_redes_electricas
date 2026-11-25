"""B2 -- Conservatism of what the method actually returns (the per-bus intervals of Algorithm 1).

Deterministic (no sampling). Writes results/b2_conservatism.csv
"""
import numpy as np

from common import CFG, method, pct, study
from cdhim import save_table


def stats3(x):
    return {"mean": float(np.mean(x)), "median": float(np.median(x)), "max": float(np.max(x))}


def describe(key, model, st, delta, legacy):
    res = st.run(delta, method(legacy=legacy))
    res_c = st.run(delta, method(legacy=legacy, radius="coordinate"))
    e = res.bus_ellipse
    m = res.ellipse
    row = {"system": key, "model": model, "delta": delta, "n_pq_buses": len(e), "n_ellipse_buses": int(e.sum()),
           "n_branches": res.n_branches, "n_ellipse_branches": res.n_ellipse}
    # (1) over-conservatism of lambda_max against the exact projection
    oc_pairs = res.overconservatism()                                  # every (Ellipse branch, terminal)
    oc_bus = res.r_out[e] / res.r_coord[e]                             # after keeping the smallest radius per bus
    for name, x in (("overcons_pairs", oc_pairs), ("overcons_bus", oc_bus)):
        row.update({f"{name}_{k}": v for k, v in stats3(x).items()})
    # (2) width reduction of the output intervals (1 - r_out / r_ref), Ellipse-mode buses and all PQ buses
    for label, r in (("lmax", res), ("coord", res_c)):
        for ref in ("affine", "cltm"):
            row[f"width_red_{label}_vs_{ref}_ellipse_buses_pct"] = pct(r.width_reduction(ref, "ellipse"))
            row[f"width_red_{label}_vs_{ref}_all_buses_pct"] = pct(r.width_reduction(ref, "all"))
        row[f"area_red_{label}_vs_affine_pct"] = pct(r.output_area_reduction("affine"))
    # (3) where lambda_max is worse than the deterministic box
    row["n_buses_lmax_wider_than_affine"] = int((res.r_out[e] > res.r_affine[e]).sum())
    pair_r = np.concatenate([res.r_lmax[m], res.r_lmax[m]])
    pair_aff = np.concatenate([res.r_affine[res.fi[m]], res.r_affine[res.ti[m]]])
    row["n_pairs_lmax_wider_than_affine"] = int((pair_r > pair_aff).sum())
    row["n_pairs"] = int(pair_r.size)
    row["n_buses_coord_wider_than_affine"] = int((res_c.r_out[e] > res.r_affine[e]).sum())
    # (4) the manuscript's CRI for reference, and mean absolute widths
    row["cri_area_ellipse_branches_pct"] = pct(res.cri_area())
    row["mean_width_affine_pu"] = float(2 * res.r_affine[e].mean())
    row["mean_width_cltm_pu"] = float(2 * res.r_cltm[e].mean())
    row["mean_width_lmax_pu"] = float(2 * res.r_out[e].mean())
    row["mean_width_coord_pu"] = float(2 * res_c.r_out[e].mean())
    row["affine_over_3sigma_median"] = float(np.median(res.r_affine[e] / res.r_cltm[e]))
    return row


def main():
    rows = []
    for key, sysc in CFG["systems"].items():
        rows.append(describe(key, "legacy", study(sysc["case"], legacy=True), sysc["delta"], True))
        rows.append(describe(key, "corrected", study(sysc["case"]), sysc["delta"], False))
    print(save_table("b2_conservatism", rows, {"note": "deterministic; overcons = sqrt(lambda_max/Q_ii); "
                     "width reduction = mean over buses of 1 - r_out/r_ref; negative = wider than the reference"}))


if __name__ == "__main__":
    main()
