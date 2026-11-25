"""Run CD-HIM on one system and print a summary.

    python -m cdhim ieee33
    python -m cdhim pegase2869 --radius coordinate --csv bounds.csv
    python -m cdhim case118 --delta 0.05 --legacy
"""
from __future__ import annotations

import argparse
import warnings
from time import perf_counter

import numpy as np

from .config import load_config, method_config
from .network import load_grid
from .study import make_study


def main(argv=None):
    cfg = load_config()
    ap = argparse.ArgumentParser(prog="python -m cdhim", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", help="a key of 'systems' in config.yaml, a pandapower case name, or a pglib_opf_* file name")
    ap.add_argument("--delta", type=float, help="relative half-width of every uncertain injection")
    ap.add_argument("--radius", choices=["lambda_max", "coordinate"])
    ap.add_argument("--gamma", type=float, help="gain threshold gamma_th (fraction)")
    ap.add_argument("--gain-box", choices=["covariance_1sigma", "covariance_3sigma", "affine"])
    ap.add_argument("--legacy", action="store_true", help="reproduce finalFinalCode.py (flat start, H sign error)")
    ap.add_argument("--csv", help="write the per-bus bounds to this file")
    a = ap.parse_args(argv)

    sysc = cfg["systems"].get(a.system, {"case": a.system, "delta": 0.05})
    delta = a.delta if a.delta is not None else sysc["delta"]
    over = {k: v for k, v in (("radius", a.radius), ("gamma_th", a.gamma), ("gain_box", a.gain_box)) if v is not None}
    lin = dict(point="flat", jacobian="legacy_signbug") if a.legacy else dict(cfg["linearization"])
    if a.legacy:
        over.update(min_sigma_product=1e-12, min_radius_pu=0.0)
        lin["permc_spec"] = cfg["linearization"]["permc_spec"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid = load_grid(sysc["case"], tol=cfg["powerflow"]["tol_pu"])
    t0 = perf_counter()
    st = make_study(grid, tol=cfg["powerflow"]["tol_pu"], min_injection_pu=cfg["uncertainty"]["min_injection_pu"], **lin)
    res = st.run(delta, method_config(cfg, **over))
    t = perf_counter() - t0

    e = res.bus_ellipse
    print(f"{grid.name}: {grid.n_bus} buses ({len(grid.pq)} PQ), N_u = {st.src.n} noise symbols, delta = {100 * delta:g} %")
    print(f"  linearisation      : {lin['point']} point, {lin['jacobian']} Jacobian, SuperLU ({st.lin.permc_spec}), "
          f"nnz(L+U)/nnz(J) = {(st.lin.nnz_L + st.lin.nnz_U) / st.lin.nnz_J:.2f}")
    print(f"  PQ-PQ branches     : {res.n_branches}, Ellipse mode: {res.n_ellipse} ({100 * res.n_ellipse / max(res.n_branches, 1):.1f} %)")
    print(f"  CRI (area, vs {res.cfg.gain_box} box, Ellipse branches): {100 * res.cri_area():.2f} %")
    print(f"  output width vs affine box (Ellipse buses, radius={res.cfg.radius}): {100 * res.width_reduction('affine'):+.2f} % reduction")
    print(f"  output width vs 3-sigma box                        : {100 * res.width_reduction('cltm'):+.2f} % reduction")
    print(f"  buses wider than the affine box: {int((res.r_out[e] > res.r_affine[e]).sum())} of {int(e.sum())}")
    print(f"  time (Jacobian + LU + solves + branches): {t:.3f} s")
    if a.csv:
        import pandas as pd

        v0 = np.abs(st.V_star[grid.pq])
        pd.DataFrame({"bus": grid.bus_label[grid.pq], "v0_pu": v0, "radius_pu": res.r_out, "v_min_pu": v0 - res.r_out,
                      "v_max_pu": v0 + res.r_out, "mode": np.where(e, "ellipse", "box"), "radius_affine_pu": res.r_affine,
                      "radius_3sigma_pu": res.r_cltm, "radius_coordinate_pu": res.r_coord}).to_csv(a.csv, index=False)
        print(f"  bounds written to {a.csv}")


if __name__ == "__main__":
    main()
