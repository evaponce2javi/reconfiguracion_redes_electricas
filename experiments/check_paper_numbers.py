"""Check that every number quoted in paper/main.tex is the one the CSV tables give.

    python experiments/check_paper_numbers.py

Each check recomputes a value from results/*.csv, formats it as the paper writes it, and looks for that exact string
in main.tex. Exit code 1 if any check fails. Run it after regenerating results or editing the text.
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

from common import ROOT

RES = ROOT / "results"
TEX = (ROOT / "paper" / "main.tex").read_text(encoding="utf-8")
P = "pglib_opf_case2869_pegase"


def csv(name):
    return pd.read_csv(RES / f"{name}.csv")


def checks():
    out = []

    def add(what, text):
        out.append((what, text))

    b8 = csv("b8_ppf_baselines")
    b8b = b8[b8.target == "bus"]
    b9 = csv("b9_branch_currents")
    tm = csv("b11_timing").groupby("case").median(numeric_only=True)
    c8 = csv("c8_timing").set_index("system")
    st = csv("fig1_quadratic_stats")
    q = csv("fig1_example_samples").q
    z = (q - q.mean()) / q.std()
    b4 = csv("b4_sweep")
    b4 = b4[(b4.status == "ok") & (b4.system == "pegase2869")]
    rad = csv("b5_radii")
    rp = rad[(rad.system == "pegase2869") & (rad.target == "bus") & (rad.correction == "second_order")].set_index("radius")
    b11 = csv("b11_certified")
    w = csv("b5_width")
    b6 = csv("b6_correlated")
    b6p = b6[(b6.system == "pegase2869") & (b6.target == "bus")]
    fac = b6p[b6p.model.str.match("global|zones")]
    cop = b6p[b6p.model.str.startswith("copula")]
    vx = csv("b5_vertex")
    vb = vx[(vx.system == "pegase2869") & (vx.target == "bus")]
    vi = vx[(vx.system == "ieee33") & (vx.target == "bus")]
    b7 = csv("b7_generality")
    g7 = b7[(b7.status == "ok") & (b7.load_level != "design") & (b7.level == "per_target") & (b7.radius == "bernstein")]
    b5s = csv("b5_sweep")
    b5s = b5s[b5s.status == "ok"]

    # abstract
    add("abstract: cumulant/point-estimate worst bus", f"up to {b8b[b8b.method.isin(['cumulant_cf', 'pem_cf'])].max_pct.max():.0f}\\% of")
    add("abstract: worst branch current", f"up to {b9[b9.method.isin(['cumulant_cf', 'pem_cf'])].max_pct.max():.0f}\\% outside branch-current")
    add("abstract: time with currents", f"about {tm.loc[P, 'total_with_currents_s']:.0f}\\,s on PEGASE-2869")
    # section 3
    for tg, n in (("bus", "buses"), ("drop", "branch drops")):
        s = st[st.target == tg]
        add(f"§3 median |skew| {tg}", f"{s.skewness.abs().median():.2f}")
        add(f"§3 median excess kurtosis {tg}", f"{s.excess_kurtosis.median():.2f}")
        add(f"§3 number of {tg} targets", f"{len(s)}")
    sb = st[st.target == "bus"]
    add("§3 share of buses with negative skew", f"{100 * (sb.skewness < 0).mean():.1f}\\%")
    add("§3 median minimum of buses", f"$-{-sb.min_over_sd.median():.1f}$")
    add("§3 example sampled maximum", f"{z.max():.1f} standard")
    add("§3 normal quantile", f"at {norm.ppf(1 - 0.0027 / 2):.1f}")
    add("§3 3sigma worst", f"{b4.ac_cltm_3sigma_none_max_pct.max():.1f}\\% violation")
    add("§3 3sigma+moments worst", f"(up to {b4.ac_cltm_3sigma_moment_max_pct.max():.1f}\\%)")
    add("§3 cells failing with moments (3.44 sigma)", f"fails in {int((b4.ac_cdhim_lmax_moment_n_over > 0).sum())} of the 25")
    # section 4
    add("§4 Hoeffding vs Bernstein", f"within {100 * (1 - rp.loc['bernstein', 'halfwidth_over_sigma_median'] / rp.loc['hoeffding', 'halfwidth_over_sigma_median']):.0f}\\%")
    c4 = b11[(b11.kind == "certified") & (b11.n_surrogate == 4000)]
    c1 = b11[(b11.kind == "certified") & (b11.n_surrogate == 1000)]
    add("§4 certified widening n=4000", f"widen by {100 * (c4.width_over_uncertified_median.median() - 1):.1f}\\%")
    add("§4 alpha_L for n=4000", f"$\\alpha_L = {100 * (0.0027 - 2 / 4001):.2f}\\%$")
    add("§4.4 time voltages", f"{tm.loc[P, 'total_voltages_s']:.1f}\\,s for the {int(tm.loc[P, 'n_bus'])} bus")
    add("§4.4 time with currents", f"{tm.loc[P, 'total_with_currents_s']:.1f}\\,s with the currents at the {int(tm.loc[P, 'n_current_ends'])} branch")
    add("§4.4 batched Monte Carlo", f"takes {c8.loc['pegase2869', 'ac_mc_batched_s']:.0f}\\,s")
    add("§4.4 pandapower loop", f"about {c8.loc['pegase2869', 'pandapower_loop_extrapolated_s']:.0f}\\,s")
    # section 5
    p8 = b8b[(b8b.case == P) & (b8b.load_level == "design")].set_index("method")
    add("§5.2 point-estimate AC solves", f"{int(p8.loc['pem_cf', 'n_ac_solves'])} solutions")
    add("§5.2 cumulant worst", f"{b8b[b8b.method == 'cumulant_cf'].max_pct.max():.1f}\\% (cumulants)")
    add("§5.2 point-estimate worst", f"{b8b[b8b.method == 'pem_cf'].max_pct.max():.1f}\\% (point estimate)")
    add("§5.2 buses over (cumulants, PEM)", f"{int(p8.loc['cumulant_cf', 'n_over_alpha'])} and {int(p8.loc['pem_cf', 'n_over_alpha'])} of its")
    add("§5.3 worst bus over grids", f"worst bus violation is {max(g7[g7.target == 'bus'].max_pct.max(), b5s.ac_bus_bern_max_pct.max()):.3f}\\%")
    wd = w[(w.target == "drop") & (w.relative_to == "bern") & (w.bound == "sub_bern")].set_index("system")
    add("§5.4 drop gain PEGASE", f"{wd.loc['pegase2869', 'ratio_median']:.1f} times")
    add("§5.4 drop gain IEEE-33", f"{wd.loc['ieee33', 'ratio_median']:.0f} times on IEEE-33")
    add("§5.4 current linear error", f"reaches {b9.max_err_linear_pu.max():.1f}\\,p.u.")
    add("§5.4 Bernstein without range", f"up to {b9[b9.method == 'bern_linear_only'].max_pct.max():.0f}\\% of the samples")
    add("§5.4 proposed worst branch end", f"at {b9[b9.method == 'bern'].max_pct.max():.2f}\\%")
    add("§5.4 worst case over Bernstein, PEGASE", f"{b9[(b9.case == P) & (b9.method == 'affine')].width_over_bern_median.median():.2f} times wider")
    add("§5.5 factor worst", f"(worst bus {fac[fac.bound == 'bern_model'].max_pct.max():.3f}\\%")
    ind = b6p[(b6p.model == "independent") & (b6p.bound == "bern_model")].width_over_affine_median.max()
    add("§5.5 independent width", f"from {ind:.2f} to")
    fw = fac[fac.bound == "bern_model"].width_over_affine_median
    add("§5.5 factor width range", f"between {fw.min():.2f} and {fw.max():.2f}")
    add("§5.5 ignoring correlation", f"up to {b6p[b6p.bound == 'bern_independent'].max_pct.max():.1f}\\%")
    add("§5.5 matched copula", f"(up to {cop[cop.bound == 'bern_matched_factor'].max_pct.max():.1f}\\%)")
    add("§5.6 certified n=1000", f"({100 * (c1.width_over_uncertified_median.median() - 1):.1f}\\% with $n = 1000$)")
    add("§5.6 worst certified", f"at {b11[b11.kind == 'certified'].max_pct.max():.2f}\\%")
    # section 6
    add("§6 vertex buses outside", f"{int(vb.det_s_n_outside.max())} of {int(vb.n_targets.iloc[0])} buses")
    add("§6 vertex excess over width", f"{100 * vb.det_s_max_excess_over_width.max():.0f}\\% of the bound")
    m, e = f"{vi.det_s_max_excess_pu.max():.0e}".split("e")
    add("§6 IEEE-33 vertex excess", f"${m}\\cdot10^{{{int(e)}}}$")
    sysmax = max(g7.groupby(["case", "load_factor", "delta"]).any_bus_or_drop_pct.max().max(), 0)
    add("§6 per-target system-level share (grid)", f"below {sysmax:.3f}\\%")
    return out


def main():
    bad = 0
    flat = " ".join(TEX.split())
    for what, text in checks():
        ok = " ".join(text.split()) in flat
        bad += not ok
        print(f"{'ok  ' if ok else 'FAIL'}  {what:45s} -> {text!r}")
    print(f"\n{bad} failing check(s)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
