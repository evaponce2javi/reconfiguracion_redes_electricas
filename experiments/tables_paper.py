"""LaTeX tables of the affine-forms paper (paper/tables/*.tex), from the CSV tables in results/.

    python experiments/tables_paper.py

tab_baselines.tex   worst bus violation of published methods and of the proposed bound (B8)
tab_coverage.tex    coverage of the proposed bound over every load x dispersion grid (B5 sweep, B7 grid)
tab_cost.tex        single-process cost of the certified method (B11 timing)
"""
import pandas as pd

from common import ROOT

RES = ROOT / "results"
OUT = ROOT / "paper" / "tables"
NAMES = {"case33bw": "IEEE-33", "case118": "IEEE-118", "case1354pegase": "PEGASE-1354", "pglib_opf_case2869_pegase": "PEGASE-2869",
         "ieee33": "IEEE-33", "pegase2869": "PEGASE-2869"}
ORDER = ["case33bw", "case118", "case1354pegase", "pglib_opf_case2869_pegase"]
CELLS = {"design": "design point", "stress_delta": "load 1.2, $\\delta = 30\\%$", "near_limit": "95\\% of limit, $\\delta = 15\\%$"}


def f(x, d=2):
    return f"{x:.{d}f}"


def write(name, body):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(body, encoding="utf-8")
    print(OUT / name)


def tab_baselines():
    d = pd.read_csv(RES / "b8_ppf_baselines.csv")
    d = d[d.target == "bus"]
    p = d.pivot_table(index=["case", "load_level"], columns="method", values="max_pct")
    w = d.pivot_table(index=["case", "load_level"], columns="method", values="width_over_bern_median")
    lines = [r"\begin{tabular}{llrrrrr}", r"\toprule",
             r" & & \multicolumn{4}{c}{worst bus violation (\%)} & width \\ \cmidrule(lr){3-6}",
             r"system & cell & cumulants & point est. & $3\sigma$+range & proposed & cum./prop. \\", r"\midrule"]
    for c in ORDER:
        for k, cell in enumerate(("design", "stress_delta", "near_limit")):
            r = p.loc[(c, cell)]
            lines.append(f"{NAMES[c] if k == 0 else ''} & {CELLS[cell]} & {f(r.cumulant_cf)} & {f(r.pem_cf)} & {f(r.sig3)} & "
                         f"\\textbf{{{f(r.bern, 3)}}} & {f(w.loc[(c, cell)].cumulant_cf)} \\\\")
        lines.append(r"\addlinespace" if c != ORDER[-1] else r"\bottomrule")
    lines.append(r"\end{tabular}")
    write("tab_baselines.tex", "\n".join(lines) + "\n")


def tab_coverage():
    rows = []
    s5 = pd.read_csv(RES / "b5_sweep.csv")
    s5 = s5[s5.status == "ok"]
    for sysname in ("ieee33", "pegase2869"):
        g = s5[s5.system == sysname]
        rows.append((NAMES[sysname], len(g), int(((g.ac_bus_bern_n_over > 0) | (g.ac_drop_bern_n_over > 0)).sum()),
                     g.ac_bus_bern_max_pct.max(), g.ac_drop_bern_max_pct.max(), g.width_bus_bern_over_affine_median.median(),
                     g.nonconvergence_pct.max()))
    b7 = pd.read_csv(RES / "b7_generality.csv")
    b7 = b7[(b7.status == "ok") & (b7.load_level != "design") & (b7.level == "per_target") & (b7.radius == "bernstein")]
    for c in ("case118", "case1354pegase"):
        g = b7[b7.case == c]
        gb, gd = g[g.target == "bus"], g[g.target == "drop"]
        over = int(((gb.set_index(["load_factor", "delta"]).n_over_alpha > 0) | (gd.set_index(["load_factor", "delta"]).n_over_alpha > 0)).sum())
        rows.append((NAMES[c], len(gb), over, gb.max_pct.max(), gd.max_pct.max(), gb.width_over_affine_median.median(),
                     100 * gb.n_ac_not_converged.max() / gb.n_samples.max()))
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"system & cells & cells over 0.27\% & worst bus (\%) & worst drop (\%) & width / worst case & AC not converged (\%) \\",
             r"\midrule"]
    for r in rows:
        lines.append(f"{r[0]} & {r[1]} & {r[2]} & {f(r[3], 3)} & {f(r[4], 3)} & {f(r[5])} & {f(r[6], 1)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_coverage.tex", "\n".join(lines) + "\n")


def tab_cost():
    d = pd.read_csv(RES / "b11_timing.csv").groupby("case").median(numeric_only=True)
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"system & buses & branch ends & first order & second-order model & total, voltages & total, with currents \\",
             r"\midrule"]
    for c in ORDER:
        r = d.loc[c]
        lines.append(f"{NAMES[c]} & {int(r.n_bus)} & {int(r.n_current_ends)} & {f(r.linear_study_s + r.voltage_forms_s)} & "
                     f"{f(r.second_order_model_s)} & {f(r.total_voltages_s)} & {f(r.total_with_currents_s)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write("tab_cost.tex", "\n".join(lines) + "\n")


if __name__ == "__main__":
    tab_baselines()
    tab_coverage()
    tab_cost()
