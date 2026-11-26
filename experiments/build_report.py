"""Assemble RESULTS.md from the CSV tables in results/ and the commentary in results_text.py.

Every number in a table is read from a CSV; nothing is typed by hand. Run after experiments/run_all.py.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import CFG, ROOT

RES = ROOT / "results"
SYS = {"ieee33": "IEEE-33", "pegase2869": "PEGASE-2869"}


def num(x, spec=".2f"):
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    if isinstance(x, (bool, np.bool_)):
        return "sí" if x else "no"
    if spec == "d":
        return f"{int(round(float(x)))}"
    return format(float(x), spec).replace(".", ",")


def table(df, cols):
    """``cols``: list of (column, header, format). Returns a Markdown table."""
    head = "| " + " | ".join(h for _, h, _ in cols) + " |"
    sep = "|" + "|".join("---" if f is None else "---:" for _, _, f in cols) + "|"
    lines = [head, sep]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) if f is None else num(r[c], f) for c, _, f in cols) + " |")
    return "\n".join(lines)


def load(name):
    df = pd.read_csv(RES / f"{name}.csv")
    if "system" in df:
        df["system"] = df["system"].map(lambda s: SYS.get(s, s))
    return df


def meta(name):
    return json.loads((RES / f"{name}.meta.json").read_text(encoding="utf-8"))


def source(names, n=None, seed=None):
    names = [names] if isinstance(names, str) else names
    m = meta(names[0])
    bits = [f"Script: `experiments/{m['script']}`", "CSV: " + ", ".join(f"`results/{n_}.csv`" for n_ in names)]
    if n:
        bits.append(f"N = {n}")
    if seed:
        bits.append(f"semilla: {seed}")
    return "*" + " · ".join(bits) + "*"


def seeds(df):
    return ", ".join(f"{s} ({k})" for k, s in df.groupby("system")["seed"].first().items()) if "system" in df else str(df["seed"].iloc[0])


# ------------------------------------------------------------------------------------------- tables
def t_reproduction():
    df = load("a_reproduction")
    keep = ["branches", "ellipse", "ellipse_pct", "cri_pct", "n_rho_ge_0.999", "fig4_branch", "fig4_rho", "fig4_gain_pct",
            "table1_0.05_pct", "table1_0.10_pct", "table1_0.20_pct", "table1_0.30_pct", "beta_max_pct", "beta_cases_in_box_mode", "pp_at_0.5"]
    names = {"branches": "Ramas PQ-PQ analizadas", "ellipse": "Ramas en modo Elipse", "ellipse_pct": "Ramas en modo Elipse (%)",
             "cri_pct": "CRI: G medio en ramas Elipse (%)", "n_rho_ge_0.999": "Ramas con \\|ρ\\| ≥ 0,999 (a Caja)",
             "fig4_branch": "Rama de la Fig. 4", "fig4_rho": "ρ de esa rama", "fig4_gain_pct": "G de esa rama (%)",
             "table1_0.05_pct": "Tabla I, δ = 5 % (%)", "table1_0.10_pct": "Tabla I, δ = 10 % (%)",
             "table1_0.20_pct": "Tabla I, δ = 20 % (%)", "table1_0.30_pct": "Tabla I, δ = 30 % (%)",
             "beta_max_pct": "Beta: máximo de las 3 ramas (%)", "beta_cases_in_box_mode": "Beta: ramas probadas que están en modo Caja",
             "pp_at_0.5": "P-P: cobertura a confianza 0,50"}
    out = []
    for s in SYS.values():
        sub = df[(df.system == s) & df.quantity.isin(keep)].copy()
        sub["quantity"] = pd.Categorical(sub["quantity"], keep, ordered=True)
        sub = sub.sort_values("quantity")
        rows = []
        for _, r in sub.iterrows():
            def cell(v):
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    return "–" if (v is None or (isinstance(v, float) and np.isnan(v)) or v == "") else str(v)
                if np.isnan(v):
                    return "–"
                q = str(r.quantity)
                if q in ("branches", "ellipse", "n_rho_ge_0.999", "beta_cases_in_box_mode"):
                    return num(v, "d")
                return num(v, ".4f" if q.startswith(("fig4_rho", "pp_at")) else ".2f")
            rows.append({"q": names[r.quantity], "paper": cell(r.paper), "leg": cell(r.legacy_mode), "flat": cell(r.flat_exact), "cor": cell(r.corrected)})
        out.append(f"**{s}**\n\n" + table(pd.DataFrame(rows), [("q", "Cantidad", None), ("paper", "Paper", None),
                   ("leg", "Paquete nuevo, modo *legacy*", None), ("flat", "Signo corregido, arranque plano", None), ("cor", "Corregido", None)]))
    return "\n\n".join(out)


def t_b1_main():
    df = load("b1_bus_violation")
    d = df[(df.reference == "ac") & (df.bound == "cdhim_lmax") & (df.bus_set == "all_pq")]
    return table(d, [("system", "Sistema", None), ("model", "Linealización", None), ("delta", "δ", ".2f"), ("n_buses", "Barras PQ", "d"),
                     ("mean_pct", "Media (%)", ".4f"), ("p99_pct", "P99 (%)", ".3f"), ("max_pct", "Máx. (%)", ".3f"),
                     ("n_buses_over_0.27pct", "Barras > 0,27 %", "d"), ("system_level_pct", "Muestras con ≥ 1 violación (%)", ".2f"),
                     ("n_ac_not_converged", "AC sin converger", "d")])


def t_b1_bounds():
    df = load("b1_bus_violation")
    d = df[(df.reference == "ac") & (df.model == "corrected") & (df.bus_set != "box_buses")].copy()
    d["bus_set"] = d.bus_set.map({"all_pq": "todas las PQ", "ellipse_buses": "solo modo Elipse"})
    return table(d, [("system", "Sistema", None), ("bound", "Cota", None), ("bus_set", "Barras", None), ("n_buses", "n", "d"),
                     ("mean_pct", "Media (%)", ".4f"), ("p99_pct", "P99 (%)", ".3f"), ("max_pct", "Máx. (%)", ".3f"),
                     ("n_buses_over_0.27pct", "Barras > 0,27 %", "d")])


def t_b1_branch():
    d = load("b1_branch_violation").copy()
    d["reference"] = d.reference.map({"ac": "AC no lineal", "lin": "modelo lineal"})
    return table(d, [("system", "Sistema", None), ("model", "Linealización", None), ("reference", "Muestras", None),
                     ("n_ellipse_branches", "Ramas Elipse", "d"), ("mean_pct", "Media (%)", ".3f"), ("p99_pct", "P99 (%)", ".3f"),
                     ("max_pct", "Máx. (%)", ".3f"), ("n_branches_over_0.27pct", "Ramas > 0,27 %", "d"),
                     ("system_level_pct", "Muestras con ≥ 1 elipse violada (%)", ".2f")])


def t_b2():
    df = load("b2_conservatism")
    spec = [("n_ellipse_buses", "Barras en modo Elipse", "d"), ("n_pq_buses", "Barras PQ", "d"),
            ("overcons_pairs_mean", "√(λ_max/Q_ii) por (rama, extremo): media", ".3f"),
            ("overcons_pairs_median", "√(λ_max/Q_ii) por (rama, extremo): mediana", ".3f"),
            ("overcons_pairs_max", "√(λ_max/Q_ii) por (rama, extremo): máximo", ".2f"),
            ("overcons_bus_mean", "Ídem por barra (radio mínimo): media", ".3f"),
            ("overcons_bus_median", "Ídem por barra: mediana", ".3f"), ("overcons_bus_max", "Ídem por barra: máximo", ".2f"),
            ("width_red_lmax_vs_affine_ellipse_buses_pct", "Reducción de ancho, λ_max vs caja afín, barras Elipse (%)", ".2f"),
            ("width_red_lmax_vs_affine_all_buses_pct", "Reducción de ancho, λ_max vs caja afín, todas las PQ (%)", ".2f"),
            ("width_red_coord_vs_affine_ellipse_buses_pct", "Reducción de ancho, proyección exacta vs caja afín, barras Elipse (%)", ".2f"),
            ("width_red_coord_vs_affine_all_buses_pct", "Reducción de ancho, proyección exacta vs caja afín, todas las PQ (%)", ".2f"),
            ("width_red_lmax_vs_cltm_ellipse_buses_pct", "Reducción de ancho, λ_max vs caja 3σ, barras Elipse (%)", ".2f"),
            ("width_red_coord_vs_cltm_ellipse_buses_pct", "Reducción de ancho, proyección exacta vs caja 3σ, barras Elipse (%)", ".2f"),
            ("area_red_lmax_vs_affine_pct", "Reducción de área del rectángulo de salida vs caja afín, λ_max (%)", ".2f"),
            ("cri_area_ellipse_branches_pct", "CRI del paper (área de elipse vs caja de covarianza) (%)", ".2f"),
            ("n_buses_lmax_wider_than_affine", "Barras donde λ_max es más ancho que la caja afín", "d"),
            ("n_pairs_lmax_wider_than_affine", "(rama, extremo) donde λ_max es más ancho que la caja afín", "d"),
            ("n_pairs", "Total de (rama, extremo)", "d"),
            ("n_buses_coord_wider_than_affine", "Barras donde la proyección exacta es más ancha que la caja afín", "d"),
            ("affine_over_3sigma_median", "Caja afín / caja 3σ (mediana)", ".2f")]
    cols = [f"{r.system}, {r.model}" for _, r in df.iterrows()]
    rows = [{"q": h, **{c: num(df.iloc[i][k], f) for i, c in enumerate(cols)}} for k, h, f in spec]
    return table(pd.DataFrame(rows), [("q", "Cantidad", None)] + [(c, c, None) for c in cols])


def t_b3():
    d = load("b3_table1")
    return table(d, [("system", "Sistema", None), ("model", "Linealización", None), ("branch", "Rama", None), ("rho", "ρ", ".5f"),
                     ("delta", "δ", ".2f"), ("violation_exact_pct", "Violación exacta (%)", ".4f"),
                     ("violation_with_1e-9_jitter_pct", "Con regularización 1e-9 (%)", ".4f"), ("std_error_pct", "Error estándar (%)", ".4f"),
                     ("max_rel_change_of_d2_across_delta", "Cambio rel. máx. de d² entre δ", ".1e")])


def t_b3_all():
    d = load("b3_all_branches")
    return table(d, [("system", "Sistema", None), ("model", "Linealización", None), ("n_ellipse_branches", "Ramas Elipse", "d"),
                     ("mean_pct", "Media (%)", ".4f"), ("median_pct", "Mediana (%)", ".4f"), ("p99_pct", "P99 (%)", ".4f"),
                     ("max_pct", "Máx. (%)", ".4f"), ("n_branches_over_0.27pct", "Ramas > 0,27 %", "d")])


def t_c1(system):
    d = load("c1_ac_sweep")
    d = d[(d.system == system)].sort_values(["load_factor", "delta"])
    ok = d[d.status == "ok"]
    return table(ok, [("load_factor", "Carga", ".3f"), ("delta", "δ", ".2f"), ("nonconvergence_pct", "No conv. NR (%)", ".2f"),
                      ("a_err_max_pu", "(a) error máx. (p.u.)", ".2e"), ("a_err_p99_pu", "(a) error P99 (p.u.)", ".2e"),
                      ("a_err_over_radius_p99", "(a) error/radio P99", ".4f"), ("a_err_over_radius_max", "(a) error/radio máx.", ".3f"),
                      ("b_ellipse_lin_mean_pct", "(b) elipse lineal, media (%)", ".3f"), ("b_ellipse_lin_max_pct", "(b) máx. (%)", ".3f"),
                      ("c_bus_lmax_mean_pct", "(c) barra AC, media (%)", ".4f"), ("c_bus_lmax_max_pct", "(c) barra AC, máx. (%)", ".3f"),
                      ("c_bus_lmax_n_over", "(c) barras > 0,27 %", "d"), ("c_bus_lmax_system_pct", "(c) sistema (%)", ".2f"),
                      ("c_ellipse_ac_mean_pct", "(c) elipse AC, media (%)", ".3f"), ("c_ellipse_ac_max_pct", "(c) elipse AC, máx. (%)", ".2f"),
                      ("c_ellipse_ac_n_over", "(c) elipses > 0,27 %", "d")])


def t_c1_other(system):
    d = load("c1_ac_sweep")
    ok = d[(d.system == system) & (d.status == "ok")].sort_values(["load_factor", "delta"])
    return table(ok, [("load_factor", "Carga", ".3f"), ("delta", "δ", ".2f"),
                      ("c_bus_coord_mean_pct", "Proy. exacta: media (%)", ".4f"), ("c_bus_coord_max_pct", "máx. (%)", ".3f"), ("c_bus_coord_n_over", "> 0,27 %", "d"),
                      ("c_bus_cltm_mean_pct", "Caja 3σ: media (%)", ".4f"), ("c_bus_cltm_max_pct", "máx. (%)", ".3f"), ("c_bus_cltm_n_over", "> 0,27 %", "d"),
                      ("c_bus_affine_mean_pct", "Caja afín: media (%)", ".4f"), ("c_bus_affine_max_pct", "máx. (%)", ".3f"), ("c_bus_affine_n_over", "> 0,27 %", "d")])


def t_c1_first():
    d = load("c1_first_exceedance").fillna("")
    d["first"] = [f"carga {num(float(a), '.3f')}, δ {num(float(b), '.2f')} → {num(float(v), '.3f')} %" if a != "" else "nunca" for a, b, v in
                  zip(d.first_load_factor, d.first_delta, d.value_pct)]
    return table(d, [("system", "Sistema", None), ("criterion", "Criterio", None), ("n_cells_over_0.27pct", "Celdas > 0,27 %", "d"),
                     ("n_cells", "Celdas", "d"), ("first", "Primera celda que lo supera (orden: carga, luego δ)", None)])


def t_c2():
    d = load("c2_correlated_sources")
    return table(d, [("system", "Sistema", None), ("copula_r", "r cópula", ".1f"), ("pearson_rho_eps", "ρ entre ε", ".3f"),
                     ("covariance_model", "Q", None), ("ellipse_branches_pct", "Ramas Elipse (%)", ".1f"), ("cri_area_pct", "CRI área (%)", ".2f"),
                     ("width_red_lmax_vs_affine_pct", "Red. ancho λ_max vs afín (%)", ".2f"),
                     ("ellipse_lin_mean_pct", "Elipse, lineal: media (%)", ".3f"), ("ellipse_lin_max_pct", "máx. (%)", ".2f"),
                     ("ellipse_lin_n_over", "> 0,27 %", "d"), ("ellipse_ac_mean_pct", "Elipse, AC: media (%)", ".3f"),
                     ("bus_lmax_ac_mean_pct", "Barra λ_max, AC: media (%)", ".4f"), ("bus_lmax_ac_max_pct", "máx. (%)", ".2f"),
                     ("bus_lmax_ac_n_over", "> 0,27 %", "d"), ("bus_cltm_ac_max_pct", "Caja 3σ, AC: máx. (%)", ".2f"),
                     ("bus_affine_ac_max_pct", "Caja afín, AC: máx. (%)", ".3f")])


def t_c3_marginals():
    d = load("c3_marginals")
    return table(d, [("marginal", "Distribución", None), ("support_lo", "Soporte mín.", ".3f"), ("support_hi", "Soporte máx.", ".3f"),
                     ("mean", "Media", ".4f"), ("variance", "Varianza", ".4f"), ("inside_[-1,1]", "Dentro de [−1, 1]", ".0f"),
                     ("fraction_outside_[-1,1]_pct", "Masa fuera de [−1, 1] (%)", ".2f"), ("definition", "Definición", None)])


def t_c3(system):
    d = load("c3_nongaussian")
    d = d[d.system == system]
    return table(d, [("marginal", "Distribución", None), ("covariance_model", "Varianza en Q", None), ("variance_eps", "Var(ε) real", ".4f"),
                     ("ellipse_lin_mean_pct", "Elipse, lineal: media (%)", ".3f"), ("ellipse_lin_p99_pct", "P99 (%)", ".3f"),
                     ("ellipse_lin_max_pct", "máx. (%)", ".3f"), ("ellipse_lin_n_over", "ramas > 0,27 %", "d"),
                     ("bus_lmax_lin_max_pct", "Barra λ_max, lineal: máx. (%)", ".3f"), ("bus_coord_lin_max_pct", "Barra proy., lineal: máx. (%)", ".3f"),
                     ("ellipse_ac_mean_pct", "Elipse, AC: media (%)", ".3f"), ("bus_lmax_ac_max_pct", "Barra λ_max, AC: máx. (%)", ".3f")])


def t_c3_combo():
    d = load("c3_skewed_correlated")
    return table(d, [("system", "Sistema", None), ("covariance_model", "Q", None), ("pearson_rho_eps", "ρ entre ε", ".3f"),
                     ("cri_area_pct", "CRI área (%)", ".2f"), ("ellipse_lin_mean_pct", "Elipse, lineal: media (%)", ".3f"),
                     ("ellipse_lin_max_pct", "máx. (%)", ".2f"), ("ellipse_lin_n_over", "> 0,27 %", "d"),
                     ("ellipse_ac_mean_pct", "Elipse, AC: media (%)", ".3f"), ("bus_lmax_lin_max_pct", "Barra λ_max, lineal: máx. (%)", ".3f"),
                     ("bus_lmax_ac_max_pct", "Barra λ_max, AC: máx. (%)", ".3f"), ("bus_lmax_ac_n_over", "barras > 0,27 %", "d")])


def t_c4_kappa():
    d = load("c4_kappa")
    return table(d, [("system", "Sistema", None), ("model", "Linealización", None), ("gain_box", "Caja en G", None),
                     ("kappa_min", "κ mín.", ".4f"), ("kappa_median", "κ mediana", ".4f"), ("kappa_max", "κ máx.", ".4f"),
                     ("min_gain_pct", "G mínimo (%)", ".2f")])


def t_c4():
    d = load("c4_rho_rule")
    cols = [("system", "Sistema", None), ("model", "Linealización", None), ("gain_box", "Caja en G", None),
            ("rho_th_equivalent", "ρ_th equivalente", ".4f"), ("n_ellipse_gain_rule", "Elipse (regla G)", "d"),
            ("n_ellipse_rho_rule", "Elipse (regla ρ)", "d"), ("n_branches_changed", "Ramas que cambian", "d"),
            ("cri_area_gain_pct", "CRI regla G (%)", ".2f"), ("cri_area_rho_pct", "CRI regla ρ (%)", ".2f"),
            ("time_gain_rule_ms", "Tiempo G (ms)", ".2f"), ("time_rho_rule_ms", "Tiempo ρ (ms)", ".2f")]
    return table(d, cols)


def t_c4_effect():
    d = load("c4_rho_rule")
    if "bus_ac_max_gain_pct" not in d:
        return ""
    d = d[d.bus_ac_max_gain_pct.notna()]
    return table(d, [("system", "Sistema", None), ("gain_box", "Caja en G", None), ("n_branches_changed", "Ramas que cambian", "d"),
                     ("width_red_vs_affine_gain_pct", "Red. ancho regla G (%)", ".2f"), ("width_red_vs_affine_rho_pct", "Red. ancho regla ρ (%)", ".2f"),
                     ("bus_ac_max_gain_pct", "Barra AC máx., G (%)", ".3f"), ("bus_ac_max_rho_pct", "Barra AC máx., ρ (%)", ".3f"),
                     ("ellipse_ac_mean_gain_pct", "Elipse AC media, G (%)", ".3f"), ("ellipse_ac_mean_rho_pct", "Elipse AC media, ρ (%)", ".3f")])


def t_c5():
    df = load("c5_system_level")
    spec = [("m_ellipse_branches", "m (ramas en modo Elipse)", "d"), ("chi2_quantile", "Cuantil χ²", ".2f"),
            ("radius_factor_vs_uncorrected", "Factor sobre los radios", ".3f"),
            ("any_bus_lmax_lin_pct", "≥ 1 barra fuera, intervalo λ_max, modelo lineal (%)", ".2f"),
            ("any_bus_lmax_ac_pct", "≥ 1 barra fuera, intervalo λ_max, AC (%)", ".2f"),
            ("expected_bus_violations_lmax_ac", "Barras fuera por muestra (esperado), λ_max, AC", ".4f"),
            ("any_bus_coord_lin_pct", "≥ 1 barra fuera, proyección exacta, lineal (%)", ".2f"),
            ("any_bus_coord_ac_pct", "≥ 1 barra fuera, proyección exacta, AC (%)", ".2f"),
            ("expected_bus_violations_coord_ac", "Barras fuera por muestra, proyección exacta, AC", ".4f"),
            ("any_bus_cltm_ac_pct", "≥ 1 barra fuera, caja 3σ, AC (%)", ".2f"),
            ("expected_bus_violations_cltm_ac", "Barras fuera por muestra, caja 3σ, AC", ".3f"),
            ("any_bus_affine_ac_pct", "≥ 1 barra fuera, caja afín, AC (%)", ".2f"),
            ("any_ellipse_lin_pct", "≥ 1 elipse violada, modelo lineal (%)", ".2f"),
            ("any_ellipse_ac_pct", "≥ 1 elipse violada, AC (%)", ".2f"),
            ("expected_ellipse_violations_lin", "Elipses violadas por muestra, lineal", ".3f"),
            ("expected_ellipse_violations_ac", "Elipses violadas por muestra, AC", ".3f"),
            ("cri_area_cov_box_pct", "CRI de área vs caja de covarianza (%)", ".2f"),
            ("cri_area_vs_affine_box_pct", "CRI de área vs caja afín (%)", ".2f"),
            ("width_red_lmax_vs_affine_pct", "Reducción de ancho λ_max vs caja afín (%)", ".2f"),
            ("width_red_coord_vs_affine_pct", "Reducción de ancho proyección exacta vs caja afín (%)", ".2f"),
            ("n_buses_lmax_wider_than_affine", "Barras con λ_max más ancho que la caja afín", "d"),
            ("n_buses_coord_wider_than_affine", "Barras con proyección más ancha que la caja afín", "d")]
    names = {"per_branch_99.73": "99,73 % por rama", "bonferroni": "Bonferroni"}
    cols = [f"{r.system}, {names[r.variant]}" for _, r in df.iterrows()]
    rows = [{"q": h, **{c: num(df.iloc[i][k], f) for i, c in enumerate(cols)}} for k, h, f in spec]
    return table(pd.DataFrame(rows), [("q", "Cantidad", None)] + [(c, c, None) for c in cols])


def t_c6():
    d = load("c6_gamma_transfer")
    out = []
    for box in d.gain_box.unique():
        sub = d[d.gain_box == box]
        piv = sub.pivot(index="case", columns="gamma_th", values="ellipse_pct").loc[CFG["transfer_cases"]]
        cri = sub.pivot(index="case", columns="gamma_th", values="cri_area_pct").loc[CFG["transfer_cases"]]
        rows = []
        for case in piv.index:
            rows.append({"case": case, **{f"g{g}": f"{num(piv.loc[case, g], '.1f')} / {num(cri.loc[case, g], '.1f')}" for g in piv.columns}})
        out.append(f"Caja en G: `{box}` — ramas en modo Elipse (%) / CRI de área (%)\n\n" + table(pd.DataFrame(rows),
                   [("case", "Sistema", None)] + [(f"g{g}", f"γ_th = {int(round(100 * g))} %", None) for g in piv.columns]))
    return "\n\n".join(out)


def t_c6_opt():
    d = load("c6_optimum")
    d["bus_ac_max_pct_range"] = d["bus_ac_max_pct_range"].str.replace(".", ",").str.replace("-", " – ")
    return table(d, [("case", "Sistema", None), ("gain_box", "Caja en G", None), ("distinct_selections_over_sweep", "Selecciones distintas en el barrido", "d"),
                     ("ellipse_pct_min", "Elipse mín. (%)", ".1f"), ("ellipse_pct_max", "Elipse máx. (%)", ".1f"),
                     ("gamma_opt_cri_over_time", "γ óptimo (CRI/tiempo)", ".2f"), ("gamma_opt_cri_only", "γ óptimo (solo CRI)", ".2f"),
                     ("time_spread_over_sweep_pct", "Variación de tiempo en el barrido (%)", ".1f"),
                     ("bus_ac_max_pct_range", "Violación AC máx. por barra, rango (%)", None)])


def t_c7():
    d = load("c7_timing")
    for s in ("jacobian", "lu", "solves", "branch", "total"):
        d[s] = [f"{num(m, '.2f')} ± {num(sd, '.2f')}" for m, sd in zip(d[f"{s}_mean_ms"], d[f"{s}_std_ms"])]
    return table(d, [("case", "Sistema", None), ("n_bus", "Barras", "d"), ("N_u", "N_u", "d"), ("n_pq_pq_branches", "Ramas PQ-PQ", "d"),
                     ("jacobian", "Jacobiano (ms)", None), ("lu", "LU (ms)", None), ("solves", "Resoluciones (ms)", None),
                     ("branch", "Ramas (ms)", None), ("total", "Total (ms)", None), ("solves_share_of_total_pct", "Resoluciones / total (%)", ".1f"),
                     ("solves_default_threads_ms", "Resoluciones con hilos (ms)", ".2f"), ("base_pf_pandapower_ms", "Flujo base pandapower (ms)", ".0f"),
                     ("dense_block_MB", "Bloque denso de sensibilidades (MB)", ".1f")])


def t_c7_fill():
    d = load("c7_fill")
    return table(d, [("case", "Sistema", None), ("permc_spec", "Ordenamiento", None), ("nnz_J", "nnz(J)", "d"), ("nnz_L", "nnz(L)", "d"),
                     ("nnz_U", "nnz(U)", "d"), ("fill_ratio_(L+U)/J", "(L+U)/J", ".2f"), ("lu_time_ms", "LU (ms)", ".2f")])


def t_c7_fit():
    d = load("c7_fit")
    d["ci"] = [f"[{num(a, '.2f')}; {num(b, '.2f')}]" for a, b in zip(d.ci95_low, d.ci95_high)]
    return table(d, [("quantity", "Cantidad", None), ("systems", "Sistemas", None), ("n_points", "Puntos", "d"), ("exponent", "Exponente", ".2f"),
                     ("ci", "IC 95 %", None), ("r2", "R²", ".3f")])


def t_c8():
    d = load("c8_methods")
    return table(d, [("system", "Sistema", None), ("method", "Método", None), ("guarantee", "Tipo de garantía", None),
                     ("mean_width_pu", "Ancho medio (p.u.)", ".3e"), ("median_ratio_to_ac_99.73_range", "Ancho / rango AC 99,73 % (mediana)", ".2f"),
                     ("bus_ac_violation_mean_pct", "Violación AC media (%)", ".4f"), ("bus_ac_violation_max_pct", "máx. (%)", ".3f"),
                     ("n_buses_over_0.27pct", "Barras > 0,27 %", "d"), ("system_level_violation_pct", "Sistema (%)", ".2f"), ("time_s", "Tiempo (s)", ".4f")])


def t_c8_time():
    df = load("c8_timing")
    spec = [("linearisation_and_sensitivities_s", "Jacobiano + LU + resoluciones (s)", ".4f"), ("cdhim_branch_stage_s", "Etapa por ramas (s)", ".5f"),
            ("cdhim_total_s", "CD-HIM total, sin flujo base (s)", ".4f"), ("base_pf_pandapower_s", "Flujo base, pandapower (s)", ".3f"),
            ("ac_mc_batched_s", "Monte Carlo AC por lotes, N muestras (s)", ".1f"), ("ac_mc_batched_ms_per_sample", "ídem, ms por muestra", ".3f"),
            ("pandapower_ms_per_pf", "pandapower.runpp, ms por flujo (protocolo original)", ".1f"),
            ("pandapower_loop_extrapolated_s", "Lazo pandapower extrapolado a N muestras (s)", ".0f"),
            ("speedup_vs_batched_mc", "Aceleración frente al Monte Carlo por lotes", ".1f"),
            ("speedup_vs_pandapower_loop_extrapolated", "Aceleración frente al lazo pandapower (extrapolado)", ".0f")]
    cols = list(df.system)
    rows = [{"q": h, **{c: num(df.iloc[i][k], f) for i, c in enumerate(cols)}} for k, h, f in spec]
    return table(pd.DataFrame(rows), [("q", "Cantidad", None)] + [(c, c, None) for c in cols])


def _x1(spec):
    """Transposed table of results/x1_summary.csv: one row per quantity, one column per system."""
    df = load("x1_summary")
    cols = list(df.system)
    rows = [{"q": h, **{c: num(df.iloc[i].get(k, np.nan), f) for i, c in enumerate(cols)}} for k, h, f in spec]
    return table(pd.DataFrame(rows), [("q", "Cantidad", None)] + [(c, c, None) for c in cols])


def t_x1_form():
    return _x1([("n_symbols", "Símbolos de ruido N_u", "d"), ("n_buses_above_floor", "Barras PQ con radio afín > 10⁻⁶ p.u.", "d"),
                ("lin_err_over_r_affine_median_pct", "Error de linealización máx. / radio afín: mediana sobre barras (%)", ".2f"),
                ("lin_err_over_r_affine_max_pct", "ídem, peor barra (%)", ".2f"),
                ("lin_err_over_r_3sigma_median_pct", "Error de linealización máx. / radio 3σ: mediana (%)", ".2f"),
                ("lin_err_over_r_3sigma_max_pct", "ídem, peor barra (%)", ".1f"),
                ("err_reduction_by_quadratic_median", "Factor en que baja el error al sumar εᵀH_iε: mediana", ".0f"),
                ("err_reduction_by_quadratic_min", "ídem, mínimo", ".0f"),
                ("forms_n", "Barras con H_i ensamblada", "d"),
                ("forms_negative_share_of_spectrum_median_pct", "Parte negativa del espectro de H_i: mediana (%)", ".1f"),
                ("forms_true_rise_over_depth_max_pct", "Mayor valor positivo alcanzado / mayor valor negativo alcanzado: máximo (%)", ".2f"),
                ("forms_eig90_median", "Autovalores que juntan el 90 % del espectro: mediana", ".1f")])


def t_x1_vertex():
    return _x1([("vertex_lo_buses_outside_affine_box", "Vértice inferior: barras con la solución AC fuera de la caja afín", "d"),
                ("vertex_lo_excess_over_r_median_pct", "Exceso / radio afín: mediana (%)", ".2f"),
                ("vertex_lo_excess_over_r_p90_pct", "ídem, percentil 90 (%)", ".2f"),
                ("vertex_lo_excess_over_r_max_pct", "ídem, máximo (%)", ".2f"),
                ("vertex_lo_buses_excess_over_10pct", "Barras con exceso > 10 % del radio", "d"),
                ("vertex_lo_buses_excess_over_100pct", "Barras con exceso > 100 % del radio", "d"),
                ("vertex_lo_second_order_misses_median_pct_of_excess", "Parte del exceso que εᵀH_iε no explica: mediana (%)", ".1f"),
                ("vertex_lo_second_order_misses_p90_pct_of_excess", "ídem, percentil 90 (%)", ".1f"),
                ("vertex_hi_buses_outside_affine_box", "Vértice superior: barras con la solución AC fuera de la caja afín", "d"),
                ("vertex_hi_slack_over_r_median_pct", "Vértice superior: holgura / radio afín, mediana (%)", ".2f"),
                ("vertex_lo_not_converged", "Vértices inferiores sin solución AC", "d"),
                ("vertex_hi_not_converged", "Vértices superiores sin solución AC", "d")])


def t_x1_bounds():
    df = load("x1_summary")
    names = [("natural", "Producto afín estándar, término a término", "N·N_u"),
             ("binomial_terms", "Los mismos términos como cuadrados de binomios escalados", "N·N_u"),
             ("monomial", "H_i ensamblada, monomio a monomio (k = 1, AF2)", "N_u²"),
             ("greedy_k2", "H_i ensamblada, bloques binomiales con presupuesto diagonal (k = 2)", "N_u² log N_u"),
             ("knomial_1", "Teorema k-nomial: 1 bloque espectral + residuo monomial", "N_u² + 1 autovector"),
             ("knomial_8", "Teorema k-nomial: 8 bloques espectrales + residuo monomial", "N_u² + 8 autovectores"),
             ("spectral", "Teorema k-nomial: todos los bloques espectrales, sin residuo", "N_u³")]
    rows = []
    for key, label, cost in names:
        row = {"b": label, "cost": cost}
        for _, r in df.iterrows():
            row[f"{r.system} med"] = num(r[f"forms_{key}_over_inner_median"], ".2f")
            row[f"{r.system} max"] = num(r[f"forms_{key}_over_inner_max"], ".2f")
            row[f"{r.system} rise"] = num(r[f"forms_{key}_spurious_rise_over_depth_median_pct"], ".0f")
        rows.append(row)
    cols = [("b", "Cota de εᵀH_iε", None), ("cost", "Costo por barra", None)]
    for sname in df.system:
        cols += [(f"{sname} med", f"{sname}: ancho / alcanzado, mediana", None), (f"{sname} max", "máx.", None),
                 (f"{sname} rise", "lado positivo / profundidad (%)", None)]
    return table(pd.DataFrame(rows), cols)


def t_x1_width():
    return _x1([("natural_width_over_r_affine_median_pct", "Producto afín estándar: ancho de la cota / radio afín, mediana (%)", ".1f"),
                ("binomial_width_over_r_affine_median_pct", "Binomios escalados: ídem (%)", ".1f"),
                ("binomial_over_natural_width_median", "Binomios / producto estándar, mediana", ".3f"),
                ("natural_share_products_median_pct", "De la cota estándar: productos dV·dI (%)", ".1f"),
                ("natural_share_pv_squares_median_pct", "ídem: cuadrados en barras PV (%)", ".1f"),
                ("natural_share_curvature_median_pct", "ídem: curvatura de \\|V\\| (%)", ".1f"),
                ("angle_radius_median_rad", "Radio afín del ángulo, mediana (rad)", ".4f"),
                ("magnitude_radius_median_pu", "Radio afín de la magnitud, mediana (p.u.)", ".4f"),
                ("forms_worst_depth_over_r_affine_median_pct", "Profundidad alcanzada / radio afín, barras de mayor error: mediana (%)", ".0f"),
                ("forms_worst_depth_over_r_affine_max_pct", "ídem, máximo (%)", ".0f"),
                ("forms_random_depth_over_r_affine_median_pct", "Profundidad alcanzada / radio afín, barras al azar: mediana (%)", ".1f"),
                ("forms_random_depth_over_r_affine_max_pct", "ídem, máximo (%)", ".1f"),
                ("forms_all_depth_over_r_affine_median_pct", "Profundidad alcanzada / radio afín, todas las barras: mediana (%)", ".2f"),
                ("forms_all_depth_over_r_affine_max_pct", "ídem, máximo (%)", ".2f")])


def t_x1_prob():
    df = load("x1_summary")
    kinds = [("sigma3_plain", "Caja 3σ, centrada en V*"), ("sigma3_shift", "Caja 3σ, centro corrido en tr(H_i)/3"),
             ("coord_plain", "Proyección exacta √(χ²Q_ii), centrada en V*"), ("coord_shift", "Proyección exacta, centro corrido"),
             ("hoeffding_plain", "Radio de Hoeffding, centrado en V*"), ("affine_plain", "Caja afín, centrada en V*"),
             ("affine_shift", "Caja afín, centro corrido")]
    rows = []
    for _, r in df.iterrows():
        for key, label in kinds:
            rows.append({"system": r.system, "kind": label, "mean": r[f"viol_{key}_mean_pct"], "max": r[f"viol_{key}_max_pct"],
                         "n": r[f"viol_{key}_buses_over_nominal"], "sys": r[f"viol_{key}_system_pct"]})
    return table(pd.DataFrame(rows), [("system", "Sistema", None), ("kind", "Intervalo", None), ("mean", "Violación AC media (%)", ".4f"),
                                      ("max", "máx. (%)", ".3f"), ("n", "Barras > 0,27 %", "d"), ("sys", "Muestras con ≥ 1 barra fuera (%)", ".2f")])


def t_x1_dep():
    return _x1([("n_branches", "Ramas PQ-PQ", "d"),
                ("zonotope_gain_median_pct", "Área del zonotopo frente a la caja afín: reducción mediana (%)", ".1f"),
                ("zonotope_gain_min_pct", "ídem, mínima (%)", ".1f"),
                ("zonotope_over_ellipse_area_median", "Área del zonotopo / área de la elipse al 99,73 %: mediana", ".2f"),
                ("branches_zonotope_smaller_than_ellipse", "Ramas donde el zonotopo es más chico que la elipse", "d"),
                ("drop_reduction_median_pct", "Radio de V_i − V_j: reducción de la forma afín frente a la resta de intervalos, mediana (%)", ".1f"),
                ("drop_reduction_min_pct", "ídem, mínima (%)", ".1f"), ("drop_reduction_max_pct", "ídem, máxima (%)", ".1f")])


def t_x1_radii():
    return _x1([("hoeffding_over_affine_median", "Radio de Hoeffding / radio afín, mediana", ".3f"),
                ("hoeffding_over_3sigma", "Radio de Hoeffding / radio 3σ", ".3f"),
                ("shift_over_r_3sigma_median_pct", "Corrimiento del centro \\|tr(H_i)/3\\| / radio 3σ: mediana (%)", ".2f"),
                ("shift_over_r_3sigma_max_pct", "ídem, máximo (%)", ".2f"),
                ("shift_vs_mc_mean_max_abs_diff_over_r_3sigma_pct", "Corrimiento predicho frente a la media del Monte Carlo: diferencia máxima / radio 3σ (%)", ".4f"),
                ("cltm_buses_below_6_effective_terms_above_floor", "Barras con menos de 6 términos efectivos (regla M = 6 de la Ref. [10])", "d")])


def _x2(design=True):
    d = load("x2_summary")
    d["balancing"] = d.balancing.map({"slack": "slack única", "participation": "participación"})
    d = d[d.design_delta] if design else d[d.K == d.groupby(["system", "balancing"]).K.transform("min")]
    return d.sort_values(["system", "balancing", "K", "delta"], kind="stable")


_X2_CASE = [("system", "Sistema", None), ("balancing", "Balance", None), ("K", "K", "d")]


def t_x2_setting():
    return table(_x2(), _X2_CASE + [
        ("n_vertices", "Vértices con flujo AC", "d"), ("r_affine_median_pu", "Radio afín, mediana (p.u.)", ".4f"),
        ("max_abs_dv_at_vertices_pu", "Mayor \\|ΔV\\| en un vértice (p.u.)", ".4f"), ("angle_radius_median_rad", "Radio afín del ángulo, mediana (rad)", ".3f"),
        ("affine_over_3sigma_radius_median", "Radio afín / radio 3σ, mediana", ".2f"),
        ("sigma3_buses_wider_than_affine", "Barras donde 3σ es más ancho que la caja afín", "d"),
        ("affine_width_over_mc9973_median", "Caja afín / rango MC 99,73 %", ".2f"),
        ("sigma3_width_over_mc9973_median", "Caja 3σ / rango MC 99,73 %", ".2f"),
        ("sigma3_mc_system_pct", "Caja 3σ: muestras con ≥ 1 barra fuera (%)", ".1f")])


def t_x2_affine():
    return table(_x2(), _X2_CASE + [
        ("n_buses", "Barras", "d"), ("affine_buses_with_vertex_outside", "Barras con algún vértice fuera", "d"),
        ("affine_low_excess_over_r_median_pct", "Exceso / radio: mediana (%)", ".1f"), ("affine_low_excess_over_r_p90_pct", "P90 (%)", ".1f"),
        ("affine_low_excess_over_r_max_pct", "máx. (%)", ".1f"), ("affine_vertex_shortfall_max_pu", "Exceso máx. (p.u.)", ".4f"),
        ("affine_buses_shortfall_over_0.001pu", "Barras con exceso > 0,001 p.u.", "d"),
        ("affine_mc_buses_outside", "Barras con alguna muestra MC fuera", "d"), ("affine_mc_samples_outside_pct", "Muestras MC con ≥ 1 barra fuera (%)", ".2f")])


def t_x2_orders():
    return table(_x2(), _X2_CASE + [
        ("second_order_depth_over_r_median_pct", "2.º orden / radio: mediana (%)", ".1f"), ("second_order_depth_over_r_p90_pct", "P90 (%)", ".1f"),
        ("second_order_depth_over_r_max_pct", "máx. (%)", ".1f"), ("second_order_depth_max_pu", "2.º orden máx. (p.u.)", ".4f"),
        ("third_order_over_r_median_pct", "3.er orden / radio: mediana (%)", ".2f"), ("third_order_over_r_p90_pct", "P90 (%)", ".2f"),
        ("third_order_over_r_max_pct", "máx. (%)", ".2f"), ("third_order_max_pu", "3.er orden máx. (p.u.)", ".5f"),
        ("third_over_second_order_median_pct", "3.er / 2.º orden: mediana (%)", ".1f"), ("third_over_second_order_p90_pct", "P90 (%)", ".1f"),
        ("buses_h_negative_semidefinite_pct", "Barras con H_i semidefinida negativa (%)", ".0f")])


def t_x2_gate1():
    d = _x2()
    rows = []
    for _, r in d.iterrows():
        for key, label in (("monomial", "AF2 (monomio a monomio)"), ("raak_kn", "RAA-K (cliques de hasta K)"),
                           ("spectral", "Bloques espectrales")):
            rows.append({"system": r.system, "balancing": r.balancing, "K": r.K, "bound": label,
                         "out": r[f"{key}_buses_with_vertex_outside"], "short": r[f"{key}_vertex_shortfall_max_pu"],
                         "rel": r[f"{key}_vertex_shortfall_over_r_max_pct"], "n4": r[f"{key}_buses_shortfall_over_0.0001pu"],
                         "n3": r[f"{key}_buses_shortfall_over_0.001pu"], "mc": r[f"{key}_mc_buses_outside"],
                         "mcs": r[f"{key}_mc_samples_outside_pct"], "ok": bool(r[f"gate1_contains_all_{key}"])})
    return table(pd.DataFrame(rows), _X2_CASE + [
        ("bound", "Cota del resto", None), ("out", "Barras con algún vértice fuera", "d"), ("short", "Déficit máx. (p.u.)", ".5f"),
        ("rel", "Déficit máx. / radio (%)", ".2f"), ("n4", "Barras con déficit > 10⁻⁴ p.u.", "d"), ("n3", "> 10⁻³ p.u.", "d"),
        ("mc", "Barras con alguna muestra MC fuera", "d"), ("mcs", "Muestras MC con ≥ 1 barra fuera (%)", ".2f"), ("ok", "Criterio 1", None)])


def t_x2_width():
    return table(_x2(), _X2_CASE + [
        ("reference_width_over_affine_width_median", "Rango AC en vértices / caja afín", ".3f"),
        ("monomial_width_over_reference_median", "AF2 / rango AC: mediana", ".3f"), ("monomial_width_over_reference_max", "máx.", ".2f"),
        ("raak_kn_width_over_reference_median", "RAA-K / rango AC: mediana", ".3f"), ("raak_kn_width_over_reference_max", "máx.", ".2f"),
        ("spectral_width_over_reference_median", "Bloques espectrales / rango AC: mediana", ".3f"), ("spectral_width_over_reference_max", "máx.", ".2f"),
        ("natural_width_over_reference_median", "Producto estándar / rango AC: mediana", ".2f"),
        ("spectral_width_over_mc9973_median", "Bloques espectrales / rango MC 99,73 %", ".2f"),
        ("gate2_width_ok_spectral", "Criterio 2", None)])


def t_x2_operators():
    d = _x2()
    d = d[d.n_buses_exact_form_range > 0]
    names = [("natural", "Producto estándar"), ("monomial", "AF2"), ("raak_k2", "RAA-K, k = 2"), ("raak_k3", "RAA-K, k = 3"),
             ("raak_k6", "RAA-K, k = 6"), ("raak_kn", "RAA-K, k = K"), ("knomial_1", "Espectral, 1 bloque"),
             ("knomial_8", "Espectral, 8 bloques"), ("spectral", "Espectral, todos")]
    cols = _X2_CASE + [("n_buses_exact_form_range", "Barras con rango exacto", "d")]
    for key, label in names:
        d[f"{key}_cell"] = [f"{num(a, '.2f')} ({num(b, '.0f')} %)" for a, b in
                            zip(d[f"{key}_form_width_over_exact_median"], d[f"{key}_positive_side_over_depth_median_pct"])]
        cols.append((f"{key}_cell", label, None))
    return table(d, cols)


def t_x2_gain():
    return table(_x2(), _X2_CASE + [
        ("af2_over_raak_k2_form_width_median", "AF2 / RAA-K k = 2", ".2f"), ("af2_over_raak_k3_form_width_median", "k = 3", ".2f"),
        ("af2_over_raak_k6_form_width_median", "k = 6", ".2f"), ("af2_over_raak_kn_form_width_median", "k = K", ".2f"),
        ("af2_over_knomial_1_form_width_median", "AF2 / espectral 1 bloque", ".2f"),
        ("af2_over_knomial_8_form_width_median", "8 bloques", ".2f"), ("af2_over_spectral_form_width_median", "todos", ".2f"),
        ("monomial_remainder_over_r_median_pct", "Resto AF2 / radio, mediana (%)", ".1f"),
        ("first_order_seconds", "1.er orden (s)", ".3f"), ("second_order_seconds", "2.º orden (s)", ".2f"), ("bounds_seconds", "Cotas (s)", ".2f")])


def t_x2_sweep():
    d = _x2(design=False)
    d["delta_pct"] = 100 * d.delta
    return table(d, _X2_CASE + [
        ("delta_pct", "δ (%)", ".0f"), ("affine_low_excess_over_r_median_pct", "Exceso de la caja afín / radio: mediana (%)", ".2f"),
        ("affine_low_excess_over_r_max_pct", "máx. (%)", ".1f"), ("second_order_depth_over_r_max_pct", "2.º orden / radio, máx. (%)", ".1f"),
        ("third_order_over_r_max_pct", "3.er orden / radio, máx. (%)", ".2f"), ("third_order_max_pu", "3.er orden máx. (p.u.)", ".2e"),
        ("monomial_buses_with_vertex_outside", "AF2: barras con algún vértice fuera", "d"), ("monomial_vertex_shortfall_max_pu", "AF2: déficit máx. (p.u.)", ".2e"),
        ("spectral_buses_with_vertex_outside", "Espectral: barras fuera", "d"), ("spectral_vertex_shortfall_max_pu", "Espectral: déficit máx. (p.u.)", ".2e"),
        ("spectral_buses_shortfall_over_0.0001pu", "Espectral: barras con déficit > 10⁻⁴ p.u.", "d")])


def t_x3():
    d = load("x3_worst_vertex")
    d["balancing"] = d.balancing.map({"slack": "slack única", "participation": "participación"})
    return table(d, _X2_CASE + [
        ("n_buses", "Barras", "d"), ("n_vertices_enumerated", "Vértices de la enumeración", "d"),
        ("lo_first_buses_attaining_pct", "Voltaje mínimo, 1.er orden: barras que aciertan (%)", ".1f"), ("lo_first_gap_max_pu", "brecha máx. (p.u.)", ".5f"),
        ("lo_second_buses_attaining_pct", "Voltaje mínimo, 2.º orden: aciertan (%)", ".1f"), ("lo_second_gap_max_pu", "brecha máx. (p.u.)", ".5f"),
        ("lo_second_ac_solves", "flujos AC", "d"),
        ("hi_first_buses_attaining_pct", "Voltaje máximo, 1.er orden: aciertan (%)", ".1f"), ("hi_first_gap_max_pu", "brecha máx. (p.u.)", ".5f"),
        ("hi_second_buses_attaining_pct", "Voltaje máximo, 2.º orden: aciertan (%)", ".1f"), ("hi_second_gap_max_pu", "brecha máx. (p.u.)", ".5f"),
        ("hi_second_ac_solves", "flujos AC", "d")])


def x3_values():
    x = pd.read_csv(RES / "x3_worst_vertex.csv")
    peg = x[x.system == "pegase2869"]
    out = {}
    for tag, col, spec in (("LO1", "lo_first_buses_attaining_pct", ".0f"), ("LO2", "lo_second_buses_attaining_pct", ".0f"),
                           ("HI1", "hi_first_buses_attaining_pct", ".0f"), ("HI2", "hi_second_buses_attaining_pct", ".0f"),
                           ("SOLVES", "lo_second_ac_solves", "d"), ("GAP_LO2", "lo_second_gap_max_pu", ".4f"),
                           ("GAP_LO1", "lo_first_gap_max_pu", ".4f")):
        out[f"X3_{tag}_LO"], out[f"X3_{tag}_HI"] = num(peg[col].min(), spec), num(peg[col].max(), spec)
    return out


_BAL = {"slack": "slack única", "participation": "participación"}
_CASE = {"case33bw": "IEEE-33", "case118": "IEEE-118", "case300": "IEEE-300", "pglib_opf_case2869_pegase": "PEGASE-2869"}


def _sci(v):
    """1.23e-05 -> 1,2·10⁻⁵ (plain decimals above 1e-3)."""
    v = float(v)
    if v == 0 or abs(v) >= 1e-3:
        return num(v, ".4f")
    m, e = f"{v:.1e}".split("e")
    sup = str(int(e)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m.replace('.', ',')}·10{sup}"


def t_x4_independent():
    d = load("x4_independent")
    d["balancing"] = d.balancing.map(_BAL)
    spec = [("bus", "Barra", "d"), ("batch_solver_max_mismatch_pu", "Desajuste máximo del solver por lotes (p.u.)", "sci"),
            ("max_abs_vm_difference_pandapower_pu", "Mayor diferencia de \\|V\\| con pandapower, todas las barras (p.u.)", "sci"),
            ("r_affine", "Radio afín de la barra (p.u.)", ".4f"), ("dv_ac_batch", "ΔV, solver por lotes (p.u.)", ".6f"),
            ("dv_ac_pandapower", "ΔV, pandapower (p.u.)", ".6f"), ("affine_box_violated_by_pu_pandapower", "Caja afín: violada por (p.u.)", "sci"),
            ("af2_violated_by_pu_pandapower", "2.º orden con AF2: violada por (p.u.)", "sci"),
            ("raak_violated_by_pu_pandapower", "2.º orden con RAA-K: violada por (p.u.)", "sci"),
            ("spectral_violated_by_pu_pandapower", "2.º orden con bloques espectrales: violada por (p.u.; negativo = contiene)", "sci")]
    cols = [f"{r.system}, {r.balancing}, K = {r.K}" for r in d.itertuples()]
    rows = [{"q": h, **{c: (_sci(d.iloc[i][k]) if f == "sci" else num(d.iloc[i][k], f)) for i, c in enumerate(cols)}} for k, h, f in spec]
    return table(pd.DataFrame(rows), [("q", "Cantidad", None)] + [(c, c, None) for c in cols])


def t_x4_interior():
    d = load("x4_interior")
    d["balancing"] = d.balancing.map(_BAL)
    d["side"] = d.side.map({"lowest": "mínimo", "highest": "máximo"})
    for c in ("gain_over_vertices_max_pu", "polish_gain_max_pu"):
        d[c] = d[c].map(_sci)
    return table(d, [("system", "Sistema", None), ("balancing", "Balance", None), ("K", "K", "d"), ("side", "Voltaje", None), ("n_buses", "Barras", "d"),
                     ("buses_model_optimum_interior", "Óptimo del modelo en el interior", "d"),
                     ("buses_beating_all_vertices", "Barras que superan a todos los vértices", "d"), ("gain_over_vertices_max_pu", "por hasta (p.u.)", None),
                     ("buses_gain_over_1e-5pu", "> 10⁻⁵ p.u.", "d"), ("buses_gain_over_1e-4pu", "> 10⁻⁴ p.u.", "d"),
                     ("polish_gain_max_pu", "Tras afinar con el flujo AC exacto (p.u.)", None)])


def t_x5_table1():
    d = load("x5_liao_table1")
    iv = lambda a, b: [f"[{num(x, '.4f')}; {num(y, '.4f')}]" for x, y in zip(d[a], d[b])]
    out = pd.DataFrame({"bus": d.bus, "nom": d.nominal_here, "paper": iv("liao_paper_lower", "liao_paper_upper"),
                        "here": iv("liao_here_lower", "liao_here_upper"), "af2": iv("direct_af2_lower", "direct_af2_upper"),
                        "raak": iv("direct_raak_lower", "direct_raak_upper"), "exact": iv("exact_all_vertices_lower", "exact_all_vertices_upper"),
                        "mcp": iv("mc_paper_lower", "mc_paper_upper"), "mch": iv("mc_here_lower", "mc_here_upper")})
    return table(out, [("bus", "Barra", "d"), ("nom", "Nominal", ".4f"), ("paper", "Liao, Tabla 1 del artículo", None), ("here", "Liao reimplementado", None),
                       ("af2", "Forma directa, AF2", None), ("raak", "Forma directa, RAA-K", None), ("exact", "Exacto en los 256 vértices", None),
                       ("mcp", "Monte Carlo del artículo", None), ("mch", "Monte Carlo aquí (extremos de 20.000)", None)])


def t_x5_variants():
    d = load("x5_liao_variants")
    names = {"generator P, load P and Q": "P de generadores PV; P y Q de cargas", "load P and Q": "P y Q de cargas",
             "generator P, load P": "P de generadores PV; P de cargas", "load P": "P de cargas", "generator P": "P de generadores PV"}
    d["uncertain_variables"] = d.uncertain_variables.map(names)
    return table(d, [("uncertain_variables", "Variables inciertas", None), ("K", "K", "d"),
                     ("max_abs_diff_exact_range_vs_paper_mc", "Rango exacto aquí frente al Monte Carlo del artículo: mayor diferencia (p.u.)", ".4f"),
                     ("max_abs_diff_liao_here_vs_paper", "Liao reimplementado frente a su Tabla 1 (p.u.)", ".4f"),
                     ("paper_interval_contains_exact_range_buses", "Barras (de 6) donde el intervalo del artículo contiene el rango exacto", "d")])


def t_x5_summary():
    d = pd.read_csv(RES / "x5_liao_summary.csv")
    d["case"] = ["IEEE 9, ±10 %", "IEEE 57, ±20 %", "PEGASE-2869, 8 zonas, slack única", "PEGASE-2869, 8 zonas, participación"]
    d["state"] = [f"{num(a, '.2f')} ({num(b, 'd')} de {num(c, 'd')})" if a == a else "–" for a, b, c in
                  zip(d.liao_state_width_over_reference_median, d.liao_state_components_not_containing_reference, d.liao_state_components)]
    return table(d, [("case", "Caso", None), ("K", "K", "d"), ("angle_radius_median_rad", "Radio afín del ángulo (rad)", ".3f"),
                     ("state", "Liao sobre e y f: ancho / rango (componentes que no contienen)", None),
                     ("liao_width_over_reference_median", "\\|V\\| desde [e], [f] con intervalos: ancho / rango", ".2f"),
                     ("af2_width_over_reference_median", "Forma directa, AF2", ".3f"), ("raak_width_over_reference_median", "Forma directa, RAA-K", ".3f"),
                     ("linear_solves_liao", "Resoluciones lineales, vía de Liao", "d"), ("linear_solves_direct", "vía adjunta", "d")])


def t_x6():
    d = pd.read_csv(RES / "x6_fixed_point.csv")
    d["case"] = d.case.map(_CASE)
    d["balancing"] = d.balancing.map(_BAL)
    d["pct"] = 100 * d.delta
    d["fp"] = d.fixed_point_found.map({True: "sí", False: "no"})
    for c in ("allowance_over_r_median_pct", "raak_enclosure_width_over_reference_median", "raak_enclosure_buses_violated", "raak_second_order_only_buses_violated"):
        if c not in d:
            d[c] = np.nan
    return table(d, [("case", "Sistema", None), ("balancing", "Balance", None), ("pct", "δ (%)", ".0f"),
                     ("spectral_radius_G", "Radio de G, rectangulares", ".3f"), ("spectral_radius_G_polar_first_order", "Radio de G, polares (1.er orden)", ".3f"),
                     ("spectral_radius_true_max_over_64_vertices", "Variación real del Jacobiano, rectangulares", ".3f"),
                     ("spectral_radius_true_polar_max_over_64_vertices", "ídem, polares", ".3f"), ("fp", "Punto fijo", None),
                     ("allowance_over_r_median_pct", "Holgura añadida / radio, mediana (%)", ".3f"),
                     ("raak_enclosure_width_over_reference_median", "Encierro / rango exacto, mediana", ".3f"),
                     ("raak_enclosure_buses_violated", "Barras con algún vértice fuera", "d"),
                     ("raak_second_order_only_buses_violated", "ídem, solo 2.º orden", "d")])


def t_x7():
    d = load("x7_block_choices")
    d["balancing"] = d.balancing.map(_BAL)
    d["bound"] = d.bound.map({"af2": "AF2", "raak": "RAA-K", "spectral": "Bloques espectrales", "hybrid_1": "Híbrida, 1 autovector + RAA-K",
                              "hybrid_2": "Híbrida, 2 autovectores", "hybrid_4": "Híbrida, 4 autovectores", "sdp": "Relajación semidefinida"})
    return table(d, [("system", "Sistema", None), ("balancing", "Balance", None), ("K", "K", "d"), ("bound", "Cota", None), ("n_forms", "Formas", "d"),
                     ("width_over_vertex_range_median", "Ancho / rango en vértices: mediana", ".3f"), ("width_over_vertex_range_max", "máx.", ".2f"),
                     ("lower_side_over_depth_median", "Lado inferior / profundidad (formas semidefinidas negativas)", ".3f"),
                     ("upper_side_over_depth_median", "Lado superior / profundidad", ".3f"), ("ms_per_form", "ms por forma", ".2f")])


def x4_values():
    a, b = pd.read_csv(RES / "x4_independent.csv"), pd.read_csv(RES / "x4_interior.csv")
    s5, v5 = pd.read_csv(RES / "x5_liao_summary.csv"), pd.read_csv(RES / "x5_liao_variants.csv")
    x6, x7 = pd.read_csv(RES / "x6_fixed_point.csv"), pd.read_csv(RES / "x7_block_choices.csv")
    lo, hi = b[b.side == "lowest"], b[(b.side == "highest") & (b.system == "pegase2869")]
    r = lambda fr, col, spec: (num(fr[col].min(), spec), num(fr[col].max(), spec))
    out = {"X4_PP_MAX": _sci(a.max_abs_vm_difference_pandapower_pu.max()), "X4_LO_GAIN_MAX": _sci(lo.gain_over_vertices_max_pu.max()),
           "X4_LO_N_MAX": num(lo.buses_beating_all_vertices.max(), "d"), "X4_HI_GAIN_MAX": _sci(hi.gain_over_vertices_max_pu.max()),
           "X4_HI_N4_MAX": num(hi["buses_gain_over_1e-4pu"].max(), "d"), "X4_POLISH_EXTRA": _sci(b.polish_extra_over_model_point_max_pu.max()),
           "X5_VAR_BEST": num(v5.max_abs_diff_exact_range_vs_paper_mc.min(), ".3f"),
           "X5_STATE_OUT": num(s5.liao_state_components_not_containing_reference.iloc[0], "d"),
           "X5_STATE_N": num(s5.liao_state_components.iloc[0], "d"), "X5_STATE_SHORT": _sci(s5.liao_state_shortfall_max_pu.iloc[0]),
           "X5_LIAO_9": num(s5.liao_width_over_reference_median.iloc[0], ".1f"), "X5_LIAO_57": num(s5.liao_width_over_reference_median.iloc[1], ".1f"),
           "X5_SOLVES_L8": num(s5.linear_solves_liao.iloc[2], "d"), "X5_SOLVES_D8": f"{int(s5.linear_solves_direct.iloc[2]):,}".replace(",", "."),
           "X6_N_FOUND": num(x6.fixed_point_found.sum(), "d"), "X6_N_CASES": num(len(x6), "d")}
    out["X4_HI_PCT_LO"], out["X4_HI_PCT_HI"] = r(hi, "buses_beating_all_vertices_pct", ".0f")
    out["X5_LIAO_PEG_LO"], out["X5_LIAO_PEG_HI"] = r(s5.iloc[2:], "liao_width_over_reference_median", ".0f")
    out["X5_DIRECT_LO"], out["X5_DIRECT_HI"] = r(s5, "raak_width_over_reference_median", ".2f")
    found = x6[x6.fixed_point_found]
    small = found[found.case == "case33bw"]
    out["X6_ALLOW_33_LO"], out["X6_ALLOW_33_HI"] = r(small, "allowance_over_r_median_pct", ".2f")
    out["X6_W_33_LO"], out["X6_W_33_HI"] = r(small, "raak_enclosure_width_over_reference_median", ".3f")
    m118 = found[found.case == "case118"]
    out["X6_ALLOW_118_LO"], out["X6_ALLOW_118_HI"] = r(m118, "allowance_over_r_median_pct", ".0f")
    out["X6_W_118_LO"], out["X6_W_118_HI"] = r(m118, "raak_enclosure_width_over_reference_median", ".2f")
    out["X6_VIOL"] = num(found.raak_enclosure_buses_violated.sum() + found.af2_enclosure_buses_violated.sum() + found.states_outside_rho.sum()
                         + found.states_outside_sigma.sum(), "d")
    peg5 = x6[(x6.case == "pglib_opf_case2869_pegase") & np.isclose(x6.delta, 0.05)]
    span = lambda col, spec: (lambda a, b: a if a == b else f"{a} a {b}")(*r(peg5, col, spec))
    out["X6_PEG_G"] = span("spectral_radius_G", ".1f")
    out["X6_PEG_GP"] = span("spectral_radius_G_polar_first_order", ".1f")
    out["X6_PEG_TP"] = span("spectral_radius_true_polar_max_over_64_vertices", ".2f")
    out["X6_PEG_FACTOR"] = num((peg5.spectral_radius_G_polar_first_order / peg5.spectral_radius_true_polar_max_over_64_vertices).min(), ".0f")
    for name, tag in (("af2", "AF2"), ("raak", "RAAK"), ("spectral", "SPEC"), ("sdp", "SDP")):
        f = x7[x7.bound == name]
        out[f"X7_{tag}_LO"], out[f"X7_{tag}_HI"] = r(f, "width_over_vertex_range_median", ".2f" if name != "sdp" else ".3f")
        out[f"X7_{tag}_MS_LO"], out[f"X7_{tag}_MS_HI"] = r(f, "ms_per_form", ".1f" if f.ms_per_form.max() < 10 else ".0f")
    out["X7_SDP_MAX"] = num(x7[x7.bound == "sdp"].width_over_vertex_range_max.max(), ".2f")
    hyb = x7[x7.bound.str.startswith("hybrid")].groupby(["system", "balancing", "K"]).width_over_vertex_range_median.min()
    best = x7[x7.bound.isin(["raak", "spectral"])].groupby(["system", "balancing", "K"]).width_over_vertex_range_median.min()
    out["X7_HYB_WINS"] = num((hyb < best - 0.005).sum(), "d")
    out["X7_NCASES"] = num(len(best), "d")
    return out


X8_RADIAL = ["33", "69", "85", "141"]


def _x8():
    d = pd.read_csv(RES / "x8_feeders.csv", dtype={"network": str})
    d["pct"] = 100 * d.rho_U
    d["scen"] = d.scenario.map({"F": "directo", "R": "inverso"})
    return d


def _yes(col):
    return col.map(lambda v: "–" if pd.isna(v) else ("sí" if bool(v) else "no"))


def t_x8_validation():
    d = pd.read_csv(RES / "x8_wclf_validation.csv", dtype={"eps_ave_published_pct": str, "eps_max_published_pct": str, "network": str})
    d["pct"] = 100 * d.rho_U
    for c in ("eps_ave_published_pct", "eps_max_published_pct"):
        d[c] = d[c].str.replace(".", ",", regex=False)
    d["conv"] = _yes(d.converged)
    return table(d, [("table", "Tabla del artículo", None), ("network", "Red", None), ("pct", "ρU (%)", ".0f"),
                     ("eps_ave_published_pct", "ε medio publicado (%)", None), ("eps_ave_pct", "ε medio obtenido (%)", ".3g"),
                     ("eps_max_published_pct", "ε máximo publicado (%)", None), ("eps_max_pct", "ε máximo obtenido (%)", ".3g"),
                     ("conv", "Converge", None), ("iterations", "Iteraciones", "d")])


def t_x8_radial():
    d = _x8()
    d = d[d.network.isin(X8_RADIAL)]
    return table(d, [("network", "Red", None), ("scen", "Flujo", None), ("pct", "ρU (%)", ".0f"), ("attained_width_max_pu", "Rango alcanzado, máx. (p.u.)", ".4f"),
                     ("enclosure_width_over_attained_median", "Encierro / rango alcanzado, mediana", ".3f"),
                     ("enclosure_upper_excess_ave_pct", "Exceso de la cota superior: encierro (%)", ".4f"),
                     ("wclf_upper_excess_ave_pct", "ídem, WCLF (%)", ".4f"), ("wclf_certified_gap_ave_pct", "Error certificado por WCLF (%)", ".4f"),
                     ("enclosure_lower_slack_ave_pct", "Holgura de la cota inferior: encierro (%)", ".4f"),
                     ("discs_per_bus_width_over_attained_median", "Discos / rango alcanzado, mediana", ".2f"),
                     ("enclosure_ms", "Encierro (ms)", ".0f"), ("wclf_ms", "WCLF (ms)", ".0f")])


def t_x8_general():
    d = _x8()
    d = d[~d.network.isin(X8_RADIAL)].copy()
    d["fp"] = _yes(d.enclosure_fixed_point)
    d["dk"] = _yes(d.discs_ok)
    return table(d, [("network", "Red", None), ("n_bus", "Barras", "d"), ("loops", "Mallas", "d"), ("symbols", "Símbolos", "d"), ("scen", "Flujo", None),
                     ("pct", "ρU (%)", ".0f"), ("attained_width_max_pu", "Rango alcanzado, máx. (p.u.)", ".4f"), ("fp", "Punto fijo", None),
                     ("enclosure_width_over_attained_median", "Encierro / rango alcanzado, mediana", ".3f"),
                     ("enclosure_width_over_attained_max", "máx.", ".3f"), ("enclosure_upper_excess_ave_pct", "Exceso superior medio (%)", ".4f"),
                     ("enclosure_lower_slack_ave_pct", "Holgura inferior media (%)", ".4f"), ("enclosure_violations", "Barras con un valor alcanzado fuera", "d"),
                     ("second_order_only_lower_violations", "ídem, solo 2.º orden", "d"), ("dk", "Discos aplicables", None),
                     ("discs_per_bus_width_over_attained_median", "Discos / rango alcanzado, mediana", ".2f"), ("enclosure_ms", "Encierro (ms)", ".0f")])


def t_x8_boxes():
    d = pd.read_csv(RES / "x8_load_boxes.csv")
    d["pct"] = 100 * d.delta
    d["old"], d["new"] = _yes(d.rectangular_fixed_point), _yes(d.injection_fixed_point)
    for c in ("rectangular_width_over_vertex_range_median", "injection_width_over_vertex_range_median", "injection_width_over_vertex_range_max", "injection_violations"):
        if c not in d:
            d[c] = np.nan
    return table(d, [("network", "Red", None), ("K", "K", "d"), ("pct", "δ (%)", ".0f"), ("vertex_range_width_max_pu", "Rango exacto en vértices, máx. (p.u.)", ".4f"),
                     ("old", "Punto fijo, cota de la sección 22", None), ("rectangular_width_over_vertex_range_median", "Encierro / rango, mediana", ".4f"),
                     ("new", "Punto fijo, cota por inyección de corriente", None), ("injection_width_over_vertex_range_median", "Encierro / rango, mediana", ".4f"),
                     ("injection_width_over_vertex_range_max", "máx.", ".4f"), ("injection_violations", "Barras con algún vértice fuera", "d")])


def t_x8_remainder():
    d = pd.read_csv(RES / "x8_remainder_bounds.csv", dtype={"network": str})
    d = d[d.status == "ok"].copy()
    d["ratio"] = d.rectangular_third_order_bound_max_pu / d.injection_third_order_bound_max_pu
    order = list(dict.fromkeys(d.network))
    g = d.groupby("network", sort=False)
    out = pd.DataFrame({"network": order, "cases": [len(g.get_group(n)) for n in order],
                        "y": [_sci(1.0 / g.get_group(n).largest_branch_admittance_pu.iloc[0]) for n in order],
                        "old": [int(g.get_group(n).rectangular_fixed_point.sum()) for n in order],
                        "new": [int(g.get_group(n).injection_fixed_point.sum()) for n in order],
                        "ratio": [g.get_group(n).ratio.median() for n in order]})
    return table(out, [("network", "Red", None), ("cases", "Casos", "d"), ("y", "Rama más corta, \\|z\\| (p.u.)", None),
                       ("old", "Cierra la cota de la sección 22", "d"), ("new", "Cierra la cota por inyección de corriente", "d"),
                       ("ratio", "Cota del tercer orden: sección 22 / nueva, mediana donde ambas cierran", ".1f")])


def t_x8_monotone():
    d = pd.read_csv(RES / "x8_monotonicity.csv", dtype={"network": str})
    d = d[d.status == "ok"]
    order = list(dict.fromkeys(d.network))
    g = d.groupby("network", sort=False)
    out = pd.DataFrame({"network": order, "cases": [len(g.get_group(n)) for n in order], "buses": [int(g.get_group(n).pq_buses.iloc[0]) for n in order],
                        "bad": [int((g.get_group(n).buses_with_lower_minimum > 0).sum()) for n in order],
                        "most": [int(g.get_group(n).buses_with_lower_minimum.max()) for n in order],
                        "gap": [_sci(g.get_group(n).largest_gap_pu.max()) if g.get_group(n).buses_with_lower_minimum.max() else "–" for n in order],
                        "rel": [g.get_group(n).largest_gap_over_attained_width_pct.max() if g.get_group(n).buses_with_lower_minimum.max() else np.nan for n in order]})
    return table(out, [("network", "Red", None), ("cases", "Casos", "d"), ("buses", "Barras PQ", "d"),
                       ("bad", "Casos con alguna barra cuyo mínimo no está en (0, −C)", "d"), ("most", "Barras en esa situación, máx.", "d"),
                       ("gap", "Mayor diferencia (p.u.)", None), ("rel", "ídem, como % del rango alcanzado de la barra", ".1f")])


def t_x8_scale():
    d = load("x8_raak_scale")
    d["pct"] = 100 * d.delta
    d["entry"] = d.largest_entry_median.map(_sci)
    d["design"] = d.design_delta.map({True: "sí", False: "no"})
    for c in ("raw_over_vertex_range_median", "normalised_over_vertex_range_median"):
        if c not in d:
            d[c] = np.nan
    return table(d, [("system", "Sistema", None), ("K", "K", "d"), ("pct", "δ (%)", ".0f"), ("design", "δ de diseño", None), ("forms", "Formas", "d"),
                     ("entry", "Mayor coeficiente, mediana", None), ("raw_over_af2_median", "RAA-K sin normalizar / AF2", ".3f"),
                     ("normalised_over_af2_median", "RAA-K normalizada / AF2", ".3f"),
                     ("raw_over_vertex_range_median", "Sin normalizar / rango en vértices", ".3f"),
                     ("normalised_over_vertex_range_median", "Normalizada / rango en vértices", ".3f"),
                     ("forms_narrower_by_over_0p1pct", "Formas que mejoran más de 0,1 %", "d"), ("largest_narrowing_pct", "Mayor mejora (%)", ".0f")])


def x8_values():
    d = _x8()
    v = pd.read_csv(RES / "x8_wclf_validation.csv", dtype={"eps_ave_published_pct": str, "eps_max_published_pct": str, "network": str})
    b = pd.read_csv(RES / "x8_load_boxes.csv")
    ok = d[d.status == "ok"]
    rad, gen = ok[ok.network.isin(X8_RADIAL)], ok[~ok.network.isin(X8_RADIAL)]
    r = lambda fr, col, spec: (num(fr[col].min(), spec), num(fr[col].max(), spec))
    def printed(pub, got):
        """Does ``got`` round to the published string, at the number of significant digits it was printed with?"""
        digits = len(pub.lower().split("e")[0].replace(".", "").replace("-", "").lstrip("0"))
        return float(f"{got:.{max(digits, 1)}g}") == float(pub)

    same = lambda pub, got: all(printed(p, g) for p, g in zip(pub, got))
    t1 = v[(v.table == "I") & (v.rho_U > 0)]
    big = v[v.table != "I"]
    hit = [printed(p, a) and printed(q, m) for a, p, m, q in zip(big.eps_ave_pct, big.eps_ave_published_pct, big.eps_max_pct, big.eps_max_published_pct)]
    out = {"X8_N": num(len(d), "d"), "X8_N_OK": num(len(ok), "d"), "X8_N_FP": num(ok.enclosure_fixed_point.sum(), "d"),
           "X8_VIOL": num(ok.enclosure_violations.sum() + rad.wclf_violations.sum() + ok.discs_violations.fillna(0).sum()
                          + ok.discs_per_bus_violations.fillna(0).sum(), "d"),
           "X8_AC_FAIL": num(ok.ac_failures_in_samples.sum(), "d"), "X8_SAMPLES": num(meta("x8_feeders")["samples_per_case"], "d"),
           "X8_T1_OK": "sí" if same(t1.eps_ave_published_pct, t1.eps_ave_pct) and same(t1.eps_max_published_pct, t1.eps_max_pct) else "no",
           "X8_T67_HIT": num(sum(hit), "d"), "X8_T67_N": num(len(big), "d"),
           "X8_T67_NETS": "; ".join(f"{n} barras, Tabla {t}" for n, t, h in zip(big.network, big.table, hit) if h) or "ninguna",
           "X8_NONMONO": num(ok.buses_with_minimum_not_at_zero_p_minus_C.sum(), "d"),
           "X8_VMIN_GAP": _sci(rad.wclf_vmin_minus_attained_max_pu.max()),
           "X8_2ND_LOW": num((ok.second_order_only_lower_violations > 0).sum(), "d"),
           "X8_2ND_LOW_MAX": _sci(ok.second_order_only_lower_violation_max_pu.max()),
           "X8_2ND_UP": num((ok.second_order_only_upper_violations > 0).sum(), "d"),
           "X8_DISC_FAIL": num((~ok.discs_ok).sum(), "d"),
           "X8_RAD_N": num(len(rad), "d"), "X8_GEN_N": num(len(gen), "d"), "X8_GEN_FP": num(gen.enclosure_fixed_point.sum(), "d"),
           "X8_RAD_FP": num(rad.enclosure_fixed_point.sum(), "d"),
           "X8_WCLF_CONV": num(rad.wclf_converged.sum(), "d")}
    out["X8_W_RAD_LO"], out["X8_W_RAD_HI"] = r(rad, "enclosure_width_over_attained_median", ".3f")
    out["X8_W_GEN_LO"], out["X8_W_GEN_HI"] = r(gen, "enclosure_width_over_attained_median", ".3f")
    out["X8_WMAX_GEN"] = num(gen.enclosure_width_over_attained_max.max(), ".3f")
    out["X8_WIDTH_MAX"] = num(rad.attained_width_max_pu.max(), ".2f")
    nofp = gen[~gen.enclosure_fixed_point.astype(bool)]
    out["X8_GEN_NOFP"] = "; ".join(f"{n} con ρU = {100 * u:.0f} % (flujo {s})" for n, u, s in zip(nofp.network, nofp.rho_U, nofp.scen)) or "ninguno"
    out["X8_W_GEN_SMALL_HI"] = num(gen[gen.rho_U <= 0.2 + 1e-12].enclosure_width_over_attained_median.max(), ".3f")
    cg = gen[(gen.network == "CIGRE MV") & np.isclose(gen.rho_U, 0.5)]
    out["X8_CIGRE_WIDTH"] = num(cg.attained_width_max_pu.max(), ".2f")
    out["X8_CIGRE_W"] = num(cg.enclosure_width_over_attained_median.max(), ".2f")
    out["X8_DISC_GEN_LO"], out["X8_DISC_GEN_HI"] = r(gen[gen.discs_ok], "discs_per_bus_width_over_attained_median", ".2f")
    out["X8_DISCU_GEN_LO"], out["X8_DISCU_GEN_HI"] = r(gen[gen.discs_ok], "discs_width_over_attained_median", ".1f")
    out["X8_DISC_LO"], out["X8_DISC_HI"] = r(ok[ok.discs_ok], "discs_per_bus_width_over_attained_median", ".2f")
    out["X8_DISCU_LO"], out["X8_DISCU_HI"] = r(ok[ok.discs_ok], "discs_width_over_attained_median", ".1f")
    small = rad[rad.rho_U <= 0.1 + 1e-12]
    large = rad[rad.rho_U >= 0.5 - 1e-12]
    win = lambda fr: num((fr.enclosure_upper_excess_ave_pct < fr.wclf_upper_excess_ave_pct).sum(), "d")
    for tag, fr in (("ALL", rad), ("F", rad[rad.scenario == "F"]), ("R", rad[rad.scenario == "R"]), ("SMALL", small), ("LARGE_F", large[large.scenario == "F"]),
                    ("SMALL_R", small[small.scenario == "R"])):
        out[f"X8_WIN_{tag}"], out[f"X8_CNT_{tag}"] = win(fr), num(len(fr), "d")
        out[f"X8_LOSE_{tag}"] = num((fr.enclosure_upper_excess_ave_pct >= fr.wclf_upper_excess_ave_pct).sum(), "d")
    ratio = rad.wclf_upper_excess_ave_pct / rad.enclosure_upper_excess_ave_pct
    out["X8_BEST_FOR_US"] = num(ratio.max(), ".0f")
    out["X8_BEST_FOR_WCLF"] = num((1.0 / ratio).max(), ".0f")
    worst = rad.loc[(1.0 / ratio).idxmax()]
    out["X8_WORST_NET"], out["X8_WORST_RHO"] = str(worst.network), num(100 * worst.rho_U, ".0f")
    out["X8_WORST_EXC"], out["X8_WORST_WCLF"] = num(worst.enclosure_upper_excess_ave_pct, ".2f"), num(worst.wclf_upper_excess_ave_pct, ".3f")
    out["X8_GAP_SHARE_LO"], out["X8_GAP_SHARE_HI"] = r(rad.assign(s=100 * rad.attained_minus_wclf_lower_bound_ave_pct / rad.wclf_certified_gap_ave_pct), "s", ".0f")
    out["X8_LOW_RAD_MAX"] = num(rad.enclosure_lower_slack_max_pct.max(), ".2f")
    out["X8_LOW_SMALL_MAX"] = num(small.enclosure_lower_slack_max_pct.max(), ".4f")
    out["X8_MS_ENC_LO"], out["X8_MS_ENC_HI"] = r(ok, "enclosure_ms", ".0f")
    out["X8_MS_WCLF_LO"], out["X8_MS_WCLF_HI"] = r(rad, "wclf_ms", ".0f")
    ob = ok[ok.network == "Oberrhein"]
    out["X8_OB_BUS"], out["X8_OB_K"] = num(ob.n_bus.iloc[0], "d"), num(ob.symbols.iloc[0], "d")
    out["X8_OB_MS"] = num(ob.enclosure_ms.max() / 1e3, ".1f")
    out["X8_BOX_OLD"], out["X8_BOX_NEW"], out["X8_BOX_N"] = num(b.rectangular_fixed_point.sum(), "d"), num(b.injection_fixed_point.sum(), "d"), num(len(b), "d")
    out["X8_BOX_W_LO"], out["X8_BOX_W_HI"] = r(b[b.injection_fixed_point], "injection_width_over_vertex_range_median", ".3f")
    out["X8_BOX_VIOL"] = num(b.injection_violations.fillna(0).sum(), "d")
    out["X8_BOX_DMAX"] = num(100 * b[b.injection_fixed_point].delta.max(), ".0f")
    rb = pd.read_csv(RES / "x8_remainder_bounds.csv", dtype={"network": str})
    rb = rb[rb.status == "ok"]
    out["X8_RB_N"], out["X8_RB_OLD"], out["X8_RB_NEW"] = num(len(rb), "d"), num(rb.rectangular_fixed_point.sum(), "d"), num(rb.injection_fixed_point.sum(), "d")
    lost = rb[~rb.rectangular_fixed_point.astype(bool)]
    out["X8_RB_RADIUS_LO"] = num(lost.rectangular_spectral_radius_G.min(), ".2f")
    both = rb[rb.rectangular_fixed_point.astype(bool) & rb.injection_fixed_point.astype(bool)]
    ratio = both.rectangular_third_order_bound_max_pu / both.injection_third_order_bound_max_pu
    out["X8_RB_RATIO_LO"], out["X8_RB_RATIO_HI"] = num(ratio.min(), ".1f"), num(ratio.max(), ".0f")
    never = [n for n, fr in rb.groupby("network", sort=False) if fr.rectangular_fixed_point.sum() == 0]
    out["X8_RB_NEVER"] = ", ".join(never) or "ninguna"
    sc = pd.read_csv(RES / "x8_raak_scale.csv")
    des = sc[sc.design_delta]
    out["X8_SC_NRM_LO"], out["X8_SC_NRM_HI"] = r(sc, "normalised_over_af2_median", ".2f")
    out["X8_SC_RAW_LO"], out["X8_SC_RAW_HI"] = r(sc, "raw_over_af2_median", ".2f")
    out["X8_SC_DES_RAW_LO"], out["X8_SC_DES_RAW_HI"] = r(des, "raw_over_af2_median", ".2f")
    out["X8_SC_DES_NRM_LO"], out["X8_SC_DES_NRM_HI"] = r(des, "normalised_over_af2_median", ".2f")
    k8 = des[des.K == 8]
    out["X8_SC_VR_RAW_LO"], out["X8_SC_VR_RAW_HI"] = r(k8, "raw_over_vertex_range_median", ".2f")
    out["X8_SC_VR_NRM_LO"], out["X8_SC_VR_NRM_HI"] = r(k8, "normalised_over_vertex_range_median", ".2f")
    out["X8_SC_NARROW_MAX"] = num(sc.largest_narrowing_pct.max(), ".0f")
    out["X8_SC_WIDER"], out["X8_SC_FORMS"] = num(sc.normalised_wider_than_raw.sum(), "d"), f"{int(sc.forms.sum()):,}".replace(",", ".")
    out["X8_SC_SHORT"] = num(sc.normalised_short_of_vertex_range.fillna(0).sum(), "d")
    mo = pd.read_csv(RES / "x8_monotonicity.csv", dtype={"network": str})
    mo = mo[mo.status == "ok"]
    bad = mo[mo.buses_with_lower_minimum > 0]
    out["X8_MONO_N"], out["X8_MONO_CASES"] = num(len(mo), "d"), num(len(bad), "d")
    out["X8_MONO_NETS"] = " y ".join(dict.fromkeys(bad.network)) or "ninguna"
    out["X8_MONO_OK_NETS"] = ", ".join(n for n in dict.fromkeys(mo.network) if n not in set(bad.network))
    out["X8_MONO_GAP"] = _sci(bad.largest_gap_pu.max()) if len(bad) else "0"
    out["X8_MONO_GAP_PCT"] = num(bad.largest_gap_over_attained_width_pct.max(), ".1f") if len(bad) else "0"
    out["X8_MONO_BUSES"] = num(bad.buses_with_lower_minimum.max(), "d") if len(bad) else "0"
    out["X8_MS_ENC_HI_S"] = num(ok.enclosure_ms.max() / 1e3, ".1f")
    return out


TABLES = {"reproduction": t_reproduction, "b1_main": t_b1_main, "b1_bounds": t_b1_bounds, "b1_branch": t_b1_branch, "b2": t_b2,
          "b3": t_b3, "b3_all": t_b3_all, "c1_ieee33": lambda: t_c1("IEEE-33"), "c1_pegase": lambda: t_c1("PEGASE-2869"),
          "c1_other_ieee33": lambda: t_c1_other("IEEE-33"), "c1_other_pegase": lambda: t_c1_other("PEGASE-2869"),
          "c1_first": t_c1_first, "c2": t_c2, "c3_marginals": t_c3_marginals, "c3_ieee33": lambda: t_c3("IEEE-33"),
          "c3_pegase": lambda: t_c3("PEGASE-2869"), "c3_combo": t_c3_combo, "c4_kappa": t_c4_kappa, "c4": t_c4, "c4_effect": t_c4_effect,
          "c5": t_c5, "c6": t_c6, "c6_opt": t_c6_opt, "c7": t_c7, "c7_fill": t_c7_fill, "c7_fit": t_c7_fit, "c8": t_c8, "c8_time": t_c8_time,
          "x1_form": t_x1_form, "x1_vertex": t_x1_vertex, "x1_bounds": t_x1_bounds, "x1_width": t_x1_width,
          "x1_prob": t_x1_prob, "x1_radii": t_x1_radii, "x1_dep": t_x1_dep,
          "x2_setting": t_x2_setting, "x2_affine": t_x2_affine, "x2_orders": t_x2_orders, "x2_gate1": t_x2_gate1,
          "x2_width": t_x2_width, "x2_operators": t_x2_operators, "x2_gain": t_x2_gain, "x2_sweep": t_x2_sweep, "x3": t_x3,
          "x4_independent": t_x4_independent, "x4_interior": t_x4_interior, "x5_table1": t_x5_table1, "x5_variants": t_x5_variants,
          "x5_summary": t_x5_summary, "x6": t_x6, "x7": t_x7,
          "x8_validation": t_x8_validation, "x8_radial": t_x8_radial, "x8_general": t_x8_general, "x8_boxes": t_x8_boxes,
          "x8_scale": t_x8_scale, "x8_remainder": t_x8_remainder, "x8_monotone": t_x8_monotone}


def computed_values():
    """Numbers quoted inside the commentary, read from the CSVs."""
    t = pd.read_csv(RES / "c8_timing.csv").set_index("system")
    c7 = pd.read_csv(RES / "c7_timing.csv").set_index("case")
    fit = pd.read_csv(RES / "c7_fit.csv").set_index(["quantity", "systems"])
    peg, i33 = t.loc["pegase2869"], t.loc["ieee33"]
    big = c7.loc["case2869pegase"]

    def k(q, sub="all systems"):
        return num(fit.loc[(q, sub), "exponent"], ".2f")

    def ci(q, sub="all systems"):
        r = fit.loc[(q, sub)]
        return f"[{num(r.ci95_low, '.2f')}; {num(r.ci95_high, '.2f')}]"

    thr = lambda r: (r.jacobian_mean_ms + r.lu_mean_ms + r.solves_default_threads_ms + r.branch_mean_ms) / 1e3
    thousands = lambda x: f"{float(x):,.0f}".replace(",", ".")    # 3780 -> 3.780 in the commentary
    return {
        "T_PEG": num(peg.cdhim_total_s, ".2f"), "T_33": num(i33.cdhim_total_s, ".4f"), "T_PEG_THR": num(thr(big), ".2f"),
        "T_MC_PEG": num(peg.ac_mc_batched_s, ".0f"), "PP_PEG": num(peg.pandapower_ms_per_pf, ".0f"),
        "PP_33": num(i33.pandapower_ms_per_pf, ".0f"), "PP_PEG_TOT": thousands(peg.pandapower_loop_extrapolated_s),
        "SP_PEG": num(peg.speedup_vs_batched_mc, ".0f"), "SP_33": num(i33.speedup_vs_batched_mc, ".0f"),
        "SP_PEG_PP": thousands(peg.speedup_vs_pandapower_loop_extrapolated),
        "SP_33_PP": thousands(i33.speedup_vs_pandapower_loop_extrapolated),
        "T_BASE_PEG": num(peg.base_pf_pandapower_s, ".2f"), "BR_SHARE": num(100 * peg.cdhim_branch_stage_s / peg.cdhim_total_s, ".0f"),
        "K_LU": k("lu"), "K_LU_CI": ci("lu"), "K_NNZ": k("nnz(L+U), COLAMD"), "K_NNZ_CI": ci("nnz(L+U), COLAMD"),
        "K_SOLVE": k("solves"), "K_SOLVE_CI": ci("solves"), "K_TOT": k("total"), "K_TOT_CI": ci("total"),
        "K_TOT5": k("total", "n_bus >= 118"), "K_TOT5_CI": ci("total", "n_bus >= 118"),
        "SH_300": num(c7.loc["case300", "solves_share_of_total_pct"], ".0f"), "SH_2869": num(big.solves_share_of_total_pct, ".0f"),
        "SH_9241": num(c7.loc["case9241pegase", "solves_share_of_total_pct"], ".0f"),
        "T7_2869": num(big.total_mean_ms / 1e3, ".2f"), "T7_9241": num(c7.loc["case9241pegase", "total_mean_ms"] / 1e3, ".1f"),
        "T7_2869_THR": num(big.solves_default_threads_ms / 1e3, ".2f"),
        "T7_9241_THR": num(c7.loc["case9241pegase", "solves_default_threads_ms"] / 1e3, ".1f"),
    }


def x2_values():
    """Numbers quoted in the commentary of section 18, read from results/x2_summary.csv."""
    x = pd.read_csv(RES / "x2_summary.csv")
    d = x[x.design_delta]
    peg, i33 = d[d.system == "pegase2869"], d[d.system == "ieee33"]
    part, slack = peg[peg.balancing == "participation"], peg[peg.balancing == "slack"]
    p8, s8 = part[part.K == part.K.min()].iloc[0], slack[slack.K == slack.K.min()].iloc[0]
    thousands = lambda v: f"{float(v):,.0f}".replace(",", ".")
    out = {}

    def rng(tag, frame, col, spec):
        out[f"X2_{tag}_LO"], out[f"X2_{tag}_HI"] = num(frame[col].min(), spec), num(frame[col].max(), spec)

    out.update({
        "X2_S3_WIDER_P8": thousands(p8.sigma3_buses_wider_than_affine), "X2_NBUS_P8": thousands(p8.n_buses),
        "X2_S3_SYS_P8": num(p8.sigma3_mc_system_pct, ".0f"),
        "X2_T1_MAX": num(peg.first_order_seconds.max(), ".2f"), "X2_T2_MAX": num(peg.second_order_seconds.max(), ".1f"),
        "X2_TB_MAX": num(peg.bounds_seconds.max(), ".1f"), "X2_DV_MAX": num(peg.max_abs_dv_at_vertices_pu.max(), ".3f"),
        "X2_AFF_MCS_P8": num(p8.affine_mc_samples_outside_pct, ".0f"), "X2_AFF_MCS_S8": num(s8.affine_mc_samples_outside_pct, ".0f"),
        "X2_AFF_SHORT_MAX": num(peg.affine_vertex_shortfall_max_pu.max(), ".3f"),
        "X2_MONO_SHORT_I": f"{i33.monomial_vertex_shortfall_max_pu.max():.1e}".replace(".", ",").replace("e-0", "·10⁻").replace("⁻6", "⁻⁶").replace("⁻5", "⁻⁵"),
        "X2_SPEC_OUT_P8": thousands(p8.spectral_buses_with_vertex_outside), "X2_SPEC_N4_P8": thousands(p8["spectral_buses_shortfall_over_0.0001pu"]),
        "X2_SPEC_N3_P8": thousands(p8["spectral_buses_shortfall_over_0.001pu"]), "X2_SPEC_SHORT_P8": num(p8.spectral_vertex_shortfall_max_pu, ".4f"),
        "X2_SPEC_OUT_S8": thousands(s8.spectral_buses_with_vertex_outside), "X2_SPEC_N4_S8": thousands(s8["spectral_buses_shortfall_over_0.0001pu"]),
        "X2_SPEC_N3_S8": thousands(s8["spectral_buses_shortfall_over_0.001pu"]), "X2_SPEC_SHORT_S8": num(s8.spectral_vertex_shortfall_max_pu, ".4f"),
        "X2_SPEC_SHORT_MAX": num(peg.spectral_vertex_shortfall_max_pu.max(), ".4f"),
        "X2_SPEC_MCS_MAX": num(d.spectral_mc_samples_outside_pct.max(), ".2f"), "X2_SPEC_MCB_MAX": thousands(d.spectral_mc_buses_outside.max()),
        "X2_W_SPEC_MAX": num(d.spectral_width_over_reference_max.max(), ".2f"), "X2_W_MONO_MAX": num(d.monomial_width_over_reference_max.max(), ".2f"),
    })
    rng("AFF_MC", peg, "affine_width_over_mc9973_median", ".1f")
    rng("AFF_N3", peg, "affine_buses_shortfall_over_0.001pu", "d")
    rng("D2_MED", peg, "second_order_depth_over_r_median_pct", ".0f")
    rng("D3_MED", peg, "third_order_over_r_median_pct", ".2f")
    rng("32_MED", peg, "third_over_second_order_median_pct", ".1f")
    rng("32_P90", peg, "third_over_second_order_p90_pct", ".0f")
    rng("D3_PU", peg, "third_order_max_pu", ".4f")
    rng("NSD", peg, "buses_h_negative_semidefinite_pct", ".0f")
    rng("W_SPEC", d, "spectral_width_over_reference_median", ".2f")
    rng("W_NAT_P", part, "natural_width_over_reference_median", ".0f")
    rng("W_NAT_S", slack, "natural_width_over_reference_median", ".0f")
    rng("G", d, "af2_over_spectral_form_width_median", ".2f")
    for tag, frame in (("I", i33), ("P", peg)):
        rng(f"GS_{tag}", frame, "af2_over_spectral_form_width_median", ".2f")
        for k in ("k2", "k3", "k6", "kn"):
            rng(f"GR{k.upper()}_{tag}", frame, f"af2_over_raak_{k}_form_width_median", ".2f")
    ex = d[d.n_buses_exact_form_range > 0]
    for kk in (8, 16):                                     # per number of zones, over both systems
        frame = ex[ex.K == kk]
        rng(f"EX_RKN_K{kk}", frame, "raak_kn_form_width_over_exact_median", ".2f")
        rng(f"EX_SPEC_K{kk}", frame, "spectral_form_width_over_exact_median", ".2f")
        rng(f"EX_AF2_K{kk}", frame, "monomial_form_width_over_exact_median", ".2f")
    i4 = ex[(ex.system == "ieee33") & (ex.K == 4)].iloc[0]
    out.update({"X2_EX_RKN_I4": num(i4.raak_kn_form_width_over_exact_median, ".2f"),
                "X2_EX_SPEC_I4": num(i4.spectral_form_width_over_exact_median, ".2f"),
                "X2_EX_AF2_I4": num(i4.monomial_form_width_over_exact_median, ".2f")})
    rng("RISE_RKN", ex, "raak_kn_positive_side_over_depth_median_pct", ".0f")
    rng("GRK2", d, "af2_over_raak_k2_form_width_median", ".2f")
    for kk in (8, 64):
        rng(f"GRKN_P{kk}", peg[peg.K == kk], "af2_over_raak_kn_form_width_median", ".2f")
    out["X2_KN_WORSE"] = thousands((d.af2_over_raak_kn_form_width_median < d.af2_over_raak_k6_form_width_median - 0.005).sum())
    out["X2_NCASES"] = thousands(len(d))
    # lower side of every bound over the exact depth (K <= 16, negative semidefinite forms)
    bus = pd.read_csv(RES / "x2_bus.csv")
    bus = bus[(bus.K <= 16) & bus.h_negative_semidefinite & (bus.vertex_q_min < 0)]
    design = {(r.system, r.balancing, r.K): r.delta for r in d.itertuples()}
    bus = bus[[design.get((r.system, r.balancing, r.K)) == r.delta for r in bus.itertuples()]]
    for name, col in (("AF2", "monomial_lo"), ("RKN", "raak_kn_lo"), ("SPEC", "spectral_lo")):
        ratio = (bus[col] / bus.vertex_q_min).groupby([bus.system, bus.balancing, bus.K]).median()
        lo_, hi_ = num(ratio.min(), ".2f"), num(ratio.max(), ".2f")
        out[f"X2_LOW_{name}"] = lo_ if lo_ == hi_ else f"{lo_} a {hi_}"
    rng("RISE", d[d.n_buses_exact_form_range > 0], "monomial_positive_side_over_depth_median_pct", ".0f")
    assert x.vertices_without_ac_solution.sum() == 0 and x.mc_without_ac_solution.sum() == 0     # the text says so
    bus = pd.read_csv(RES / "x2_bus.csv")
    bus = bus[(bus.system == "pegase2869") & (bus.balancing == "participation") & (bus.K == p8.K) & (bus.delta == p8.delta)]
    wb = bus.loc[bus.vertex_third_order_max.idxmax()]
    out.update({"X2_WB_BUS": str(int(wb.bus)), "X2_WB_R": num(wb.r_affine, ".3f"), "X2_WB_D2": num(-wb.vertex_q_min, ".3f"),
                "X2_WB_D3": num(wb.vertex_third_order_max, ".4f"),
                "X2_AFF_EXC_I_MED": num(i33.affine_low_excess_over_r_median_pct.max(), ".1f"),
                "X2_AFF_EXC_I_MAX": num(i33.affine_low_excess_over_r_max_pct.max(), ".1f"),
                "X2_ANG_RATIO_LO": num((part.angle_radius_median_rad.values / slack.angle_radius_median_rad.values).min(), ".2f"),
                "X2_ANG_RATIO_HI": num((part.angle_radius_median_rad.values / slack.angle_radius_median_rad.values).max(), ".2f")})
    sw = x[(x.system == "pegase2869") & (x.balancing == "participation") & (x.K == part.K.min())].set_index("delta")
    sci = lambda v: f"{v:.1e}".replace(".", ",").replace("e-05", "·10⁻⁵").replace("e-04", "·10⁻⁴").replace("e-03", "·10⁻³")
    out.update({"X2_SW_P1": sci(sw.loc[0.01, "spectral_vertex_shortfall_max_pu"]), "X2_SW_P2": sci(sw.loc[0.02, "spectral_vertex_shortfall_max_pu"]),
                "X2_SW_P5": sci(sw.loc[0.05, "spectral_vertex_shortfall_max_pu"]),
                "X2_SW_N4_P2": thousands(sw.loc[0.02, "spectral_buses_shortfall_over_0.0001pu"])})
    return out


def x1_values():
    """Numbers quoted in the commentary of section 17, read from results/x1_summary.csv."""
    x = pd.read_csv(RES / "x1_summary.csv").set_index("system")
    forms = pd.read_csv(RES / "x1_quadratic_bounds.csv")
    thousands = lambda v: f"{float(v):,.0f}".replace(",", ".")
    out = {}
    for tag, key in (("33", "ieee33"), ("PEG", "pegase2869")):
        r = x.loc[key]
        for name, col, spec in (
                ("RED", "err_reduction_by_quadratic_median", ".0f"), ("NAT1", "forms_natural_over_inner_median", ".1f"),
                ("BIN1", "forms_binomial_terms_over_inner_median", ".1f"), ("MONO", "forms_monomial_over_inner_median", ".2f"),
                ("K2", "forms_greedy_k2_over_inner_median", ".2f"), ("KN1", "forms_knomial_1_over_inner_median", ".2f"),
                ("KN8", "forms_knomial_8_over_inner_median", ".2f"), ("SPEC", "forms_spectral_over_inner_median", ".2f"),
                ("VX_N", "vertex_lo_buses_outside_affine_box", "d"), ("N_ACT", "n_buses_above_floor", "d"),
                ("VX_MED", "vertex_lo_excess_over_r_median_pct", ".1f"), ("VX_P90", "vertex_lo_excess_over_r_p90_pct", ".1f"),
                ("VX_MAX", "vertex_lo_excess_over_r_max_pct", ".0f"), ("VX_MAX1", "vertex_lo_excess_over_r_max_pct", ".1f"),
                ("VX_10", "vertex_lo_buses_excess_over_10pct", "d"), ("VX_100", "vertex_lo_buses_excess_over_100pct", "d"),
                ("VX_MISS", "vertex_lo_second_order_misses_median_pct_of_excess", ".1f"),
                ("VXH_N", "vertex_hi_buses_outside_affine_box", "d"),
                ("ERR_AFF_MAX", "lin_err_over_r_affine_max_pct", ".1f"), ("ERR_3S_MAX", "lin_err_over_r_3sigma_max_pct", ".0f"),
                ("ZON", "zonotope_gain_median_pct", ".0f"), ("ZON_ELL", "zonotope_over_ellipse_area_median", ".1f"),
                ("ZON_LT", "branches_zonotope_smaller_than_ellipse", "d"), ("NBR", "n_branches", "d"),
                ("DROP", "drop_reduction_median_pct", ".0f"), ("HOEF", "hoeffding_over_affine_median", ".2f"),
                ("S3_MAX", "viol_sigma3_plain_max_pct", ".2f"), ("S3S_MAX", "viol_sigma3_shift_max_pct", ".2f"),
                ("S3_SYS", "viol_sigma3_plain_system_pct", ".1f"), ("S3S_SYS", "viol_sigma3_shift_system_pct", ".1f"),
                ("SH_PV", "natural_share_pv_squares_median_pct", ".0f"), ("SH_CURV", "natural_share_curvature_median_pct", ".0f"),
                ("HOEF_SYS", "viol_hoeffding_plain_system_pct", ".2f"), ("HOEF_MAX", "viol_hoeffding_plain_max_pct", ".3f"),
                ("ANG", "angle_radius_median_rad", ".3f"), ("MAG", "magnitude_radius_median_pu", ".4f"),
                ("CLTM6", "cltm_buses_below_6_effective_terms_above_floor", "d"), ("NFORMS", "forms_n", "d"),
                ("RISE_MONO", "forms_monomial_spurious_rise_over_depth_median_pct", ".0f"),
                ("RISE_SPEC", "forms_spectral_spurious_rise_over_depth_median_pct", ".1f"),
                ("NEG", "forms_negative_share_of_spectrum_median_pct", ".1f"), ("EIG90", "forms_eig90_median", ".0f")):
            out[f"X_{name}_{tag}"] = thousands(r[col]) if spec == "d" else num(r[col], spec)
        f = forms[forms.system == key]                       # same bus, same attained range: ratio of the two widths
        for name, col in (("GAIN1", "knomial_1_over_inner"), ("GAIN8", "knomial_8_over_inner"), ("GAINK2", "greedy_k2_over_inner")):
            out[f"X_{name}_{tag}"] = num((f.monomial_over_inner / f[col]).median(), ".2f")
    peg = x.loc["pegase2869"]
    out["X_NAT_PEG_K"] = thousands(peg.forms_natural_over_inner_median)
    out["X_NATW_PEG_K"] = thousands(peg.natural_width_over_r_affine_median_pct / 100.0)
    out["X_DEPTH_WORST_MED"] = num(peg.forms_worst_depth_over_r_affine_median_pct / 100.0, ".1f")
    out["X_DEPTH_WORST_MAX"] = num(peg.forms_worst_depth_over_r_affine_max_pct / 100.0, ".1f")
    out["X_DEPTH_RND_MED"] = num(peg.forms_random_depth_over_r_affine_median_pct / 100.0, ".2f")
    out["X_DEPTH_RND_MAX"] = num(peg.forms_random_depth_over_r_affine_max_pct / 100.0, ".2f")
    out["X_DEPTH_33"] = num(x.loc["ieee33", "forms_all_depth_over_r_affine_median_pct"], ".1f")
    out["X_HOEF_3S"] = num(x.loc["ieee33", "hoeffding_over_3sigma"], ".2f")
    return out


def main():
    import results_text
    import results_text_c
    import results_text_end as end

    c_text = (results_text_c.C_TEXT.replace("<<C1_LECTURA>>", end.C1_LECTURA.strip())
              .replace("<<C7_LECTURA>>", end.C7_LECTURA.strip()).replace("<<C8_LECTURA>>", end.C8_LECTURA.strip()))
    import results_text_x

    import results_text_x2
    import results_text_x4
    import sys

    import report_x12
    import results_text_x8
    import results_text_x12
    import report_b45
    import results_text_b45

    here = sys.modules[__name__]
    tables = {**TABLES, **report_x12.tables(here), **report_b45.tables(here)}
    text = results_text.TEXT.replace("<<C_SECTIONS>>", c_text.strip() + "\n\n" + end.FINAL_TEXT.strip() + "\n\n"
                                     + results_text_x.X_TEXT.strip() + "\n\n" + results_text_x2.X2_TEXT.strip() + "\n\n"
                                     + results_text_x2.X3_TEXT.strip() + "\n\n" + results_text_x4.X4_TEXT.strip() + "\n\n"
                                     + results_text_x8.X8_TEXT.strip() + "\n\n" + results_text_x12.X12_TEXT.strip() + "\n\n"
                                     + results_text_b45.B45_TEXT.strip() + "\n")
    for key, fn in tables.items():
        tag = "{{" + key + "}}"
        if tag in text:
            text = text.replace(tag, fn())
    for key, value in {**computed_values(), **x1_values(), **x2_values(), **x3_values(), **x4_values(), **x8_values(), **report_x12.values(here), **report_b45.values(here)}.items():
        text = text.replace(f"<<{key}>>", value)
    v = meta("b2_conservatism")["versions"]
    text = text.replace("{{versions}}", f"Python {v['python']}, pandapower {v['pandapower']}, NumPy {v['numpy']}, SciPy {v['scipy']} ({v['platform']})")
    text = text.replace("{{seed}}", str(CFG["experiments"]["seed"]))
    assert "{{" not in text and "<<" not in text, [w for w in text.split() if "<<" in w or "{{" in w][:10]
    (ROOT / "RESULTS.md").write_text(text, encoding="utf-8")
    print(ROOT / "RESULTS.md", len(text), "chars")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:                       # print single tables for inspection
        for k in sys.argv[1:]:
            print(f"\n## {k}\n")
            print(TABLES[k]())
    else:
        main()
