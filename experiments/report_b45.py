"""Tables and numbers of sections 27 to 29 of RESULTS.md (B4, B5 and the list of claims). ``h`` is the
build_report module: its ``num``, ``table``, ``RES`` and ``SYS`` keep the formatting of the other sections."""
import numpy as np
import pandas as pd


def _csv(h, name):
    return pd.read_csv(h.RES / f"{name}.csv")


def _sys(s):
    return {"ieee33": "IEEE-33", "pegase2869": "PEGASE-2869"}.get(s, s)


# ----------------------------------------------------------------------------- B4
B4_BUS = [("cdhim_lmax", "none", "Intervalo CD-HIM (λ_max), sin corregir"),
          ("cdhim_lmax", "moment", "  + media y varianza de q"),
          ("cdhim_lmax", "range", "  + rango de q (Caja: q en el vértice)"),
          ("cdhim_coord", "none", "Intervalo CD-HIM (proyección exacta), sin corregir"),
          ("cdhim_coord", "range", "  + rango de q"),
          ("cltm_3sigma", "none", "Caja 3σ, sin corregir"),
          ("cltm_3sigma", "range", "  + rango de q"),
          ("affine_box", "none", "Caja afín, sin corregir"),
          ("affine_box", "box2", "  + cota rigurosa de q (término a término)"),
          ("affine_box", "box2v", "  + q en el vértice extremo")]
B4_ELL = [("none", "Elipse, sin corregir"), ("shift", "  + media de q"), ("moment", "  + media y covarianza de q"),
          ("minkowski", "  + rango de q en los ejes principales")]


def t_b4_design(h):
    bus, br = _csv(h, "b4_design_bus"), _csv(h, "b4_design_branch")
    rows = []
    for s in ("ieee33", "pegase2869"):
        for bound, kind, label in B4_BUS:
            g = bus[(bus.system == s) & (bus.bound == bound) & (bus.correction == kind)].set_index("reference")
            rows.append({"sys": _sys(s), "what": label, "unit": "barra", "ac": g.loc["ac", "max_pct"], "n": g.loc["ac", "n_buses_over_0.27pct"],
                         "lin": g.loc["lin", "max_pct"], "so": g.loc["so", "max_pct"], "w": g.loc["ac", "width_ratio_median"],
                         "wmax": g.loc["ac", "width_ratio_max"]})
        for kind, label in B4_ELL:
            g = br[(br.system == s) & (br.correction == kind)].set_index("reference")
            rows.append({"sys": _sys(s), "what": label, "unit": "rama", "ac": g.loc["ac", "max_pct"], "n": g.loc["ac", "n_branches_over_0.27pct"],
                         "lin": g.loc["lin", "max_pct"], "so": g.loc["so", "max_pct"], "w": g.loc["ac", "area_ratio_median"],
                         "wmax": g.loc["ac", "area_ratio_max"]})
    return h.table(pd.DataFrame(rows), [("sys", "Sistema", None), ("what", "Cota", None), ("unit", "Por", None),
                                        ("ac", "Violación máx. contra AC (%)", ".3f"), ("n", "Sobre 0,27 %", "d"),
                                        ("lin", "Máx. contra el modelo lineal (%)", ".3f"), ("so", "Máx. contra el sustituto de 2.º orden (%)", ".3f"),
                                        ("w", "Ancho (área) / sin corregir, mediana", ".3f"), ("wmax", "máximo", ".2f")])


def _sweep_rows(h, d, items, prefix_fmt, width_fmt=None):
    rows = []
    for s in ("ieee33", "pegase2869"):
        g = d[(d.system == s) & (d.status == "ok")]
        for key, label in items:
            nov = g[prefix_fmt.format(ref="ac", key=key) + "_n_over"]
            mx = g[prefix_fmt.format(ref="ac", key=key) + "_max_pct"]
            wcol = width_fmt.format(key=key) if width_fmt else None
            w = g[wcol] if wcol and wcol in g else None
            rows.append({"sys": _sys(s), "what": label, "cells": len(g), "bad": int((nov > 0).sum()), "nmax": int(nov.max()),
                         "worst": float(mx.max()), "w": float(w.median()) if w is not None else np.nan,
                         "wmax": float(w.max()) if w is not None else np.nan})
    return rows


SWEEP_COLS = [("sys", "Sistema", None), ("what", "Cota", None), ("cells", "Celdas", "d"),
              ("bad", "Celdas con alguna barra o rama sobre 0,27 %", "d"), ("nmax", "Máx. barras o ramas sobre 0,27 % en una celda", "d"),
              ("worst", "Peor violación (%)", ".3f"), ("w", "Ancho / sin corregir, mediana de las celdas", ".3f"), ("wmax", "máximo", ".3f")]


def t_b4_sweep(h):
    d = _csv(h, "b4_sweep")
    items = [("cdhim_lmax_none", "Intervalo CD-HIM, sin corregir"), ("cdhim_lmax_moment", "  + media y varianza de q"),
             ("cdhim_lmax_range", "  + rango de q"), ("cltm_3sigma_none", "Caja 3σ, sin corregir"),
             ("cltm_3sigma_moment", "  + media y varianza de q"), ("affine_box_none", "Caja afín, sin corregir"),
             ("affine_box_box2v", "  + q en el vértice extremo"), ("ellipse_none", "Elipse, sin corregir"),
             ("ellipse_moment", "  + media y covarianza de q"), ("ellipse_minkowski", "  + rango de q en los ejes principales")]
    return h.table(pd.DataFrame(_sweep_rows(h, d, items, "{ref}_{key}", "width_{key}_median")), SWEEP_COLS)


def t_b4_timing(h):
    d = _csv(h, "b4_timing").groupby("system").mean(numeric_only=True).reset_index()
    d["sys"] = d.system.map(_sys)
    return h.table(d, [("sys", "Sistema", None), ("linear_method_s", "Método lineal (s)", ".3f"),
                       ("second_order_build_s", "Modelo de 2.º orden (s)", ".3f"), ("mean_trace_s", "Media exacta de q (s)", ".3f"),
                       ("vertex_box_buses_s", "q en los vértices de las barras Caja (s)", ".3f"),
                       ("surrogate_1000_s", "1.000 muestras de q (s)", ".3f"), ("surrogate_4000_s", "4.000 muestras de q (s)", ".3f")])


# ----------------------------------------------------------------------------- B5
B5_BUS = [("det_lin", "Caja afín (peor caso, 1.er orden)"), ("det_v", "  + q en el vértice extremo"),
          ("det_s", "  + extremo del 2.º orden hallado por búsqueda"), ("bern", "Bernstein + rango de q"),
          ("sig3", "3σ + rango de q (supone normalidad)"), ("cdhim", "Intervalo CD-HIM + rango de q (B4)")]
B5_DROP = [("det_lin", "Forma afín de la caída, peor caso"), ("det_v", "  + q en el vértice extremo"),
           ("det_s", "  + extremo hallado por búsqueda"), ("bern", "Forma afín de la caída, Bernstein + rango de q"),
           ("sig3", "Forma afín de la caída, 3σ + rango de q"), ("sub_det_v", "Resta de las cajas afines corregidas de las dos barras"),
           ("sub_bern", "Resta de los intervalos Bernstein de las dos barras"), ("sub_cdhim", "Resta de los intervalos CD-HIM corregidos")]


def t_b5_design(h):
    mc, w = _csv(h, "b5_mc"), _csv(h, "b5_width")
    rows = []
    for s in ("ieee33", "pegase2869"):
        for tg, items in (("bus", B5_BUS), ("drop", B5_DROP)):
            for k, label in items:
                g = mc[(mc.system == s) & (mc.target == tg) & (mc.bound == k)].set_index("reference")
                ww = w[(w.system == s) & (w.target == tg) & (w.bound == k)].set_index("relative_to")
                rows.append({"sys": _sys(s), "tg": "barra" if tg == "bus" else "caída", "what": label,
                             "ac": g.loc["ac", "max_pct"], "n": g.loc["ac", "n_over_alpha"], "sysl": g.loc["ac", "system_level_pct"],
                             "lin": g.loc["lin", "max_pct"], "waff": ww.loc["det_lin", "ratio_median"],
                             "wsig": ww.loc["sig3", "ratio_median"], "wpu": ww.loc["det_lin", "median_width_pu"]})
    return h.table(pd.DataFrame(rows), [("sys", "Sistema", None), ("tg", "Objetivo", None), ("what", "Cota", None),
                                        ("ac", "Violación máx. contra AC (%)", ".3f"), ("n", "Sobre 0,27 %", "d"),
                                        ("sysl", "Muestras con alguna fuera (%)", ".2f"), ("lin", "Máx. contra el modelo lineal (%)", ".3f"),
                                        ("waff", "Ancho / caja afín, mediana", ".2f"), ("wsig", "Ancho / 3σ, mediana", ".2f"),
                                        ("wpu", "Ancho mediano (p.u.)", ".5f")])


def t_b5_vertex(h):
    v = _csv(h, "b5_vertex")
    names = {"lin-": "1.er orden, mínimo", "lin+": "1.er orden, máximo", "so-": "búsqueda de 2.º orden, mínimo", "so+": "búsqueda de 2.º orden, máximo"}
    v["sys"] = v.system.map(_sys)
    v["tg"] = v.target.map({"bus": "barra", "drop": "caída"})
    v["vx"] = v.vertex.map(names)
    return h.table(v, [("sys", "Sistema", None), ("tg", "Objetivo", None), ("vx", "Vértice", None), ("n_targets", "Objetivos", "d"),
                       ("n_not_converged", "AC sin converger", "d"), ("det_lin_n_outside", "Fuera de la caja afín", "d"),
                       ("det_v_n_outside", "Fuera de la caja + q en el vértice", "d"), ("det_s_n_outside", "Fuera de la caja + búsqueda", "d"),
                       ("det_s_max_excess_pu", "Exceso máx. sobre esta última (p.u.)", ".6f"),
                       ("det_s_max_excess_over_width", "Exceso máx. / ancho", ".3f")])


def t_b5_sweep(h):
    d = _csv(h, "b5_sweep")
    items = [("bus_det_lin", "Barra: caja afín"), ("bus_det_v", "Barra: caja afín + q en el vértice"),
             ("bus_bern", "Barra: Bernstein + rango de q"), ("bus_sig3", "Barra: 3σ + rango de q"),
             ("bus_cdhim", "Barra: CD-HIM + rango de q"), ("drop_det_v", "Caída: forma afín, peor caso + q en el vértice"),
             ("drop_bern", "Caída: forma afín, Bernstein + rango de q"), ("drop_sig3", "Caída: forma afín, 3σ + rango de q"),
             ("drop_sub_bern", "Caída: resta de intervalos Bernstein"), ("drop_sub_cdhim", "Caída: resta de intervalos CD-HIM")]
    rows = []
    for s in ("ieee33", "pegase2869"):
        g = d[(d.system == s) & (d.status == "ok")]
        for key, label in items:
            nov, mx = g[f"ac_{key}_n_over"], g[f"ac_{key}_max_pct"]
            tg, k = key.split("_", 1)
            w = g[f"width_{tg}_{k}_over_affine_median"] if f"width_{tg}_{k}_over_affine_median" in g else None
            rows.append({"sys": _sys(s), "what": label, "cells": len(g), "bad": int((nov > 0).sum()), "nmax": int(nov.max()),
                         "worst": float(mx.max()), "w": float(w.median()) if w is not None else np.nan,
                         "wmax": float(w.max()) if w is not None else np.nan})
    cols = [c if c[0] not in ("w", "wmax") else (c[0], "Ancho / caja afín, mediana de las celdas" if c[0] == "w" else "máximo", c[2]) for c in SWEEP_COLS]
    return h.table(pd.DataFrame(rows), cols)


def t_b5_surrogate(h):
    d = _csv(h, "b5_surrogate")
    g = d.groupby(["system", "target", "n_surrogate"]).agg(mx=("max_pct", "max"), n=("n_over_alpha", "max"),
                                                          w=("width_over_uncorrected_median", "median"), t=("t_range_s", "mean")).reset_index()
    g["sys"] = g.system.map(_sys)
    g["tg"] = g.target.map({"bus": "barra", "drop": "caída"})
    g["m"] = g.n_surrogate.map(lambda x: "0 (sin corregir)" if x == 0 else f"{x}")
    return h.table(g, [("sys", "Sistema", None), ("tg", "Objetivo", None), ("m", "Muestras de q", None),
                       ("mx", "Peor violación, 3 semillas (%)", ".3f"), ("n", "Sobre 0,27 %, peor semilla", "d"),
                       ("w", "Ancho / sin corregir, mediana", ".4f"), ("t", "Tiempo del muestreo (s)", ".2f")])


def t_b5_radii(h):
    d = _csv(h, "b5_radii")
    names = {"gauss_3sigma": "3σ (normal)", "bernstein": "Bernstein (soporte y varianza)", "hoeffding": "Hoeffding (soporte)",
             "unimodal_gauss_ineq": "Gauss, unimodal (varianza)", "chebyshev": "Chebyshev (varianza)", "affine": "Caja afín (peor caso)"}
    d = d[d.radius.isin(names)].copy()
    d["sys"] = d.system.map(_sys)
    d["tg"] = d.target.map({"bus": "barra", "drop": "caída"})
    d["r"] = d.radius.map(names)
    piv = d.pivot_table(index=["sys", "tg", "r"], columns="correction",
                        values=["max_pct", "n_over_alpha"], aggfunc="first").reset_index()
    piv.columns = ["_".join([c for c in col if c]) for col in piv.columns]
    base = d[d.correction == "second_order"][["sys", "tg", "r", "halfwidth_over_sigma_median", "width_over_affine_median", "capped_by_affine_pct"]]
    out = base.merge(piv, on=["sys", "tg", "r"])
    order = {v: i for i, v in enumerate(names.values())}
    out = out.sort_values(["sys", "tg", "r"], key=lambda s: s.map(order) if s.name == "r" else s)
    return h.table(out, [("sys", "Sistema", None), ("tg", "Objetivo", None), ("r", "Radio (qué supone)", None),
                         ("halfwidth_over_sigma_median", "Semiancho / σ, mediana", ".2f"), ("width_over_affine_median", "Ancho / caja afín, mediana", ".2f"),
                         ("capped_by_affine_pct", "Topado por la caja afín (% de objetivos)", ".0f"),
                         ("max_pct_none", "Violación máx. AC sin corregir (%)", ".3f"), ("n_over_alpha_none", "sobre 0,27 %", "d"),
                         ("max_pct_second_order", "con el rango de q (%)", ".3f"), ("n_over_alpha_second_order", "sobre 0,27 %", "d")])


_B6_NAMES = {"independent": "independientes", "global_rho0.155": "factor global, ρ = 0,155", "global_rho0.5": "factor global, ρ = 0,5",
             "global_rho0.845": "factor global, ρ = 0,845", "zones_rho0.155": "8 zonas, ρ = 0,155", "zones_rho0.5": "8 zonas, ρ = 0,5",
             "zones_rho0.845": "8 zonas, ρ = 0,845", "copula_r0.3": "cópula de C2, r = 0,3", "copula_r0.6": "cópula de C2, r = 0,6",
             "copula_r0.9": "cópula de C2, r = 0,9"}


def t_b6(h, target):
    d = _csv(h, "b6_correlated")
    d = d[d.target == target]
    rows = []
    for (s, model, bal), g in d.groupby(["system", "model", "balancing"], sort=False):
        g = g.set_index("bound")
        own = "bern_model" if "bern_model" in g.index else "bern_matched_factor"
        rows.append({"sys": _sys(s), "model": _B6_NAMES[model], "bal": "slack" if bal == "slack" else "generadores",
                     "ind": g.loc["bern_independent", "max_pct"], "own": g.loc[own, "max_pct"], "ownn": g.loc[own, "n_over_alpha"],
                     "w": g.loc[own, "width_over_affine_median"], "s3": g.loc["sig3_true", "max_pct"],
                     "aff": g.loc["affine", "max_pct"], "nc": g.loc["affine", "n_ac_not_converged"]})
    order = list(_B6_NAMES.values())
    df = pd.DataFrame(rows)
    df["o"] = df.model.map(order.index)
    df = df.sort_values(["sys", "o", "bal"], kind="stable")
    return h.table(df, [("sys", "Sistema", None), ("model", "Correlación de las fuentes", None), ("bal", "Balance", None),
                        ("ind", "Bernstein suponiendo independencia: violación máx. (%)", ".3f"),
                        ("own", "Bernstein con el modelo de factores (cópula: factor con el mismo ρ): violación máx. (%)", ".3f"),
                        ("ownn", "sobre 0,27 %", "d"), ("w", "su ancho / caja afín, mediana", ".2f"),
                        ("s3", "3σ con la Σ_ε verdadera: violación máx. (%)", ".3f"),
                        ("aff", "Caja afín (peor caso lineal): violación máx. (%)", ".3f"), ("nc", "AC sin converger", "d")])


_CASES = {"case33bw": "IEEE-33", "pglib_opf_case2869_pegase": "PEGASE-2869", "case118": "IEEE-118", "case1354pegase": "PEGASE-1354"}
_LEVEL = {"per_target": "por objetivo (α)", "bonferroni_own": "simultáneo, Bonferroni", "bonferroni_both": "simultáneo barras y caídas"}
_RAD = {"bernstein": "Bernstein", "gauss": "Gauss", "chebyshev": "Chebyshev", "affine": "caja afín"}


def t_b7_design(h):
    d = _csv(h, "b7_generality")
    d = d[(d.status == "ok") & ((d.load_level == "design") | ((d.load_factor == 1.0) & (d.delta == 0.05)))]
    d = d[d.level != "bonferroni_both"].copy()
    d["sys"] = d.case.map(_CASES)
    d["tg"] = d.target.map({"bus": "barra", "drop": "caída"})
    d["lv"] = d.level.map(_LEVEL)
    d["r"] = d.radius.map(_RAD)
    d["o"] = d.case.map(list(_CASES).index)
    d = d.sort_values(["o", "target", "level"], kind="stable")
    return h.table(d, [("sys", "Sistema", None), ("delta", "δ", ".2f"), ("tg", "Objetivo", None), ("lv", "Nivel", None), ("r", "Radio", None),
                       ("max_pct", "Peor objetivo (%)", ".3f"), ("any_target_pct", "Muestras con alguno fuera (%)", ".3f"),
                       ("any_bus_or_drop_pct", "con alguna barra o caída fuera (%)", ".3f"),
                       ("width_over_affine_median", "Ancho / caja afín, mediana", ".2f"),
                       ("capped_by_affine_pct", "Topado por la caja afín (%)", ".0f")])


def t_b7_grid(h):
    d = _csv(h, "b7_generality")
    d = d[(d.status == "ok") & (d.load_level != "design") & (d.level != "bonferroni_both")]
    rows = []
    for (case, tg, lv, r), g in d.groupby(["case", "target", "level", "radius"], sort=False):
        rows.append({"sys": _CASES[case], "tg": "barra" if tg == "bus" else "caída", "lv": _LEVEL[lv], "r": _RAD[r], "cells": len(g),
                     "bad": int((g.n_over_alpha > 0).sum()), "worst": g.max_pct.max(), "anyw": g.any_bus_or_drop_pct.max(),
                     "w": g.width_over_affine_median.median(), "nc": int(g.n_ac_not_converged.max())})
    return h.table(pd.DataFrame(rows), [("sys", "Sistema", None), ("tg", "Objetivo", None), ("lv", "Nivel", None), ("r", "Radio", None),
                                        ("cells", "Celdas", "d"), ("bad", "Celdas con algún objetivo sobre 0,27 %", "d"),
                                        ("worst", "Peor objetivo (%)", ".3f"), ("anyw", "Peor celda: muestras con alguna barra o caída fuera (%)", ".3f"),
                                        ("w", "Ancho / caja afín, mediana de las celdas", ".2f"), ("nc", "Máx. muestras AC sin converger en una celda", "d")])


_METH = {"bern": "Bernstein + rango de q", "sig3": "3σ + rango de q", "cumulant_cf": "cumulantes + Cornish-Fisher",
         "pem_cf": "estimación por puntos (Hong 2m+1, AC) + Cornish-Fisher", "affine": "caja afín + rango de q",
         "bern_linear_only": "Bernstein sin el rango de q"}
_LVL = {"design": "diseño", "stress_delta": "carga 1,2, δ = 30 %", "near_limit": "95 % del límite, δ = 15 %"}


def t_b8(h):
    d = _csv(h, "b8_ppf_baselines")
    d = d[d.target == "bus"].copy()
    d["sys"] = d.case.map(_CASES)
    d["cell"] = d.load_level.map(_LVL)
    d["m"] = d.method.map(_METH)
    d["o"] = d.case.map(list(_CASES).index)
    d["mo"] = d.method.map(list(_METH).index)
    d = d.sort_values(["o", "load_factor", "mo"], kind="stable")
    return h.table(d, [("sys", "Sistema", None), ("cell", "Celda", None), ("m", "Método", None), ("max_pct", "Peor barra (%)", ".3f"),
                       ("n_over_alpha", "Barras sobre 0,27 %", "d"), ("any_target_pct", "Muestras con alguna barra fuera (%)", ".2f"),
                       ("width_over_bern_median", "Ancho / Bernstein, mediana", ".2f"), ("n_ac_solves", "Flujos AC del método", "d"),
                       ("t_method_s", "Tiempo del método (s)", ".2f"), ("n_ac_not_converged", "AC sin converger (de 20.000)", "d")])


def t_b9(h):
    d = _csv(h, "b9_branch_currents").copy()
    d["sys"] = d.case.map(_CASES)
    d["cell"] = d.load_level.map(_LVL)
    d["m"] = d.method.map(_METH)
    d["o"] = d.case.map(list(_CASES).index)
    d["mo"] = d.method.map(list(_METH).index)
    d = d.sort_values(["o", "load_factor", "mo"], kind="stable")
    return h.table(d, [("sys", "Sistema", None), ("cell", "Celda", None), ("n_branch_ends", "Extremos de rama", "d"), ("m", "Método", None),
                       ("max_pct", "Peor extremo (%)", ".3f"), ("n_over_alpha", "Sobre 0,27 %", "d"),
                       ("any_end_pct", "Muestras con alguno fuera (%)", ".2f"), ("width_over_bern_median", "Ancho / Bernstein, mediana", ".2f"),
                       ("rel_halfwidth_median_pct", "Semiancho / |I₀|, mediana (%)", ".1f"),
                       ("max_err_linear_pu", "Error máx. del modelo lineal (p.u.)", ".4f"),
                       ("max_err_second_order_pu", "del estado de 2.º orden (p.u.)", ".5f")])


_MARG = ["uniform", "triangular", "truncnorm_3s", "beta25_support", "beta52_support", "beta0505_ushape", "student_t3_trunc",
         "lognormal_trunc", "beta25_var13", "beta25_paper", "student_t3", "lognormal"]
_RAD10 = {"hoeffding": "Hoeffding", "bern_var13": "Bernstein, var 1/3", "bern_true_var": "Bernstein, var verdadera",
          "sig3_true_var": "3σ, var verdadera", "affine": "caja afín"}


def t_b10(h):
    d = _csv(h, "b10_marginals")
    d = d[d.target == "bus"]
    rows = []
    for (case, m), g in d.groupby(["case", "marginal"], sort=False):
        g = g.set_index("radius")
        row = {"sys": _CASES[case], "marg": m, "inm": "sí" if bool(g.in_model.iloc[0]) else "no", "var": g.eps_variance.iloc[0],
               "o": list(_CASES).index(case), "mo": _MARG.index(m) if m in _MARG else 99}
        for r in _RAD10:
            row[r] = g.loc[r, "max_pct"]
            row[r + "_n"] = g.loc[r, "n_over_alpha"]
        rows.append(row)
    df = pd.DataFrame(rows).sort_values(["o", "mo"])
    cols = [("sys", "Sistema", None), ("marg", "Distribución (C3)", None), ("inm", "Media 0 y soporte en [−1, 1]", None), ("var", "Varianza", ".3f")]
    for r, lab in _RAD10.items():
        cols += [(r, f"{lab}: peor barra (%)", ".3f"), (r + "_n", "sobre 0,27 %", "d")]
    return h.table(df, cols)


def t_b5_n250(h):
    d4, d2 = _csv(h, "b5_sweep"), _csv(h, "b5_sweep_n250")
    rows = []
    for tag, d in (("4.000", d4), ("250", d2)):
        d = d[d.status == "ok"]
        for s in ("ieee33", "pegase2869"):
            g = d[d.system == s]
            for tg in ("bus", "drop"):
                rows.append({"n": tag, "sys": _sys(s), "tg": "barra" if tg == "bus" else "caída", "cells": len(g),
                             "bad": int((g[f"ac_{tg}_bern_n_over"] > 0).sum()), "worst": g[f"ac_{tg}_bern_max_pct"].max(),
                             "w": g[f"width_{tg}_bern_over_affine_median"].median()})
    return h.table(pd.DataFrame(rows), [("n", "Muestras de q", None), ("sys", "Sistema", None), ("tg", "Objetivo", None), ("cells", "Celdas", "d"),
                                        ("bad", "Celdas con alguno sobre 0,27 %", "d"), ("worst", "Peor violación (%)", ".3f"),
                                        ("w", "Ancho / caja afín, mediana", ".3f")])


def t_b11(h):
    d = _csv(h, "b11_certified")
    rows = []
    for (case, tg, kind, m), g in d.groupby(["case", "target", "kind", "n_surrogate"], sort=False):
        rows.append({"sys": _CASES[case], "o": list(_CASES).index(case), "tg": {"bus": "barra", "drop": "caída", "current": "corriente"}[tg],
                     "to": ["bus", "drop", "current"].index(tg), "kind": "sin certificar (B5)" if kind == "uncertified" else "certificada",
                     "ko": 0 if kind == "uncertified" else 1, "m": int(m), "al": g.alpha_linear.iloc[0], "cells": len(g),
                     "bad": int((g.n_over_alpha > 0).sum()), "worst": g.max_pct.max(), "w": g.width_over_uncertified_median.median(),
                     "wmax": g.width_over_uncertified_max.max()})
    df = pd.DataFrame(rows).sort_values(["o", "to", "ko", "m"])
    return h.table(df, [("sys", "Sistema", None), ("tg", "Objetivo", None), ("kind", "Cota", None), ("m", "Muestras de q", "d"),
                        ("al", "Nivel de la parte lineal α_L", ".5f"), ("cells", "Celdas", "d"), ("bad", "Celdas con alguno sobre 0,27 %", "d"),
                        ("worst", "Peor violación (%)", ".3f"), ("w", "Ancho / sin certificar, mediana", ".3f"), ("wmax", "máximo", ".2f")])


def t_b11_timing(h):
    d = _csv(h, "b11_timing").groupby("case").median(numeric_only=True).reset_index()
    d["sys"] = d.case.map(_CASES)
    d["o"] = d.case.map(list(_CASES).index)
    d = d.sort_values("o")
    return h.table(d, [("sys", "Sistema", None), ("n_bus", "Barras PQ", "d"), ("n_current_ends", "Extremos de rama", "d"),
                       ("linear_study_s", "Estudio lineal (s)", ".2f"), ("second_order_model_s", "Modelo de 2.º orden (s)", ".2f"),
                       ("samples_voltage_targets_s", "4.000 muestras, barras y caídas (s)", ".2f"),
                       ("samples_all_targets_s", "4.000 muestras, con corrientes (s)", ".2f"),
                       ("total_voltages_s", "Total, voltajes (s)", ".2f"), ("total_with_currents_s", "Total, con corrientes (s)", ".2f")])


def tables(h):
    return {"b4_design": lambda: t_b4_design(h), "b4_sweep": lambda: t_b4_sweep(h), "b4_timing": lambda: t_b4_timing(h),
            "b5_design": lambda: t_b5_design(h), "b5_vertex": lambda: t_b5_vertex(h), "b5_sweep": lambda: t_b5_sweep(h),
            "b5_surrogate": lambda: t_b5_surrogate(h), "b5_radii": lambda: t_b5_radii(h),
            "b6_bus": lambda: t_b6(h, "bus"), "b6_drop": lambda: t_b6(h, "drop"),
            "b7_design": lambda: t_b7_design(h), "b7_grid": lambda: t_b7_grid(h),
            "b8": lambda: t_b8(h), "b11": lambda: t_b11(h), "b11_timing": lambda: t_b11_timing(h), "b9": lambda: t_b9(h), "b10": lambda: t_b10(h), "b5_n250": lambda: t_b5_n250(h)}


# ----------------------------------------------------------------------------- numbers quoted in the text
def values(h):
    n = h.num
    bus, br = _csv(h, "b4_design_bus"), _csv(h, "b4_design_branch")
    pb = bus[(bus.system == "pegase2869") & (bus.reference == "ac")]
    pr = br[(br.system == "pegase2869") & (br.reference == "ac")].set_index("correction")
    sel = lambda b, k: pb[(pb.bound == b) & (pb.correction == k)].iloc[0]
    s4 = _csv(h, "b4_sweep")
    s4 = s4[s4.status == "ok"]
    tm = _csv(h, "b4_timing").groupby("system").mean(numeric_only=True)
    c8 = _csv(h, "c8_timing").set_index("system")
    mc, w, vx = _csv(h, "b5_mc"), _csv(h, "b5_width"), _csv(h, "b5_vertex")
    pm = mc[(mc.system == "pegase2869") & (mc.reference == "ac")]
    pw = w[(w.system == "pegase2869")]
    get = lambda df, **kw: df.loc[np.logical_and.reduce([df[k] == v for k, v in kw.items()])].iloc[0]
    vb = vx[(vx.system == "pegase2869") & (vx.target == "bus")]
    vd = vx[(vx.system == "pegase2869") & (vx.target == "drop")]
    vi = vx[(vx.system == "ieee33") & (vx.target == "bus")]
    s5 = _csv(h, "b5_sweep")
    s5 = s5[s5.status == "ok"]
    sur = _csv(h, "b5_surrogate")
    sp = sur[(sur.system == "pegase2869") & (sur.target == "bus")]
    rad = _csv(h, "b5_radii")
    rp = rad[(rad.system == "pegase2869") & (rad.target == "bus") & (rad.correction == "second_order")].set_index("radius")
    t_method = float(sp.t_study_s.iloc[0] + sp.t_run_and_second_order_s.iloc[0] + sp[sp.n_surrogate == 250].t_range_s.mean())
    t_mc = float(c8.loc["pegase2869", "ac_mc_batched_s"])
    drop_ratio = lambda k: get(pw, target="drop", bound=k, relative_to="bern").ratio_median
    b6 = _csv(h, "b6_correlated")
    p6 = b6[(b6.system == "pegase2869") & (b6.target == "bus")]
    fac = p6[p6.model.str.startswith(("global", "zones"))]
    cop = p6[p6.model.str.startswith("copula")]
    b6v = {
        "B6_IND_MIN": n(fac[fac.bound == "bern_independent"].max_pct.min(), ".1f"),
        "B6_IND_MAX": n(p6[p6.bound == "bern_independent"].max_pct.max(), ".0f"),
        "B6_FAC_MAX": n(fac[fac.bound == "bern_model"].max_pct.max(), ".3f"),
        "B6_FAC_WMIN": n(fac[fac.bound == "bern_model"].width_over_affine_median.min(), ".2f"),
        "B6_FAC_WMAX": n(fac[fac.bound == "bern_model"].width_over_affine_median.max(), ".2f"),
        "B6_IND_W": n(p6[(p6.model == "independent") & (p6.bound == "bern_model")].width_over_affine_median.max(), ".2f"),
        "B6_COP_MATCH": n(cop[cop.bound == "bern_matched_factor"].max_pct.max(), ".1f"),
        "B6_COP_AFF_SLACK": n(cop[(cop.bound == "affine") & (cop.balancing == "slack")].max_pct.max(), ".2f"),
        "B6_COP_AFF_PART": n(cop[(cop.bound == "affine") & (cop.balancing == "participation")].max_pct.max(), ".3f"),
        "B6_COP_S3": n(cop[cop.bound == "sig3_true"].max_pct.max(), ".2f"),
    }
    b7 = _csv(h, "b7_generality")
    b7 = b7[b7.status == "ok"]
    dz = b7[(b7.load_level == "design") | ((b7.load_factor == 1.0) & (b7.delta == 0.05))]
    pick = lambda case, tg, lv, r, col: dz[(dz.case == case) & (dz.target == tg) & (dz.level == lv) & (dz.radius == r)][col].iloc[0]
    P = "pglib_opf_case2869_pegase"
    grid7 = b7[b7.load_level != "design"]
    b7v = {
        "B7_BERN_PT": n(pick(P, "bus", "per_target", "bernstein", "width_over_affine_median"), ".2f"),
        "B7_BERN_BF": n(pick(P, "bus", "bonferroni_own", "bernstein", "width_over_affine_median"), ".2f"),
        "B7_GAUSS_BF": n(pick(P, "bus", "bonferroni_own", "gauss", "width_over_affine_median"), ".2f"),
        "B7_GAUSS_BF_ANY": n(pick(P, "bus", "bonferroni_own", "gauss", "any_target_pct"), ".3f"),
        "B7_GAUSS_PT_ANY": n(pick(P, "bus", "per_target", "gauss", "any_target_pct"), ".0f"),
        "B7_BERN_PT_ANY": n(pick(P, "bus", "per_target", "bernstein", "any_bus_or_drop_pct"), ".3f"),
        "B7_BERN_PT_ANY_GRID": n(grid7[(grid7.level == "per_target") & (grid7.radius == "bernstein")].any_bus_or_drop_pct.max(), ".3f"),
        "B7_W118": n(pick("case118", "bus", "per_target", "bernstein", "width_over_affine_median"), ".2f"),
        "B7_W1354": n(pick("case1354pegase", "bus", "per_target", "bernstein", "width_over_affine_median"), ".2f"),
        "B7_CAP118": n(pick("case118", "bus", "per_target", "bernstein", "capped_by_affine_pct"), ".0f"),
        "B7_GAUSS_BAD1354": n(int((grid7[(grid7.case == "case1354pegase") & (grid7.target == "bus") & (grid7.level == "per_target")
                                         & (grid7.radius == "gauss")].n_over_alpha > 0).sum()), "d"),
        "B7_NC_MAX": n(100 * grid7.n_ac_not_converged.max() / grid7.n_samples.max(), ".0f"),
    }
    b8 = _csv(h, "b8_ppf_baselines")
    b8b = b8[b8.target == "bus"]
    b9 = _csv(h, "b9_branch_currents")
    b10 = _csv(h, "b10_marginals")
    b10b = b10[b10.target == "bus"]
    n250 = _csv(h, "b5_sweep_n250")
    n250 = n250[n250.status == "ok"]
    p8 = b8b[(b8b.case == P) & (b8b.load_level == "design")].set_index("method")
    b8v = {
        "B8_BERN_MAX": n(b8b[b8b.method == "bern"].max_pct.max(), ".2f"),
        "B8_CUM_MAX": n(b8b[b8b.method == "cumulant_cf"].max_pct.max(), ".1f"),
        "B8_PEM_MAX": n(b8b[b8b.method == "pem_cf"].max_pct.max(), ".1f"),
        "B8_CUM_N": n(p8.loc["cumulant_cf", "n_over_alpha"], "d"), "B8_PEM_N": n(p8.loc["pem_cf", "n_over_alpha"], "d"),
        "B8_NB": n(int(dz[(dz.case == P) & (dz.target == "bus") & (dz.level == "per_target") & (dz.radius == "bernstein")].n_bus_targets.iloc[0]), "d"),
        "B8_W": n(b8b[b8b.method.isin(["cumulant_cf", "pem_cf"])].width_over_bern_median.median(), ".2f"),
        "B8_PEM_SOLVES": n(p8.loc["pem_cf", "n_ac_solves"], "d"), "B8_PEM_T": n(p8.loc["pem_cf", "t_method_s"], ".0f"),
        "B9_BERN_MAX": n(b9[b9.method == "bern"].max_pct.max(), ".2f"),
        "B9_LIN_MAX": n(b9[b9.method == "bern_linear_only"].max_pct.max(), ".0f"),
        "B9_CUM_MAX": n(b9[b9.method.isin(["cumulant_cf", "pem_cf"])].max_pct.max(), ".0f"),
        "B9_ERRLIN": n(b9.max_err_linear_pu.max(), ".1f"),
        "B9_AFF_W": n(b9[(b9.case == P) & (b9.method == "affine")].width_over_bern_median.median(), ".2f"),
        "B10_IN_MAX": n(b10b[b10b.in_model & b10b.radius.isin(["hoeffding", "bern_var13", "bern_true_var"])].max_pct.max(), ".3f"),
        "B10_U_MAX": n(b10b[(b10b.marginal == "beta0505_ushape") & (b10b.radius == "bern_var13")].max_pct.max(), ".3f"),
        "B10_OUT_MAX": n(b10b[~b10b.in_model & (b10b.radius == "affine")].max_pct.max(), ".1f"),
        "B10_S3_MAX": n(b10b[b10b.in_model & (b10b.radius == "sig3_true_var")].max_pct.max(), ".1f"),
        "B5_N250_MAX": n(n250.ac_bus_bern_max_pct.max(), ".3f"),
    }
    b11 = _csv(h, "b11_certified")
    c4 = b11[(b11.kind == "certified") & (b11.n_surrogate == 4000)]
    c1 = b11[(b11.kind == "certified") & (b11.n_surrogate == 1000)]
    tq = b11[(b11.case == P) & (b11.n_surrogate == 4000)].t_samples_s.max()
    b11v = {"B11_W4": n(100 * (c4.width_over_uncertified_median.median() - 1), ".1f"),
            "B11_W1": n(100 * (c1.width_over_uncertified_median.median() - 1), ".1f"),
            "B11_WMAX4": n(c4.width_over_uncertified_max.max(), ".2f"),
            "B11_WORST": n(b11[b11.kind == "certified"].max_pct.max(), ".3f"), "B11_T4": n(tq, ".0f")}
    tm11 = _csv(h, "b11_timing").groupby("case").median(numeric_only=True)
    c8t = _csv(h, "c8_timing").set_index("system")
    tv, tc = tm11.loc[P, "total_voltages_s"], tm11.loc[P, "total_with_currents_s"]
    mc, pp = c8t.loc["pegase2869", "ac_mc_batched_s"], c8t.loc["pegase2869", "pandapower_loop_extrapolated_s"]
    b11v.update({"B11_TV": n(tv, ".1f"), "B11_TC": n(tc, ".1f"), "B11_SPV": n(mc / tv, ".0f"), "B11_SPC": n(mc / tc, ".0f"),
                 "B11_SPPP": n(pp / tc, ".0f"), "B11_MC": n(mc, ".0f")})
    return {**b6v, **b7v, **b8v, **b11v,
        "B4_LMAX_0": n(sel("cdhim_lmax", "none").max_pct, ".2f"), "B4_LMAX_R": n(sel("cdhim_lmax", "range").max_pct, ".2f"),
        "B4_LMAX_RW": n(100 * (sel("cdhim_lmax", "range").width_ratio_median - 1), ".1f"),
        "B4_ELL_0": n(pr.loc["none", "max_pct"], ".2f"), "B4_ELL_0N": n(pr.loc["none", "n_branches_over_0.27pct"], "d"),
        "B4_ELL_M": n(pr.loc["moment", "max_pct"], ".2f"), "B4_ELL_MN": n(pr.loc["moment", "n_branches_over_0.27pct"], "d"),
        "B4_ELL_K": n(pr.loc["minkowski", "max_pct"], ".3f"), "B4_ELL_KA": n(100 * (pr.loc["minkowski", "area_ratio_median"] - 1), ".1f"),
        "B4_ELL_KAM": n(pr.loc["minkowski", "area_ratio_max"], ".1f"),
        "B4_SW_CELLS": n(len(s4), "d"),
        "B4_SW_MOM_BAD": n(int((s4[s4.system == "pegase2869"].ac_cdhim_lmax_moment_n_over > 0).sum()), "d"),
        "B4_T_LIN": n(tm.loc["pegase2869", "linear_method_s"], ".1f"),
        "B4_T_COR": n(tm.loc["pegase2869", ["second_order_build_s", "mean_trace_s", "vertex_box_buses_s", "surrogate_4000_s"]].sum(), ".1f"),
        "B5_VX_BUS": n(int(vb.det_s_n_outside.max()), "d"), "B5_VX_NBUS": n(int(vb.n_targets.iloc[0]), "d"),
        "B5_VX_EXC": n(vb.det_s_max_excess_pu.max(), ".4f"), "B5_VX_REL": n(100 * vb.det_s_max_excess_over_width.max(), ".0f"),
        "B5_VX_DROP": n(int(vd.det_s_n_outside.max()), "d"), "B5_VX_NDROP": n(int(vd.n_targets.iloc[0]), "d"),
        "B5_VX_33": h._sci(vi.det_s_max_excess_pu.max()) if hasattr(h, "_sci") else n(vi.det_s_max_excess_pu.max(), ".6f"),
        "B5_BERN_BUS": n(get(pm, target="bus", bound="bern").max_pct, ".3f"),
        "B5_BERN_DROP": n(get(pm, target="drop", bound="bern").max_pct, ".3f"),
        "B5_BERN_AFF": n(get(pw, target="bus", bound="bern", relative_to="det_lin").ratio_median, ".2f"),
        "B5_BERN_SIG": n(get(pw, target="bus", bound="bern", relative_to="sig3").ratio_median, ".2f"),
        "B5_CDHIM_BERN": n(get(pw, target="bus", bound="cdhim", relative_to="bern").ratio_median, ".2f"),
        "B5_SUB_CDHIM": n(drop_ratio("sub_cdhim"), ".1f"), "B5_SUB_BERN": n(drop_ratio("sub_bern"), ".1f"),
        "B5_SUB_DET": n(drop_ratio("sub_det_v"), ".1f"),
        "B5_SUB_CDHIM_33": n(get(w[w.system == "ieee33"], target="drop", bound="sub_cdhim", relative_to="bern").ratio_median, ".1f"),
        "B5_SW_CELLS": n(len(s5), "d"),
        "B5_SW_BERN_BUS": n(s5.ac_bus_bern_max_pct.max(), ".3f"), "B5_SW_BERN_DROP": n(s5.ac_drop_bern_max_pct.max(), ".3f"),
        "B5_SW_SIG_BAD": n(int((s5.ac_bus_sig3_n_over > 0).sum()), "d"),
        "B5_T_METHOD": n(t_method, ".1f"), "B5_T_MC": n(t_mc, ".0f"), "B5_SPEEDUP": n(t_mc / t_method, ".0f"),
        "B5_R_BERN": n(rp.loc["bernstein", "halfwidth_over_sigma_median"], ".1f"), "B5_R_HOEF": n(rp.loc["hoeffding", "halfwidth_over_sigma_median"], ".1f"),
        "B5_R_CHEB": n(rp.loc["chebyshev", "halfwidth_over_sigma_median"], ".1f"),
        "B5_R_CHEB_CAP": n(rp.loc["chebyshev", "capped_by_affine_pct"], ".0f"), "B5_R_UNI_CAP": n(rp.loc["unimodal_gauss_ineq", "capped_by_affine_pct"], ".0f"),
        "B5_R_BH": n(100 * (1 - rp.loc["bernstein", "halfwidth_over_sigma_median"] / rp.loc["hoeffding", "halfwidth_over_sigma_median"]), ".0f"),
    }
