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


def tables(h):
    return {"b4_design": lambda: t_b4_design(h), "b4_sweep": lambda: t_b4_sweep(h), "b4_timing": lambda: t_b4_timing(h),
            "b5_design": lambda: t_b5_design(h), "b5_vertex": lambda: t_b5_vertex(h), "b5_sweep": lambda: t_b5_sweep(h),
            "b5_surrogate": lambda: t_b5_surrogate(h), "b5_radii": lambda: t_b5_radii(h)}


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
    return {
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
