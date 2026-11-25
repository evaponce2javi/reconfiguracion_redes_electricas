"""Tables and numbers of section 26 of RESULTS.md (X9, X10, X12, X13, X14). ``h`` is the build_report module:
its ``num``, ``table``, ``_sci`` and ``RES`` are used so that the formatting is the one of the other sections."""
import numpy as np
import pandas as pd


def _span(h, s, spec):
    """'a' if all equal after formatting, else 'a a b'."""
    lo, hi = h.num(s.min(), spec), h.num(s.max(), spec)
    return lo if lo == hi else f"{lo} a {hi}"


THREE_PHASE = "European LV, 3 phases"
_NET = {THREE_PHASE: "European LV, trifásica"}


def _x12(h, three_phase=None):
    """Cases of x12_verified.csv with an AC solution at the centre. ``three_phase``: None = all, False = the
    networks of X8 only, True = the three-phase feeder only."""
    d = pd.read_csv(h.RES / "x12_verified.csv", dtype={"network": str})
    d = d[d.status == "ok"].copy()
    if three_phase is not None:
        d = d[(d.network == THREE_PHASE) == three_phase]
    return d


_REASON = {"the tube reaches zero voltage": "el tubo alcanza voltaje cero", "the model reaches zero voltage": "el modelo alcanza voltaje cero",
           "the iteration diverged": "la iteración diverge", "the iteration did not settle": "la iteración no se asienta",
           "the inflated fixed point is not mapped into itself": "el punto fijo inflado no se mapea en sí mismo",
           "the bound of the modulus breaks down": "la cota del módulo deja de valer"}


def _share(h, s):
    lo, hi = s.min(), s.max()
    if lo == hi:
        return "todas las barras" if lo == 1.0 else ("ninguna" if lo == 0.0 else f"{h.num(100 * lo, '.0f')} % de las barras")
    return f"{h.num(100 * lo, '.0f')} a {h.num(100 * hi, '.0f')} % de las barras"


def t_x12_verified(h):
    d = _x12(h)
    rows = []
    for net, g in d.groupby("network", sort=False):
        v = g[g.verified.astype(bool)]
        rows.append({"network": _NET.get(net, net), "nodes": int(g.pq_nodes.iloc[0]), "K": int(g.symbols.iloc[0]), "cases": len(g),
                     "ver": int(g.verified.astype(bool).sum()), "uni": int(g.unique.astype(bool).sum()),
                     "sep": _span(h, v.separate_bound_width_over_certified_inner_median, ".3f") if len(v) else "–",
                     "wc": _span(h, v.width_over_certified_inner_median, ".3f") if len(v) else "–",
                     "wmax": h.num(v.width_over_certified_inner_max.max(), ".3f") if len(v) else "–",
                     "contra": int(v.contradictions_with_certified_points.sum()) if len(v) else 0,
                     "pts": int(g.points_not_certified.sum()), "ms": h.num(g.verified_ms.median(), ".0f")})
    return h.table(pd.DataFrame(rows), [("network", "Red", None), ("nodes", "Barras PQ", "d"), ("K", "K", "d"), ("cases", "Casos", "d"),
                                        ("ver", "Verificado", "d"), ("uni", "Con unicidad", "d"),
                                        ("sep", "Ancho / rango exacto, a lo más: cota separada", None), ("wc", "con la cota conjunta", None),
                                        ("wmax", "peor barra", None), ("contra", "Barras en contradicción con un punto certificado", "d"),
                                        ("pts", "Puntos sin certificar", "d"), ("ms", "Verificación, mediana (ms)", None)])


def _excess_ratio(g):
    """(Krawczyk - 1) / (enclosure - 1), both as width over the certified inner range, per case."""
    return (g.krawczyk_width_over_certified_inner_median - 1.0) / (g.width_over_certified_inner_median - 1.0)


def t_x12_krawczyk(h):
    d = _x12(h)
    rows = []
    for net, g in d.groupby("network", sort=False):
        both = g[g.verified.astype(bool) & g.krawczyk_closes.astype(bool)]
        rows.append({"network": _NET.get(net, net), "cases": len(g), "k": int(g.krawczyk_closes.astype(bool).sum()), "v": int(g.verified.astype(bool).sum()),
                     "kw": _span(h, both.krawczyk_width_over_certified_inner_median, ".3f") if len(both) else "–",
                     "vw": _span(h, both.width_over_certified_inner_median, ".3f") if len(both) else "–",
                     "ratio": _span(h, both.krawczyk_width_over_verified_median, ".3f") if len(both) else "–",
                     "exc": _span(h, _excess_ratio(both), ".0f") if len(both) else "–",
                     "tight": h.num(100 * max(both.krawczyk_lower_bound_tighter_share.max(), both.krawczyk_upper_bound_tighter_share.max()), ".0f") if len(both) else "–",
                     "kms": h.num(g.krawczyk_ms.median(), ".0f"), "vms": h.num(g.verified_ms.median(), ".0f")})
    return h.table(pd.DataFrame(rows), [("network", "Red", None), ("cases", "Casos", "d"), ("k", "Krawczyk cierra", "d"), ("v", "Encierro verifica", "d"),
                                        ("kw", "Krawczyk / rango exacto, a lo más (mediana)", None), ("vw", "Encierro / rango exacto, a lo más (mediana)", None),
                                        ("ratio", "Ancho Krawczyk / ancho encierro", None), ("exc", "Exceso Krawczyk / exceso encierro", None),
                                        ("tight", "Barras en que Krawczyk es más ajustado, máx. (%)", None),
                                        ("kms", "Krawczyk, mediana (ms)", None), ("vms", "Encierro, mediana (ms)", None)])


def _exact(h):
    d = pd.read_csv(h.RES / "x12_exact_range.csv", dtype={"network": str})
    d["scenario"] = d.scenario.fillna("")
    d["kind"] = np.where(d["set"].str.startswith("load"), "Demanda en 8 zonas", "Cajas envolventes de los semidiscos")
    return d


def t_x12_exact(h):
    d = _exact(h)
    rows = []
    for (kind, net, scen), g in d.groupby(["kind", "network", "scenario"], sort=False):
        g = g.sort_values("size")
        rows.append({"kind": kind, "network": net, "scen": scen or "–", "sizes": " · ".join(h.num(100 * s, ".0f") for s in g["size"]),
                     "boxes": _span(h, g.enclosures, "d"), "gap": h.num(g.bracket_gap_max_pct.max(), ".3f"),
                     "ratio": " · ".join(h.num(r, ".4f") for r in g.enclosure_over_exact_at_most_median),
                     "worst": h.num(g.enclosure_over_exact_at_most_max.max(), ".4f"), "contra": int(g.contradictions.sum())})
    return h.table(pd.DataFrame(rows), [("kind", "Conjunto", None), ("network", "Red", None), ("scen", "Flujo", None), ("sizes", "Tamaños: δ o ρU (%)", None),
                                        ("boxes", "Encierros usados", None), ("gap", "Error del rango exacto, máx. (%)", None),
                                        ("ratio", "Encierro / rango exacto, mediana, por tamaño", None), ("worst", "peor barra", None),
                                        ("contra", "Contradicciones", "d")])


def _sub(h):
    d = pd.read_csv(h.RES / "x12_subdivision.csv", dtype={"network": str})
    d["kind"] = np.where(d["set"].str.startswith("load"), "Demanda en 8 zonas", "Cajas envolventes de los semidiscos")
    return d


def t_x12_subdivision(h):
    d = _sub(h)
    rows = []
    for (kind, net), g in d.groupby(["kind", "network"], sort=False):
        row = {"kind": kind, "network": net, "cases": len(g)}
        for tag in ("verified", "krawczyk", "discs"):
            root = g[f"{tag}_root_over_exact_median"].dropna()
            row[f"{tag}_root"] = _span(h, root, ".3f") if len(root) else "–"
            row[f"{tag}_n"] = _span(h, g[f"{tag}_enclosures"], "d")
            row[f"{tag}_fail"] = int((~g[f"{tag}_reached_target"].astype(bool)).sum())
        rows.append(row)
    return h.table(pd.DataFrame(rows), [("kind", "Conjunto", None), ("network", "Red", None), ("cases", "Casos", "d"),
                                        ("verified_root", "Sin subdividir / rango exacto: encierro", None), ("krawczyk_root", "Krawczyk", None), ("discs_root", "disco", None),
                                        ("verified_n", "Encierros hasta 1 %: encierro", None), ("krawczyk_n", "Krawczyk", None), ("krawczyk_fail", "no llega", "d"),
                                        ("discs_n", "disco", None), ("discs_fail", "no llega", "d")])


def x12_values(h):
    num, out = h.num, {}
    d = _x12(h, False)
    v = d[d.verified.astype(bool)]
    out["X12_N"], out["X12_VER"], out["X12_UNI"], out["X12_X8"] = num(len(d), "d"), num(d.verified.astype(bool).sum(), "d"), num(d.unique.astype(bool).sum(), "d"), num(d.x8_fixed_point.astype(bool).sum(), "d")
    bad = d[~d.verified.astype(bool)]
    out["X12_FAIL"] = ("No verifica: " + "; ".join(f"{r.network}, {r.scenario}, ρU = {100 * r.rho_U:.0f} % ({_REASON.get(r.reason, r.reason)})" for r in bad.itertuples()) + ".") if len(bad) else "No falla en ninguno."
    sepw = v.separate_bound_width_over_certified_inner_median / v.width_over_certified_inner_median * v.width_over_x8_width_median
    out["X12_SEPW_LO"], out["X12_SEPW_HI"] = num(sepw.min(), ".3f"), num(sepw.max(), ".3f")
    out["X12_TUBE_LO"], out["X12_TUBE_HI"] = num(v.tube_radius_over_x8_remainder_bound_median.min(), ".2f"), num(v.tube_radius_over_x8_remainder_bound_median.max(), ".2f")
    out["X12_JOINT_UP"], out["X12_JOINT_LO"] = _share(h, v.joint_bound_share_upper), _share(h, v.joint_bound_share_lower)
    out["X12_SEP_LO"], out["X12_SEP_HI"] = num(v.separate_bound_width_over_certified_inner_median.min(), ".4f"), num(v.separate_bound_width_over_certified_inner_median.max(), ".3f")
    out["X12_ROUND"] = num(100 * v.rounding_share_of_tube.max(), ".1f")
    out["X12_TUBE_W"] = num(100 * (v.tube_radius_max_pu / v.attained_width_max_pu).max(), ".1f")
    out["X12_MS_LO"], out["X12_MS_HI"] = num(d.verified_ms.min(), ".0f"), num(d.verified_ms.max() / 1e3, ".1f")
    out["X12_PTS_BAD"] = num(d.points_not_certified.sum(), "d")
    out["X12_CERT_GAP"] = h._sci(d.certified_minus_attained_max_pu.max())
    out["X12_WC_LO"], out["X12_WC_HI"] = num(v.width_over_certified_inner_median.min(), ".4f"), num(v.width_over_certified_inner_median.max(), ".3f")
    out["X12_WC_MAX"] = num(v.width_over_certified_inner_max.max(), ".3f")
    out["X12_WC_20"] = num(v[v.rho_U <= 0.2].width_over_certified_inner_median.max(), ".4f")
    out["X12_CONTRA"] = num(v.contradictions_with_certified_points.sum(), "d")
    # Krawczyk against the verified enclosure
    both = d[d.verified.astype(bool) & d.krawczyk_closes.astype(bool)]
    tight = max(both.krawczyk_lower_bound_tighter_share.max(), both.krawczyk_upper_bound_tighter_share.max())
    out["X12_KR_TIGHTER"] = "0 barras en los %d casos en que ambos cierran" % len(both) if tight == 0 else f"hasta {num(100 * tight, '.0f')} % de las barras de un caso"
    out["X12_KR_OVER_LO"], out["X12_KR_OVER_HI"] = num(both.krawczyk_width_over_verified_median.min(), ".3f"), num(both.krawczyk_width_over_verified_median.max(), ".3f")
    ratio = _excess_ratio(both)
    out["X12_KR_EXC_LO"], out["X12_KR_EXC_HI"] = num(ratio.min(), ".0f"), num(ratio.max(), ".0f")
    out["X12_KR_MS"], out["X12_VER_MS"] = num(d.krawczyk_ms.median(), ".0f"), num(d.verified_ms.median(), ".0f")
    out["X12_KR_CLOSES"] = num(d.krawczyk_closes.astype(bool).sum(), "d")
    out["X12_KR_WC_LO"], out["X12_KR_WC_HI"] = num(both.krawczyk_width_over_certified_inner_median.min(), ".3f"), num(both.krawczyk_width_over_certified_inner_median.max(), ".3f")
    # the three-phase feeder
    t = _x12(h, True)
    tv = t[t.verified.astype(bool)]
    out["X12_3P_N"], out["X12_3P_VER"], out["X12_3P_UNI"] = num(len(t), "d"), num(t.verified.astype(bool).sum(), "d"), num(t.unique.astype(bool).sum(), "d")
    out["X12_3P_WC_LO"], out["X12_3P_WC_HI"] = (num(tv.width_over_certified_inner_median.min(), ".4f"), num(tv.width_over_certified_inner_median.max(), ".3f")) if len(tv) else ("–", "–")
    tk = t[t.krawczyk_closes.astype(bool)]
    out["X12_3P_KR_LO"], out["X12_3P_KR_HI"] = (num(tk.krawczyk_width_over_certified_inner_median.min(), ".3f"), num(tk.krawczyk_width_over_certified_inner_median.max(), ".3f")) if len(tk) else ("–", "–")
    e = _exact(h)
    out["X12_EX_N"], out["X12_EX_GAP"] = num(len(e), "d"), num(e.bracket_gap_max_pct.max(), ".3f")
    out["X12_EX_BOX_LO"], out["X12_EX_BOX_HI"], out["X12_EX_SEC"] = num(e.enclosures.min(), "d"), num(e.enclosures.max(), "d"), num(e.seconds.max(), ".0f")
    lo, bb = e[e.kind.str.startswith("Demanda")], e[~e.kind.str.startswith("Demanda")]
    out["X12_EX_LOAD_LO"], out["X12_EX_LOAD_HI"] = num(lo.enclosure_over_exact_at_most_median.min(), ".4f"), num(lo.enclosure_over_exact_at_most_median.max(), ".4f")
    out["X12_EX_BB_LO"], out["X12_EX_BB_HI"] = num(bb.enclosure_over_exact_at_most_median.min(), ".4f"), num(bb.enclosure_over_exact_at_most_median.max(), ".4f")
    s = _sub(h)
    out["X12_SUB_N"], out["X12_SUB_D_FAIL"] = num(len(s), "d"), num((~s.discs_reached_target.astype(bool)).sum(), "d")
    out["X12_SUB_V"] = _span(h, s.verified_enclosures, "d")
    kok = s[s.krawczyk_reached_target.astype(bool)]
    out["X12_SUB_K"] = _span(h, kok.krawczyk_enclosures, "d")
    out["X12_SUB_K_FAIL"] = num(len(s) - len(kok), "d")
    out["X12_SUB_K_END"] = num(100 * s[~s.krawczyk_reached_target.astype(bool)].krawczyk_largest_excess.max(), ".1f") if len(kok) < len(s) else "–"
    ok = s[s.discs_reached_target.astype(bool)]
    out["X12_SUB_D_LO"], out["X12_SUB_D_HI"] = (num(ok.discs_enclosures.min(), "d"), num(ok.discs_enclosures.max(), "d")) if len(ok) else ("–", "–")
    for tag, key, spec in (("D", "discs", ".2f"), ("K", "krawczyk", ".3f"), ("V", "verified", ".4f")):
        col = s[f"{key}_root_over_exact_median"].dropna()
        out[f"X12_SUB_{tag}_ROOT_LO"], out[f"X12_SUB_{tag}_ROOT_HI"] = num(col.min(), spec), num(col.max(), spec)
    return out


# ----------------------------------------------------------------------------- X9: Krawczyk
_FORM = {"power mismatch": "desbalance de potencia", "current injection": "inyección de corriente"}
_MAT = {"derivative": "derivadas", "slope": "pendientes"}


def _x9(h):
    d = pd.read_csv(h.RES / "x9_krawczyk.csv", dtype={"network": str})
    return d[d.status == "ok"].copy()


def t_x9_variants(h):
    d = _x9(h)
    rows = []
    for (f, mat, rot), g in d.groupby(["formulation", "matrix", "rotated"], sort=False):
        c = g[g.half_disk_closes.astype(bool)]
        rows.append({"f": _FORM.get(f, f), "mat": _MAT.get(mat, mat), "rot": "rotados" if bool(rot) else "globales", "n": len(g),
                     "closes": int(g.half_disk_closes.astype(bool).sum()), "box": int(g.box_closes.astype(bool).sum()),
                     "viol": int(c.half_disk_violations.sum()) if len(c) else 0,
                     "w": _span(h, c.half_disk_width_over_reference_median, ".3f") if len(c) else "–",
                     "ms": h.num(g.krawczyk_ms.median(), ".0f")})
    return h.table(pd.DataFrame(rows), [("f", "Formulación", None), ("mat", "Matriz", None), ("rot", "Ejes", None), ("n", "Casos", "d"),
                                        ("closes", "Cierra (semidiscos)", "d"), ("box", "Cierra (cajas)", "d"), ("viol", "Violaciones", "d"),
                                        ("w", "Ancho / rango alcanzado, mediana", None), ("ms", "Tiempo, mediana (ms)", None)])


def x9_values(h):
    num, out = h.num, {}
    d = _x9(h)
    key = lambda f, mat: d[d.formulation.str.startswith(f) & (d.matrix == mat)]
    out["X9_N"] = num(d.groupby(["network", "scenario", "rho_U"]).ngroups, "d")
    out["X9_POWER"] = num(key("power", "derivative").half_disk_closes.astype(bool).sum() + key("power", "slope").half_disk_closes.astype(bool).sum(), "d")
    for tag, mat in (("DER", "derivative"), ("SLOPE", "slope")):
        g = key("current", mat)
        out[f"X9_CUR_{tag}"] = num(g[g.rotated.astype(bool)].half_disk_closes.astype(bool).sum(), "d")
    out["X9_VIOL"] = num(d[d.half_disk_closes.astype(bool)].half_disk_violations.sum() + d[d.box_closes.astype(bool)].box_violations.sum(), "d")
    b = pd.read_csv(h.RES / "x9_krawczyk_load_boxes.csv", dtype={"network": str})
    best = b[b.formulation.str.startswith("current") & (b.matrix == "slope") & b.rotated.astype(bool)]
    out["X9_LB_N"], out["X9_LB_CLOSES"] = num(len(best), "d"), num(best.krawczyk_closes.astype(bool).sum(), "d")
    out["X9_LB_LO"], out["X9_LB_HI"] = num(best.krawczyk_width_over_reference_median.min(), ".3f"), num(best.krawczyk_width_over_reference_median.max(), ".3f")
    out["X9_LB_POWER"] = num(b[b.formulation.str.startswith("power")].krawczyk_closes.astype(bool).sum(), "d")
    return out


# ----------------------------------------------------------------------------- X10: three phases
def t_x10_validation(h):
    d = pd.read_csv(h.RES / "x10_model_validation.csv")
    d = d[d.status == "ok"].copy()
    d["diff"] = d.max_abs_voltage_difference_pu.map(h._sci)
    d["res"] = d.magnitude_difference_after_removing_that_residual_pu.map(h._sci)
    return h.table(d, [("snapshot", "Instantánea de carga", None), ("model", "Modelo", None), ("buses", "Barras", "d"), ("node_phases", "Nodos-fase", "d"),
                       ("total_load_kva", "Carga (kVA)", ".1f"), ("diff", "Diferencia máx. con pandapower (p.u.)", None),
                       ("res", "ídem, descontado el residuo de pandapower", None), ("phase_spread_max_pu", "Diferencia entre fases, máx. (p.u.)", ".3f"),
                       ("negative_sequence_over_positive_max_pct", "Secuencia negativa / positiva, máx. (%)", ".2f")])


def t_x10_cases(h):
    d = pd.read_csv(h.RES / "x10_three_phase.csv")
    d = d[d.status == "ok"].copy()
    d["pct"] = 100 * d.rho_U
    d["fp"], d["dok"] = d.enclosure_fixed_point.map(lambda v: "sí" if v else "no"), d.discs_ok.map(lambda v: "sí" if v else "no")
    return h.table(d, [("scenario", "Flujo", None), ("pct", "ρU (%)", ".0f"), ("attained_width_max_pu", "Rango alcanzado, máx. (p.u.)", ".4f"),
                       ("fp", "Cierra", None), ("enclosure_width_over_attained_median", "Encierro / rango, mediana", ".4f"),
                       ("enclosure_width_over_attained_max", "máx.", ".4f"), ("enclosure_violations", "Violaciones", "d"),
                       ("dok", "Discos: cierra", None), ("discs_per_bus_width_over_attained_median", "Discos por nodo / rango, mediana", ".2f"),
                       ("node_phases_with_minimum_not_at_zero_p_minus_C", "Nodos-fase con el mínimo fuera de (0, −C)", "d")])


def t_x10_limits(h):
    d = pd.read_csv(h.RES / "x10_limits.csv")
    d["fp"] = d.enclosure_fixed_point.map(lambda v: "sí" if v else "no")
    d["dok"] = d.discs_ok.map(lambda v: "sí" if v else "no")
    d["v"] = [f"{h.num(a, '.3f')} a {h.num(b, '.3f')}" for a, b in zip(d.attained_vmin_pu, d.attained_vmax_pu)]
    for c in ("enclosure_width_over_attained_median", "discs_per_bus_width_over_attained_median"):
        if c not in d:
            d[c] = np.nan
    return h.table(d, [("total_rating_over_load", "Capacidad total / carga", ".0f"), ("v", "\\|V\\| alcanzado (p.u.)", None), ("fp", "Encierro: cierra", None),
                       ("enclosure_width_over_attained_median", "Encierro / rango, mediana", ".3f"), ("dok", "Discos: cierra", None),
                       ("discs_per_bus_width_over_attained_median", "Discos por nodo / rango, mediana", ".2f"), ("ac_failures_in_samples", "Muestras sin solución AC", "d")])


def x10_values(h):
    num, out = h.num, {}
    v = pd.read_csv(h.RES / "x10_model_validation.csv")
    v = v[v.status == "ok"]
    out["X10_VAL_DIFF"] = h._sci(v.max_abs_voltage_difference_pu.max())
    out["X10_VAL_RES"] = h._sci(v.magnitude_difference_after_removing_that_residual_pu.max())
    out["X10_SPREAD"] = num(v.phase_spread_max_pu.max(), ".3f")
    out["X10_NEG"], out["X10_ZERO"] = num(v.negative_sequence_over_positive_max_pct.max(), ".2f"), num(v.zero_sequence_over_positive_max_pct.max(), ".2f")
    out["X10_TRAFO"] = h._sci(v.leaving_the_transformer_out_changes_magnitudes_by_max_pu.max())
    d = pd.read_csv(h.RES / "x10_three_phase.csv")
    d = d[d.status == "ok"]
    out["X10_N"], out["X10_FP"] = num(len(d), "d"), num(d.enclosure_fixed_point.astype(bool).sum(), "d")
    out["X10_BUSES"], out["X10_NODES"], out["X10_K"], out["X10_SRC"] = num(d.buses.iloc[0], "d"), num(d.pq_node_phases.iloc[0], "d"), num(d.symbols.iloc[0], "d"), num(d.sources.iloc[0], "d")
    out["X10_PHASES"] = str(d.sources_on_phases_a_b_c.iloc[0])
    out["X10_W_LO"], out["X10_W_HI"] = num(d.enclosure_width_over_attained_median.min(), ".3f"), num(d.enclosure_width_over_attained_median.max(), ".3f")
    out["X10_D_LO"], out["X10_D_HI"] = num(d.discs_per_bus_width_over_attained_median.min(), ".2f"), num(d.discs_per_bus_width_over_attained_median.max(), ".2f")
    out["X10_VIOL"] = num(d.enclosure_violations.sum() + d.discs_per_bus_violations.sum(), "d")
    out["X10_MONO"] = num((d.node_phases_with_minimum_not_at_zero_p_minus_C == d.pq_node_phases).sum(), "d")
    out["X10_MONO_GAP"] = _span(h, d.largest_gap_over_attained_width_pct, ".0f")
    out["X10_SO2"] = num(((d.second_order_only_upper_violations > 0) | (d.second_order_only_lower_violations > 0)).sum(), "d")
    f = pd.read_csv(h.RES / "x10_full_model.csv")
    out["X10_FULL_BUSES"], out["X10_FULL_NODES"] = num(f.buses.iloc[0], "d"), num(f.pq_node_phases.iloc[0], "d")
    out["X10_FULL_FP"], out["X10_FULL_N"] = num(f.enclosure_fixed_point.astype(bool).sum(), "d"), num(len(f), "d")
    out["X10_FULL_S_LO"], out["X10_FULL_S_HI"] = num(f.total_ms.min() / 1e3, ".1f"), num(f.total_ms.max() / 1e3, ".1f")
    out["X10_FULL_DIFF"] = h._sci(f.bounds_minus_merged_model_max_pu.max())
    out["X10_FULL_VIOL"] = num(f.violations_in_samples.sum(), "d")
    lim = pd.read_csv(h.RES / "x10_limits.csv")
    ok = lim[lim.enclosure_fixed_point.astype(bool)]
    out["X10_LIM_LAST"] = num(ok.total_rating_over_load.max(), ".0f")
    out["X10_LIM_FIRST_FAIL"] = num(lim[~lim.enclosure_fixed_point.astype(bool)].total_rating_over_load.min(), ".0f")
    out["X10_LIM_DISCS_LAST"] = num(lim[lim.discs_ok.astype(bool)].total_rating_over_load.max(), ".0f")
    return out


# ----------------------------------------------------------------------------- X13: scale
_DESIGN = {"every tenth bus": "una cada 10 barras", "8 sources": "8", "32 sources": "32"}


def _x13(h):
    return pd.read_csv(h.RES / "x13_scale.csv", dtype={"network": str})


def t_x13(h):
    d = _x13(h)
    rows = []
    for (net, des), g in d.groupby(["network", "design"], sort=False):
        ok = g[g.status == "ok"]
        row = {"network": net, "bus": int(g.n_bus.iloc[0]), "des": _DESIGN.get(des, des), "K": int(g.symbols.iloc[0]), "cases": len(g)}
        if len(ok):
            row.update({"ver": f"{int(ok.verified.astype(bool).sum())} / {int(ok.unique.astype(bool).sum())}", "model": h.num(ok.model_s.median(), ".2f"),
                        "tube": h.num(ok.tube_s.median(), ".2f"), "term": h.num(ok.magnitudes_termwise_s.median(), ".2f"),
                        "full": h.num(ok.magnitudes_full_s.median(), ".1f"), "wterm": h.num(ok.width_over_attained_termwise_median.max(), ".3f"),
                        "wk3": h.num(ok.width_over_attained_k3_median.max(), ".3f"), "wfull": h.num(ok.width_over_attained_full_median.max(), ".3f"),
                        "rss": h.num(ok.peak_rss_gb.max(), ".1f")})
        else:
            row.update({k: "–" for k in ("model", "tube", "term", "full", "wterm", "wk3", "wfull", "rss")})
            row["ver"] = "no se intentó"
        rows.append(row)
    return h.table(pd.DataFrame(rows), [("network", "Red", None), ("bus", "Barras", "d"), ("des", "Fuentes", None), ("K", "K", "d"), ("cases", "Casos", "d"),
                                        ("ver", "Verificado / con unicidad", None), ("model", "Modelo (s)", None), ("tube", "Tubo (s)", None),
                                        ("term", "Módulos, término a término (s)", None), ("full", "Módulos, RAA-K completa (s)", None),
                                        ("wterm", "Ancho / rango alcanzado, peor caso: término a término", None), ("wk3", "bloques de hasta 3", None),
                                        ("wfull", "RAA-K completa", None), ("rss", "Memoria pico (GB)", None)])


def x13_values(h):
    num, out = h.num, {}
    d = _x13(h)
    ok = d[d.status == "ok"]
    out["X13_N"], out["X13_OK"] = num(len(d), "d"), num(len(ok), "d")
    out["X13_VER"], out["X13_UNI"] = num(ok.verified.astype(bool).sum(), "d"), num(ok.unique.astype(bool).sum(), "d")
    out["X13_BUS_MAX"] = f"{int(ok.n_bus.max()):,}".replace(",", ".")
    out["X13_K_AT_BUS_MAX"] = num(ok[ok.n_bus == ok.n_bus.max()].symbols.max(), "d")
    out["X13_K_MAX"], out["X13_BUS_AT_K_MAX"] = num(ok.symbols.max(), "d"), num(ok[ok.symbols == ok.symbols.max()].n_bus.max(), "d")
    bad = d[d.status != "ok"]
    out["X13_BLOCKED_N"] = num(len(bad), "d")
    out["X13_BLOCKED"] = "; ".join(f"{net}, K = {int(g.symbols.iloc[0])}: {g.status.iloc[0]}" for net, g in bad.groupby("network", sort=False)) or "ninguno"
    out["X13_BLOCKED_K"] = num(bad.symbols.max(), "d") if len(bad) else "–"
    out["X13_BLOCKED_GB"] = num(bad.tensor_gb.max(), ".0f") if len(bad) and "tensor_gb" in bad else "–"
    big = ok[ok.n_bus == ok.n_bus.max()]
    out["X13_TUBE_BIG_LO"], out["X13_TUBE_BIG_HI"] = num(big.tube_s.min(), ".0f"), num(big.tube_s.max(), ".0f")
    out["X13_RSS_MAX"] = num(ok.peak_rss_gb.max(), ".1f")
    out["X13_FULL_MAX"], out["X13_TERM_AT_FULL_MAX"] = num(ok.magnitudes_full_s.max(), ".0f"), num(ok.loc[ok.magnitudes_full_s.idxmax(), "magnitudes_termwise_s"], ".1f")
    out["X13_VIOL"] = num(ok.violations_full.sum() + ok.violations_k3.sum() + ok.violations_termwise.sum(), "d")
    ob = ok[(ok.network == "Oberrhein") & (ok.design == "every tenth bus")]
    base = ob.model_s.median() + ob.tube_s.median()
    out["X13_OB_K"] = num(ob.symbols.iloc[0], "d")
    out["X13_OB_TERM"], out["X13_OB_FULL"] = num(base + ob.magnitudes_termwise_s.median(), ".1f"), num(base + ob.magnitudes_full_s.median(), ".1f")
    out["X13_FULL_GAIN"] = _span(h, ok.width_full_over_termwise_median, ".4f")
    out["X13_K3_GAIN"] = _span(h, ok.width_k3_over_termwise_median, ".4f")
    out["X13_K3_AT_FULL_MAX"] = num(ok.loc[ok.magnitudes_full_s.idxmax(), "magnitudes_k3_s"], ".1f")
    out["X13_W_FULL_LO"], out["X13_W_FULL_HI"] = num(ok.width_over_attained_full_median.min(), ".3f"), num(ok.width_over_attained_full_median.max(), ".3f")
    out["X13_FLOOR"] = h._sci(ok.tube_radius_without_uncertainty_max_pu.max())
    out["X13_FLOOR_SHARE"] = num(100 * ok.floor_share_at_largest_radius.max(), ".1f")
    return out


# ----------------------------------------------------------------------------- X14: PV buses
def t_x14_distribution(h):
    d = pd.read_csv(h.RES / "x14_pv_distribution.csv", dtype={"network": str})
    d = d[d.status == "ok"]
    rows = []
    for (net, scen), g in d.groupby(["network", "scenario"], sort=False):
        g = g.sort_values("rho_U")
        rows.append({"network": net, "scen": scen, "pv": int(g.pv_buses.iloc[0]), "K": int(g.symbols.iloc[0]),
                     "ver": f"{int(g.verified.astype(bool).sum())} / {int(g.unique.astype(bool).sum())} de {len(g)}",
                     "g": _span(h, g.linear_radius_G, ".2f"), "w": " · ".join(h.num(x, ".3f") for x in g.width_over_vertex_range_median),
                     "wmax": " · ".join(h.num(x, ".2f") for x in g.width_over_vertex_range_max), "viol": int(g.violations.sum())})
    return h.table(pd.DataFrame(rows), [("network", "Red", None), ("scen", "Flujo", None), ("pv", "Barras PV", "d"), ("K", "K", "d"),
                                        ("ver", "Verificado / con unicidad", None), ("g", "Radio de la parte lineal", None),
                                        ("w", "Encierro / rango exacto en vértices, mediana, con ρU = 5 · 10 · 20 · 50 · 80 %", None),
                                        ("wmax", "peor barra", None), ("viol", "Violaciones", "d")])


def t_x14_transmission(h):
    d = pd.read_csv(h.RES / "x14_pv_transmission.csv")
    d = d[d.status == "ok"].copy()
    d["pct"] = 100 * d.delta
    d["ver"] = [("sí" if v else "no") + (" (única)" if u else "") for v, u in zip(d.verified.astype(bool), d.unique.astype(bool))]
    for c in ("width_over_vertex_range_median", "width_over_vertex_range_max"):
        if c not in d:
            d[c] = np.nan
    return h.table(d, [("network", "Red", None), ("pv_buses", "Barras PV", "d"), ("pq_buses", "Barras PQ", "d"), ("pct", "δ (%)", ".0f"), ("ver", "Verificado", None),
                       ("linear_radius_G", "Radio de la parte lineal", ".3f"), ("majorant_radius", "Radio del mayorante de la forma V-I", ".2f"),
                       ("width_over_vertex_range_median", "Encierro / rango exacto en vértices, mediana", ".3f"), ("width_over_vertex_range_max", "peor barra", ".2f")])


def t_x14_forms(h):
    d = pd.read_csv(h.RES / "x14_forms.csv", dtype={"network": str})
    d = d[d.status == "ok"]
    rows = []
    for net, g in d.groupby("network", sort=False):
        both = g[g.current_form_verified.astype(bool) & g.vi_form_verified.astype(bool)]
        rows.append({"network": net, "cases": len(g), "y": h._sci(1.0 / g.largest_admittance_pu.iloc[0]), "cur": int(g.current_form_verified.astype(bool).sum()),
                     "vi": int(g.vi_form_verified.astype(bool).sum()), "maj": h.num(g.majorant_radius.max(), ".2f"),
                     "majnr": h.num(g.majorant_radius_without_rounding_terms.max(), ".2f"),
                     "ratio": _span(h, both.vi_over_current_radius_median, ".2f") if len(both) else "–"})
    return h.table(pd.DataFrame(rows), [("network", "Red", None), ("cases", "Casos", "d"), ("y", "Rama más corta, \\|z\\| (p.u.)", None),
                                        ("cur", "Verifica la forma de corriente", "d"), ("vi", "Verifica la forma V-I", "d"),
                                        ("maj", "Radio del mayorante V-I, máx.", None), ("majnr", "ídem sin los términos de redondeo", None),
                                        ("ratio", "Radio del tubo V-I / corriente, mediana", None)])


def x14_values(h):
    num, out = h.num, {}
    d = pd.read_csv(h.RES / "x14_pv_distribution.csv", dtype={"network": str})
    d = d[d.status == "ok"]
    out["X14_D_N"], out["X14_D_VER"], out["X14_D_UNI"] = num(len(d), "d"), num(d.verified.astype(bool).sum(), "d"), num(d.unique.astype(bool).sum(), "d")
    out["X14_D_PV"], out["X14_D_K"] = num(d.pv_buses.iloc[0], "d"), num(d.symbols.iloc[0], "d")
    small, big = d[d.rho_U <= 0.2], d[d.rho_U > 0.2]
    out["X14_D_W20_LO"], out["X14_D_W20_HI"] = num(small.width_over_vertex_range_median.min(), ".3f"), num(small.width_over_vertex_range_median.max(), ".3f")
    out["X14_D_WBIG_LO"], out["X14_D_WBIG_HI"] = num(big.width_over_vertex_range_median.min(), ".2f"), num(big.width_over_vertex_range_median.max(), ".2f")
    out["X14_D_WORST"] = num(d.width_over_vertex_range_max.max(), ".1f")
    out["X14_D_VIOL"] = num(d.violations.sum(), "d")
    t = pd.read_csv(h.RES / "x14_pv_transmission.csv")
    t = t[t.status == "ok"]
    ok = t[t.verified.astype(bool)]
    out["X14_T_N"], out["X14_T_VER"] = num(len(t), "d"), num(len(ok), "d")
    out["X14_T_NETS_OK"] = ", ".join(dict.fromkeys(ok.network)) or "ninguna"
    out["X14_T_NETS_BAD"] = " y ".join(n for n in dict.fromkeys(t.network) if n not in set(ok.network)) or "ninguna"
    out["X14_T_W_LO"], out["X14_T_W_HI"] = (num(ok.width_over_vertex_range_median.min(), ".3f"), num(ok.width_over_vertex_range_median.max(), ".2f")) if len(ok) else ("–", "–")
    for net in ("case118", "case300"):
        g = t[t.network == net]
        tag = net[4:]
        if len(g):
            out[f"X14_T_G{tag}_LO"], out[f"X14_T_G{tag}_HI"] = num(g.linear_radius_G.min(), ".2f"), num(g.linear_radius_G.max(), ".2f")
            out[f"X14_T_M{tag}_LO"] = num(g.majorant_radius.min(), ".1f")
    f = pd.read_csv(h.RES / "x14_forms.csv", dtype={"network": str})
    f = f[f.status == "ok"]
    out["X14_F_N"], out["X14_F_CUR"], out["X14_F_VI"] = num(len(f), "d"), num(f.current_form_verified.astype(bool).sum(), "d"), num(f.vi_form_verified.astype(bool).sum(), "d")
    only = f[f.current_form_verified.astype(bool) & ~f.vi_form_verified.astype(bool)]
    out["X14_F_ONLY_CUR"] = ", ".join(f"{n} ({len(g)})" for n, g in only.groupby("network", sort=False)) or "ninguno"
    out["X14_F_ONLY_VI"] = num((~f.current_form_verified.astype(bool) & f.vi_form_verified.astype(bool)).sum(), "d")
    both = f[f.current_form_verified.astype(bool) & f.vi_form_verified.astype(bool)]
    out["X14_F_RATIO_LO"], out["X14_F_RATIO_HI"] = num(both.vi_over_current_radius_median.min(), ".2f"), num(both.vi_over_current_radius_median.max(), ".2f")
    return out


def tables(h):
    return {"x12_verified": lambda: t_x12_verified(h), "x12_krawczyk": lambda: t_x12_krawczyk(h), "x12_exact": lambda: t_x12_exact(h),
            "x12_subdivision": lambda: t_x12_subdivision(h),
            "x9_variants": lambda: t_x9_variants(h), "x10_validation": lambda: t_x10_validation(h), "x10_cases": lambda: t_x10_cases(h),
            "x10_limits": lambda: t_x10_limits(h), "x13": lambda: t_x13(h), "x14_distribution": lambda: t_x14_distribution(h),
            "x14_transmission": lambda: t_x14_transmission(h), "x14_forms": lambda: t_x14_forms(h)}


def values(h):
    return {**x12_values(h), **x9_values(h), **x10_values(h), **x13_values(h), **x14_values(h)}
