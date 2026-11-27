"""Figures of the affine-forms paper (paper/figures/*.pdf), from the CSV tables in results/.

    python experiments/figures_paper.py compute   # Fig. 1 data: samples and shape of the quadratic term (writes CSVs)
    python experiments/figures_paper.py           # draw the three figures from the CSVs

Fig. 1  Quadratic term of a PEGASE-2869 branch drop against the normal law with its mean and variance, and the
        skewness of the quadratic term over all bus and drop targets   (results/fig1_*.csv)
Fig. 2  Worst AC violation over load x dispersion on PEGASE-2869: 3 sigma, moment correction, proposed bound
        (results/b4_sweep.csv, results/b5_sweep.csv)
Fig. 3  Worst AC violation of branch-current bounds per system and method   (results/b9_branch_currents.csv)

Colours: the first slots of the reference categorical palette (validated for colour-vision deficiency); bars carry
hatching and value labels so that the figures read in grayscale.
"""
import sys

import numpy as np
import pandas as pd

from common import CFG, ROOT, method, seed_for, study
from cdhim import save_table

RES = ROOT / "results"
FIG = ROOT / "paper" / "figures"
BLUE, ORANGE, AQUA, YELLOW = "#257de0", "#f47230", "#00b894", "#f2a600"
TEXT, MUTED, GRID = "#191919", "#666666", "#e6e6e6"
ALPHA = 1.0 - CFG["method"]["confidence"]


# ----------------------------------------------------------------------------- Fig. 1 data
def compute_fig1(n_stats=4000, n_example=20000):
    from cdhim import SecondOrder, sample_eps
    from cdhim.affine_forms import NonlinearSampler, bus_targets, drop_targets, sigma
    key = "pegase2869"
    sysc = CFG["systems"][key]
    delta = sysc["delta"]
    seed = seed_for("fig1", key)
    st = study(sysc["case"])
    res = st.run(delta, method())
    so = SecondOrder(st, delta)
    v = res.valid
    n_pq = st.C1.shape[0]
    Ws = {"bus": bus_targets(n_pq), "drop": drop_targets(n_pq, res.fi[v], res.ti[v])}
    ns = NonlinearSampler(so, Ws)
    rng = np.random.default_rng(seed)
    Q = {k: [] for k in Ws}
    for s in range(0, n_stats, 500):
        for k, q in ns.batch(sample_eps(rng, st.src.n, 500)).items():
            Q[k].append(q)
    rows = []
    for k, W in Ws.items():
        q = np.hstack(Q[k])
        sl = sigma(np.asarray(W.T @ (st.C1 * delta)))
        m, sd = q.mean(axis=1), q.std(axis=1)
        z = (q - m[:, None]) / np.where(sd > 0, sd, 1)[:, None]
        for t in range(q.shape[0]):
            if sd[t] <= 0 or sl[t] <= 0:
                continue
            rows.append({"system": key, "delta": delta, "seed": seed, "n_samples": n_stats, "target": k, "index": t,
                         "skewness": float((z[t] ** 3).mean()), "excess_kurtosis": float((z[t] ** 4).mean() - 3),
                         "max_over_sd": float(z[t].max()), "min_over_sd": float(z[t].min()),
                         "sd_q_over_sigma_linear": float(sd[t] / sl[t])})
    stats = pd.DataFrame(rows)
    print(save_table("fig1_quadratic_stats", rows, {"note": "shape of the quadratic term per bus and drop target, PEGASE-2869 "
                                                            "design point; skewness and excess kurtosis from n_samples draws"}))
    # example: the drop with the largest |skewness| among those where the quadratic term is not negligible
    cand = stats[(stats.target == "drop") & (stats.sd_q_over_sigma_linear > 0.05)]
    ex = cand.loc[cand.skewness.abs().idxmax()]
    t = int(ex["index"])
    rng = np.random.default_rng(seed + 1)
    qs = []
    for s in range(0, n_example, 500):
        qs.append(ns.batch(sample_eps(rng, st.src.n, 500))["drop"][t])
    q = np.concatenate(qs)
    print(save_table("fig1_example_samples", [{"q": float(x)} for x in q],
                     {"system": key, "target": "drop", "index": t, "branch_from_pq": int(res.fi[v][t]), "branch_to_pq": int(res.ti[v][t]),
                      "n_samples": n_example, "seed": seed + 1, "note": "samples of the quadratic term of the example drop"}))


# ----------------------------------------------------------------------------- drawing
def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": TEXT, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                         "axes.titlesize": 8, "legend.frameon": False, "lines.linewidth": 1.5})
    FIG.mkdir(parents=True, exist_ok=True)
    return plt


def fig1(plt):
    from scipy.stats import norm
    q = pd.read_csv(RES / "fig1_example_samples.csv").q.to_numpy()
    st = pd.read_csv(RES / "fig1_quadratic_stats.csv")
    m, sd = q.mean(), q.std()
    z = (q - m) / sd
    fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.6), constrained_layout=True)
    a = ax[0]
    a.hist(z, bins=120, density=True, color=BLUE, alpha=0.85, edgecolor="white", linewidth=0.2, label="quadratic term")
    xs = np.linspace(z.min() - 0.5, z.max() + 0.5, 400)
    a.plot(xs, norm.pdf(xs), color=ORANGE, label="normal, same mean and variance")
    a.axvline(z.max(), color=BLUE, linestyle=":", linewidth=1.2)
    a.axvline(norm.ppf(1 - ALPHA / 2), color=ORANGE, linestyle=":", linewidth=1.2)
    top = a.get_ylim()[1]
    a.set_ylim(0, top * 1.35)
    a.text(z.max(), top * 0.7, f"sampled maximum\n{z.max():.1f} sd ", color=TEXT, fontsize=7, ha="right")
    a.text(norm.ppf(1 - ALPHA / 2), top * 0.4, " normal\n 99.865%:\n 3.0 sd", color=TEXT, fontsize=7)
    a.set_xlabel("standardized value of the quadratic term")
    a.set_ylabel("density")
    a.set_title("(a) one branch drop, PEGASE-2869")
    a.legend(loc="upper right", fontsize=7)
    a.grid(axis="y", color=GRID, linewidth=0.5)
    b = ax[1]
    for tg, col, lab in (("bus", BLUE, "bus voltages"), ("drop", AQUA, "branch drops")):
        s = st[st.target == tg].skewness
        b.hist(s, bins=np.linspace(-6, 6, 61), color=col, alpha=0.7, edgecolor="white", linewidth=0.2,
               label=f"{lab} (median |skew| {s.abs().median():.1f})")
    b.axvline(0, color=MUTED, linewidth=0.8)
    b.set_xlabel("skewness of the quadratic term")
    b.set_ylabel("number of targets")
    b.set_title("(b) all targets, PEGASE-2869")
    b.set_ylim(0, b.get_ylim()[1] * 1.3)
    b.legend(fontsize=7, loc="upper right")
    b.grid(axis="y", color=GRID, linewidth=0.5)
    fig.savefig(FIG / "fig1_quadratic_term.pdf")
    plt.close(fig)


def fig2(plt):
    from matplotlib.colors import LogNorm
    b4 = pd.read_csv(RES / "b4_sweep.csv")
    b5 = pd.read_csv(RES / "b5_sweep.csv")
    b4 = b4[(b4.system == "pegase2869") & (b4.status == "ok")]
    b5 = b5[(b5.system == "pegase2869") & (b5.status == "ok")]
    panels = [(b4, "ac_cltm_3sigma_none_max_pct", "(a) 3$\\sigma$, first order"),
              (b4, "ac_cltm_3sigma_moment_max_pct", "(b) 3$\\sigma$ with mean and variance\nof the quadratic term"),
              (b5, "ac_bus_bern_max_pct", "(c) proposed: Bernstein + range")]
    lfs = sorted(b4.load_factor.unique())
    ds = sorted(b4.delta.unique())
    fig, ax = plt.subplots(1, 3, figsize=(6.3, 2.5), constrained_layout=True, sharey=True)
    norm = LogNorm(vmin=0.005, vmax=10)
    for a, (d, col, title) in zip(ax, panels):
        M = np.full((len(lfs), len(ds)), np.nan)
        for _, r in d.iterrows():
            M[lfs.index(r.load_factor), ds.index(r.delta)] = r[col]
        im = a.imshow(np.maximum(M, 0.005), cmap="Blues", norm=norm, origin="lower", aspect="auto")
        for i in range(len(lfs)):
            for j in range(len(ds)):
                val = M[i, j]
                over = val > 100 * ALPHA
                a.text(j, i, f"{val:.2f}" if val < 10 else f"{val:.0f}", ha="center", va="center", fontsize=6,
                       color="white" if np.maximum(val, 0.005) > 0.5 else TEXT, fontweight="bold" if over else "normal")
                if over:
                    a.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=ORANGE, linewidth=1.2))
        a.set_xticks(range(len(ds)), [f"{100 * x:.0f}" for x in ds])
        a.set_xlabel("dispersion $\\delta$ (%)")
        a.set_title(title)
        a.spines[:].set_visible(False)
    ax[0].set_yticks(range(len(lfs)), [f"{x:.2f}" for x in lfs])
    ax[0].set_ylabel("load factor")
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("worst bus violation (%), log scale")
    cb.outline.set_visible(False)
    fig.savefig(FIG / "fig2_sweep.pdf")
    plt.close(fig)


def fig3(plt):
    d = pd.read_csv(RES / "b9_branch_currents.csv")
    order = ["case33bw", "case118", "case1354pegase", "pglib_opf_case2869_pegase"]
    names = {"case33bw": "IEEE-33", "case118": "IEEE-118", "case1354pegase": "PEGASE-1354", "pglib_opf_case2869_pegase": "PEGASE-2869"}
    meths = [("cumulant_cf", "cumulants + Cornish-Fisher", ORANGE, "//"), ("pem_cf", "point estimate (AC) + Cornish-Fisher", YELLOW, ".."),
             ("bern_linear_only", "Bernstein, first order", AQUA, "\\\\"), ("bern", "proposed: Bernstein + range", BLUE, "")]
    g = d.groupby(["case", "method"]).max_pct.max()
    fig, a = plt.subplots(figsize=(6.3, 2.8), constrained_layout=True)
    w = 0.2
    for k, (m, lab, col, hatch) in enumerate(meths):
        vals = [max(g.loc[(c, m)], 1e-3) for c in order]
        xs = np.arange(len(order)) + (k - 1.5) * w
        a.bar(xs, vals, width=w - 0.02, color=col, hatch=hatch, edgecolor="white", linewidth=0.5, label=lab)
        for x, v in zip(xs, vals):
            lab_v = "0" if v <= 1e-3 else (f"{v:.0f}" if v >= 10 else f"{v:.2g}")
            a.text(x, v * 1.15, lab_v, ha="center", va="bottom", fontsize=6, color=TEXT)
    a.axhline(100 * ALPHA, color=MUTED, linestyle="--", linewidth=0.8, label="nominal level 0.27%")
    a.set_yscale("log")
    a.set_ylim(1e-3, 400)
    a.set_xticks(range(len(order)), [names[c] for c in order])
    a.set_ylabel("worst branch-end violation (%)\nover 3 cells, log scale")
    a.legend(ncol=3, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    a.grid(axis="y", color=GRID, linewidth=0.5)
    fig.savefig(FIG / "fig3_currents.pdf")
    plt.close(fig)


def main():
    if "compute" in sys.argv[1:]:
        compute_fig1()
        return
    plt = setup()
    fig1(plt)
    fig2(plt)
    fig3(plt)
    print("figures in", FIG)


if __name__ == "__main__":
    main()
