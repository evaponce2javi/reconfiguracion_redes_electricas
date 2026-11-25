"""Figures with every uncertainty set drawn under its own name (corrected model).

  fig_gain_vs_rho_<system>.png   G against |rho| for the three boxes that can sit inside G
  fig_sets_<system>.png          one branch: AC samples, covariance ellipse, the interval the method
                                 returns (lambda_max), its exact projection, the 3-sigma box, the affine box
  fig_pp_calibration.png         P-P calibration of the ellipse, with and without the 1e-9 jitter

Writes PNG files to results/figures/. Colours: first slots of a colour-blind-validated categorical
palette, always paired with a line style so that no set is identified by colour alone.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Rectangle
from scipy import stats

from common import CFG, method, seed_for, study
from cdhim import geometry as geo
from cdhim import sample_eps
from cdhim.config import RESULTS_DIR

BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, GRID, SAMPLE = "#0b0b0b", "#52514e", "#e3e2de", "#9a9992"
OUT = RESULTS_DIR / "figures"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 300})


def paper_branch(res):
    ok = np.flatnonzero(res.valid & (np.abs(res.rho) < 1 - res.cfg.delta_stab))
    return ok[np.argmax(res.gain[ok])]


def fig_gain_vs_rho(key, st, delta):
    fig, ax = plt.subplots(figsize=(5.2, 4.3))
    rho = np.linspace(0, 1, 400)
    chi2 = geo.chi2_quantile(CFG["method"]["confidence"])
    series = (("covariance_1sigma", BLUE, "-", r"covariance box, same scale ($\kappa=\pi/4$)"),
              ("covariance_3sigma", ORANGE, "--", r"3$\sigma$ box vs. $\chi^2$ ellipse ($\kappa=1.03$)"))
    for box, col, ls, label in series:
        ax.plot(rho, 100 * geo.gain_closed_form(rho, geo.kappa_constant(box, chi2)), color=col, ls=ls, lw=1.8, label=label)
    res = st.run(delta, method(gain_box="affine"))
    v = res.valid
    ax.scatter(np.abs(res.rho[v]), 100 * res.gain[v], s=7, color=AQUA, alpha=0.55, linewidths=0, zorder=3,
               label=r"affine box $\sum|c_k|$ ($\kappa$ varies by branch)")
    g = 100 * CFG["method"]["gamma_th"]
    ax.axhline(g, color=INK2, lw=0.9, ls=":")
    ax.text(0.01, g + 1.5, rf"$\gamma_{{th}}$ = {g:.0f} %", color=INK2, fontsize=8)
    ax.axvline(1 - CFG["method"]["delta_stab"], color=INK2, lw=0.9, ls=":")
    ax.set(xlabel=r"correlation magnitude $|\rho|$", ylabel="coupling gain G (%)", xlim=(0, 1.005), ylim=(0, 102),
           title=f"{CFG['systems'][key]['label']}: G depends on which box it is measured against")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=7.5, ncol=1)
    fig.tight_layout()
    fig.savefig(OUT / f"fig_gain_vs_rho_{key}.png")
    plt.close(fig)


def fig_sets(key, st, delta, n=5000):
    res = st.run(delta, method())
    b = paper_branch(res)
    i, j = res.fi[b], res.ti[b]
    eps = sample_eps(np.random.default_rng(seed_for("fig_sets", key)), st.src.n, n)
    V, conv, _ = st.solver.solve(st.injections(eps, delta))
    dv = (np.abs(V[st.grid.pq][:, conv]) - st.vm_star[:, None])[[i, j]]
    Q = np.array([[res.q11[b], res.q12[b]], [res.q12[b], res.q22[b]]])
    lam, vec = np.linalg.eigh(Q)
    ang = np.degrees(np.arctan2(vec[1, 1], vec[0, 1]))
    s = np.sqrt(res.chi2)
    fig, ax = plt.subplots(figsize=(5.2, 6.4))
    ax.scatter(dv[0], dv[1], s=2, color=SAMPLE, linewidths=0, rasterized=True, label=f"AC power-flow samples (n = {dv.shape[1]})")
    ax.add_patch(Ellipse((0, 0), 2 * s * np.sqrt(lam[1]), 2 * s * np.sqrt(lam[0]), angle=ang, fill=False, ec=BLUE, lw=1.8,
                         label=r"covariance ellipse, $\chi^2_{2,0.9973}$"))
    boxes = ((res.r_lmax[b], res.r_lmax[b], ORANGE, "-", r"returned interval, $\sqrt{\chi^2\lambda_{max}}$ (Eq. 11)"),
             (res.r_coord[i], res.r_coord[j], AQUA, "--", r"exact projection, $\sqrt{\chi^2 Q_{ii}}$"),
             (res.r_cltm[i], res.r_cltm[j], YELLOW, ":", r"3$\sigma$ box ('CLTM box' of the manuscript)"),
             (res.r_affine[i], res.r_affine[j], MAGENTA, "-.", r"affine box, $\sum|c_k|$"))
    for rx, ry, col, ls, label in boxes:
        ax.add_patch(Rectangle((-rx, -ry), 2 * rx, 2 * ry, fill=False, ec=col, lw=1.6, ls=ls, label=label))
    lim = 1.08 * max(res.r_affine[i], res.r_affine[j], res.r_lmax[b])
    lab = st.grid.bus_label
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal",
           xlabel=rf"$\Delta V$ at bus {lab[st.grid.pq[i]]} (p.u.)", ylabel=rf"$\Delta V$ at bus {lab[st.grid.pq[j]]} (p.u.)",
           title=f"{CFG['systems'][key]['label']}, $\\delta$ = {100 * delta:.0f} %, $\\rho$ = {res.rho[b]:.4f}")
    ax.ticklabel_format(style="sci", scilimits=(-2, 2))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=7.5, ncol=2, handlelength=2.6)
    fig.tight_layout()
    fig.savefig(OUT / f"fig_sets_{key}.png")
    plt.close(fig)


def fig_pp(studies):
    levels = np.array([0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 0.9973])
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.plot([0.45, 1], [0.45, 1], color=INK2, lw=0.9, ls=":", label="ideal (coverage = confidence)")
    styles = {"ieee33": (BLUE, "-", "o"), "pegase2869": (ORANGE, "--", "s")}
    for key, (st, delta) in studies.items():
        res = st.run(delta, method())
        b = paper_branch(res)
        C = st.coefficients(delta)[[res.fi[b], res.ti[b]]]
        eps = sample_eps(np.random.default_rng(seed_for("fig_pp", key)), st.src.n, 100_000)
        pts = C @ eps
        Q = C @ C.T / 3.0
        col, ls, mk = styles[key]
        d2 = np.sum(pts * np.linalg.solve(Q, pts), axis=0)
        ax.plot(levels, [np.mean(d2 <= stats.chi2.ppf(p, 2)) for p in levels], color=col, ls=ls, marker=mk, ms=4, lw=1.6,
                label=f"{CFG['systems'][key]['label']}, exact distance")
        if key == "pegase2869":
            d2j = np.sum(pts * np.linalg.solve(Q + 1e-9 * np.eye(2), pts), axis=0)
            ax.plot(levels, [np.mean(d2j <= stats.chi2.ppf(p, 2)) for p in levels], color=AQUA, ls="-.", marker="^", ms=4, lw=1.6,
                    label="PEGASE-2869, with the 1e-9 jitter\nof the original code")
    ax.set(xlabel="theoretical confidence", ylabel="empirical coverage (linear model)", xlim=(0.45, 1.01), ylim=(0.40, 1.01),
           title="Ellipse calibration: the jitter creates the margin")
    ax.legend(loc="lower right", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_pp_calibration.png")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    studies = {}
    for key, sysc in CFG["systems"].items():
        st = study(sysc["case"])
        studies[key] = (st, sysc["delta"])
        fig_gain_vs_rho(key, st, sysc["delta"])
        fig_sets(key, st, sysc["delta"])
    fig_pp(studies)
    print(OUT)


if __name__ == "__main__":
    main()
