"""Figures for the PRD network draft, built only from versioned results:
fig_bracket.pdf: the scaling attractor bracketed from below (U11, drained sparse starts) and from above (U10b, late dense
              starts); same plot as fig_bracket() of paper/zenodo_whale/figures/make_figures.py, written here so that the
              network paper does not depend on the Zenodo preprint folder;
fig_bath.pdf: (a) emission efficiency of the U11 drained networks and their undrained U6 twins (seed means);
              (b) q_RP3(ell=5) for each study under three conventions for the momentum window (audit/u4_u11/a1_network.json).

    python3 paper/prd_network/make_figures.py"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "audit" / "u4_u11"))
import a_lib as L                                                     # noqa: E402  (the audit's independent estimator)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"                   # categorical slots 1-3 (validated)
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2})

fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw=dict(width_ratios=[1.15, 1]))
g = np.linspace(L.math.log(288) - 1.3, L.math.log(288) - 0.05, 40)
for mod, seeds, col in (("RP3", range(61, 67), BLUE), ("U1", range(61, 65), ORANGE)):
    for suf, ls, lab in (("_dr", "-", "drained (U11)"), ("", "--", "undrained (U6)")):
        E = np.array([L.on_grid(*L.emission(L.load(f"{mod}_N576_s{s}_xi0.5{suf}.json"))[::2], g) for s in seeds])
        a.plot(g, np.nanmean(E, 0), ls=ls, color=col, label=f"{'$\\mathbb{{RP}}^3$' if mod == 'RP3' else 'U(1)'}, {lab}")
a.set_xlabel(r"$\ell=\ln(m_r/H)$"); a.set_ylabel(r"emission efficiency $\varepsilon$")
a.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2); a.set_title("(a) the initial massive bath counted as emission", fontsize=9, loc="left", color=INK)

A = json.loads((ROOT / "audit" / "u4_u11" / "a1_network.json").read_text())["q_variants"]
studies = [("U11 q_RP3", "U11\ndrained (6)"), ("U6 q_RP3", "U6\n(6)"), ("U9 q_RP3 densas", "U9\ndense (2)"),
           ("U10a q_RP3 densas", "U10a\ndense (2)"), ("U10a q_RP3 poco densas", "U10a\nsparse (2)")]
conv = [("a,H del instante final del par (como el publicado)", "end of pair (published)", BLUE, "o"),
        ("a,H medios geométricos (base del auditor)", "geometric mean", ORANGE, "s"),
        ("a,H del instante inicial del par", "start of pair", AQUA, "^")]
for i, (key, lab) in enumerate(studies):
    for j, (ck, cl, col, mk) in enumerate(conv):
        v = A[key]["variants"][ck]
        b.errorbar(i + (j - 1) * 0.18, v["q0_at_ell5"], yerr=v["q0_std"], fmt=mk, ms=6, color=col, mec="white", mew=1, elinewidth=1.2,
                   capsize=0, label=cl if i == 0 else None)
b.axhline(1.0, color=MUTED, lw=1, ls=":"); b.text(4.45, 1.02, "$q=1$", color=MUTED, fontsize=7, ha="right")
b.set_xticks(range(len(studies))); b.set_xticklabels([s[1] for s in studies], fontsize=7)
b.set_ylabel(r"$q_{\mathbb{RP}^3}(\ell=5)$"); b.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, title="momentum window evaluated at", title_fontsize=7)
b.set_title("(b) spectral slope vs estimator convention", fontsize=9, loc="left", color=INK)
fig.tight_layout(); fig.savefig(HERE / "fig_bath.pdf", metadata={"CreationDate": None}); fig.savefig(HERE / "fig_bath.png", dpi=160)
print("wrote", HERE / "fig_bath.pdf")


# fig_bracket.pdf (U10b): identical to the Zenodo v1 figure (paper/zenodo_whale/figures/make_figures.py, fig_bracket)
RUNS = ROOT / "network" / "results" / "u4" / "runs"
C_RP3, C_U1 = "#b2182b", "#2166ac"


def _load(pattern):
    return [json.loads(p.read_text()) for p in sorted(RUNS.glob(pattern))]


def _series(run, key="xi"):
    return np.array([(r["ell"], r[key]) for r in run["rows"] if key in r]).T


with plt.rc_context({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "figure.dpi": 150, "savefig.bbox": "tight", "axes.spines.top": True, "axes.spines.right": True,
                     "lines.linewidth": 1.5, "axes.edgecolor": "black", "axes.labelcolor": "black",
                     "xtick.color": "black", "ytick.color": "black"}):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharex=True)
    for a, mod, title in ((ax[0], "RP3", r"(a) RP$^3$"), (ax[1], "U1", "(b) U(1)")):
        lo, hi = _load(f"{mod}_N576_s6*_xi0.5_dr.json"), _load(f"{mod}_N576_s7*_ab.json")
        for i, r in enumerate(lo):
            a.plot(*_series(r), color=C_U1, lw=0.7, alpha=0.7, label=r"from below ($\xi_i=0.5$, U11)" if i == 0 else None)
        for i, r in enumerate(hi):
            a.plot(*_series(r), color=C_RP3, lw=0.9, ls="-" if r["xi_init"] == 3 else "--",
                   label=(r"from above ($\xi_i=3$, $\ell_i=3.5$, U10b)" if r["xi_init"] == 3 else r"from above ($\xi_i=6$)")
                   if (i == 0 or r["xi_init"] == 6) else None)
        a.axvspan(4.36, 5.61, color="0.92", zorder=0)
        a.set_xlim(3.4, 5.8); a.set_ylim(0, 3.2 if mod == "RP3" else 2.4); a.set_xlabel(r"$\ell$"); a.set_ylabel(r"$\xi$")
        a.set_title(title, fontsize=9)
    ax[0].legend(loc="upper right", frameon=False, fontsize=6.5)
    fig.tight_layout(); fig.savefig(HERE / "fig_bracket.pdf", metadata={"CreationDate": None}); plt.close(fig)
print("wrote", HERE / "fig_bracket.pdf")
