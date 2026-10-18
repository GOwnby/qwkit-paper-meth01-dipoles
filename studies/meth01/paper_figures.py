"""Figures for the METH-01 note from scan.json -> ../../paper_meth01/figures/*.pdf"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
outdir = os.path.join(here, "..", "..", "paper_meth01", "figures"); os.makedirs(outdir, exist_ok=True)
BLUE, ORANGE, AQUA, PURPLE, RED = "#2a78d6", "#eb6834", "#1baf7a", "#8a4fbf", "#c4314b"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "legend.fontsize": 6.6, "axes.linewidth": 0.6, "axes.edgecolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 1.3, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02, "pdf.fonttype": 42})
R = json.load(open(os.path.join(here, "scan.json")))
M = {"raw": (INK2, "no orthogonalization"), "norm": ("#9a9a9a", "two-band norm only"), "GS": (PURPLE, "Gram–Schmidt"),
     "L-bound": (RED, "Löwdin, bound states"), "L-b2": (ORANGE, "Löwdin + 2 states"), "L-b10": (BLUE, "Löwdin + 10 states"),
     "L-E0.5": (AQUA, "Löwdin + 0.5 eV"), "L-all": ("#0b6e4f", "Löwdin, all 120 states")}


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)


def worst(r, k):
    return 100 * max(abs(v) for v in r["err"][k].values())


fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4), sharey=True, gridspec_kw=dict(wspace=0.12))
for ax, name, title in zip(axs, ("ADQW", "TW", "SQW"), ("double well 64/10/29 Å", "triple well, 55 kV/cm", "single well 60 Å")):
    rr = sorted([r for r in R if r["name"] == name and not r["matched"]], key=lambda r: r["pad"])
    for k, (c, lab) in M.items():
        ax.plot([r["pad"] for r in rr], [max(worst(r, k), 0.05) for r in rr], color=c, marker="o", ms=2.4, label=lab,
                ls="--" if k == "norm" else "-")
    ax.set_yscale("log"); ax.set_xscale("log"); ax.set_xlabel("barrier padding (Å)")
    ax.set_xticks([100, 200, 400, 800]); ax.set_xticklabels(["100", "200", "400", "800"]); ax.minorticks_off(); ax.set_title(title, fontsize=7.8); clean(ax)
axs[0].set_ylabel("largest dipole error vs Kane (%)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.24), fontsize=6.4)
fig.savefig(os.path.join(outdir, "fig1_errors.pdf")); plt.close(fig)
print("written", outdir)
