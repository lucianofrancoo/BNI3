#!/usr/bin/env python3
"""
Figure of the nodes selected for the Boolean networks (selected_nodes.tsv).

Per node (Ath gene + tomato ortholog):
  - functional group
  - SSD control|drought binary pattern in Arabidopsis and in tomato (1 = filled cell)
  - maximum log2FC (WD vs CTRL) of each species: the time point with the largest |log2FC|
    among the significant ones (padj < 0.05); if none is significant, the largest overall
    |log2FC|. Filled bar = significant, outline = not significant. Axis truncated at
    ±LFC_LIM; out-of-range values are annotated.
Output: fig_boolean_nodes.pdf / .png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LFC_LIM = 6

INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
ON, OFF = "#3d3d3a", "#e9e8e4"
UP, DOWN = "#D73027", "#4575B4"   # same colors as the DGE dot plots
GROUP_ORDER = ["TF", "Metabolic", "Catabolic", "Transport", "Receptor",
               "Phosphatase clade A", "Phosphatase clade B", "Kinase", "Ligase", "Targets"]
TIMES = {"Ath": ["0", "7200", "15840", "20160"], "Sly": ["d4", "d5", "d8", "d10"]}
SPECIES = {"Ath": "Arabidopsis", "Sly": "Tomato"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "pdf.fonttype": 42,
})

nodes = pd.read_csv(HERE / "selected_nodes.tsv", sep="\t",
                    dtype={"ath_pattern": str, "sly_pattern": str})
nodes["_g"] = nodes.functional_group.map({g: i for i, g in enumerate(GROUP_ORDER)})
nodes = nodes.sort_values("_g", kind="stable").reset_index(drop=True)


def max_lfc(table, gene):
    t = table[table.gene_id == gene]
    sig = t[t.padj < 0.05]
    pool = sig if len(sig) else t
    r = pool.loc[pool.log2FoldChange.abs().idxmax()]
    return r.log2FoldChange, len(sig) > 0, r.comparison


dge = {
    "Ath": pd.read_csv(ROOT / "Ath_Drought_Leaves/Selected_genes/dotplot_aba_genes_table.tsv", sep="\t"),
    "Sly": pd.read_csv(ROOT / "Sly_Drought_Leaves/Selected_genes/dotplot_aba_genes_table.tsv", sep="\t"),
}
for sp, col in [("Ath", "arabidopsis_locus"), ("Sly", "tomato_gene")]:
    vals = [max_lfc(dge[sp], g) for g in nodes[col]]
    nodes[f"{sp}_lfc"] = [v[0] for v in vals]
    nodes[f"{sp}_sig"] = [v[1] for v in vals]
    nodes[f"{sp}_time"] = [v[2] for v in vals]

# vertical positions: an extra gap between groups
ys, y, prev = [], 0.0, None
for g in nodes.functional_group:
    if prev is not None and g != prev:
        y += 0.6
    ys.append(y)
    y += 1
    prev = g
nodes["y"] = ys
ymax = y - 1

fig = plt.figure(figsize=(12, 0.34 * (ymax + 1) + 1.6))
gs = fig.add_gridspec(1, 5, width_ratios=[2.0, 1.5, 1.5, 1.25, 1.25], wspace=0.12)
ax_lab, ax_pa, ax_ps, ax_la, ax_ls = [fig.add_subplot(gs[0, i]) for i in range(5)]
for ax in (ax_lab, ax_pa, ax_ps, ax_la, ax_ls):
    ax.set_ylim(ymax + 0.7, -2.9)

# labels: group, Ath gene (italics), tomato gene
ax_lab.set_xlim(0, 1)
ax_lab.axis("off")
for g, d in nodes.groupby("functional_group", sort=False):
    y0, y1 = d.y.min(), d.y.max()
    ax_lab.plot([0.36, 0.36], [y0 - 0.35, y1 + 0.35], color=INK2, lw=1)
    ax_lab.text(0.33, (y0 + y1) / 2, g, ha="right", va="center", fontsize=8.5,
                fontweight="bold", color=INK)
for r in nodes.itertuples():
    ax_lab.text(0.40, r.y, r.arabidopsis_gene, ha="left", va="center", fontsize=9,
                style="italic", color=INK)
    ax_lab.text(0.66, r.y, r.tomato_gene, ha="left", va="center", fontsize=8, color=INK2)
ax_lab.text(0.40, -1.2, "Arabidopsis", ha="left", va="bottom", fontsize=8.5, color=INK2)
ax_lab.text(0.66, -1.2, "Tomato ortholog", ha="left", va="bottom", fontsize=8.5, color=INK2)


def draw_patterns(ax, sp, col):
    # 8 cells: 4 control + gap + 4 drought
    xs = [0, 1, 2, 3, 4.5, 5.5, 6.5, 7.5]
    for r in nodes.itertuples():
        bits = getattr(r, col).replace("|", "")
        for x, b in zip(xs, bits):
            ax.add_patch(Rectangle((x - 0.42, r.y - 0.36), 0.84, 0.72,
                                   facecolor=ON if b == "1" else OFF, edgecolor="none"))
    for x, t in zip(xs, TIMES[sp] * 2):
        ax.text(x, -0.75, t, ha="center", va="bottom", fontsize=6.5, color=INK2,
                rotation=90)
    ax.text(1.5, -2.0, "Control", ha="center", va="bottom", fontsize=8, color=INK2)
    ax.text(6.0, -2.0, "Drought", ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_xlim(-0.7, 8.2)
    ax.axis("off")
    ax.set_title(f"{SPECIES[sp]} binary pattern", fontsize=9.5, fontweight="bold",
                 color=INK, pad=4)


draw_patterns(ax_pa, "Ath", "ath_pattern")
draw_patterns(ax_ps, "Sly", "sly_pattern")


def draw_lfc(ax, sp):
    for r in nodes.itertuples():
        v, sig = getattr(r, f"{sp}_lfc"), getattr(r, f"{sp}_sig")
        vc = max(min(v, LFC_LIM), -LFC_LIM)
        color = UP if v > 0 else DOWN
        ax.barh(r.y, vc, height=0.62, color=color if sig else SURFACE,
                edgecolor=color, linewidth=1.2)
        if abs(v) > LFC_LIM:
            ax.text(vc + (0.2 if v > 0 else -0.2), r.y, f"{v:.1f}", va="center",
                    ha="left" if v > 0 else "right", fontsize=7, color=INK)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlim(-LFC_LIM - 1.5, LFC_LIM + 1.5)
    ax.set_xticks(range(-LFC_LIM, LFC_LIM + 1, 3))
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0, labelsize=7.5)
    ax.set_xlabel("max log2FC (WD vs CTRL)", fontsize=8)
    ax.set_title(SPECIES[sp], fontsize=9.5, fontweight="bold", color=INK, pad=4)


draw_lfc(ax_la, "Ath")
draw_lfc(ax_ls, "Sly")

# legend
leg = [Rectangle((0, 0), 1, 1, facecolor=ON), Rectangle((0, 0), 1, 1, facecolor=OFF),
       Rectangle((0, 0), 1, 1, facecolor=UP), Rectangle((0, 0), 1, 1, facecolor=DOWN),
       Rectangle((0, 0), 1, 1, facecolor=SURFACE, edgecolor=INK2)]
fig.legend(leg, ["Active (1)", "Inactive (0)", "Up", "Down", "Not significant (FDR ≥ 0.05)"],
           loc="lower center", ncol=5, frameon=False, fontsize=8,
           bbox_to_anchor=(0.55, -0.005))
fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.12)
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"fig_boolean_nodes.{ext}", dpi=300, bbox_inches="tight")

nodes[["functional_group", "arabidopsis_gene", "tomato_gene", "Ath_lfc", "Ath_sig",
       "Ath_time", "Sly_lfc", "Sly_sig", "Sly_time"]].to_csv(
    HERE / "fig_boolean_nodes_lfc.tsv", sep="\t", index=False, float_format="%.2f")
print(nodes[["arabidopsis_gene", "Ath_lfc", "Ath_sig", "Sly_lfc", "Sly_sig"]].round(2).to_string(index=False))
