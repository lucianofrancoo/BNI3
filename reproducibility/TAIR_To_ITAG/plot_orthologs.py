#!/usr/bin/env python3
"""
Figures of the Arabidopsis -> tomato ortholog comparison (ABA genes).

Input:    Ath_ABA_orthologs_combined.tsv (from combine_orthologs.py)
Methods:  OrthoFinder and Ensembl Compara (both based on gene trees).
Universe: pairs (Ath gene, Sly gene) proposed by at least one of the two methods.

Outputs (vector PDF + PNG):
  fig_orthologs_summary   A: Ath genes with an ortholog according to each method (Venn)
                          B: Ath genes with at least one ortholog, by group: union (OR) and each method
  fig_orthologs_detail    pair-by-pair matrix (Ath gene -> Solyc) x method, by group
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle

HERE = Path(__file__).resolve().parent
GENES = HERE.parent / "NCBI_To_ITAG" / "Ath" / "Ath_ABA_genes.tsv"

METHODS = ["OrthoFinder", "Compara"]
COLOR = {"OrthoFinder": "#2a78d6", "Compara": "#eb6834"}
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
EMPTY = "#dcdbd6"
NONE = "No ortholog"
GROUP_LABEL = {"TF": "TF", "metabolic": "Metabolic", "catabolic": "Catabolic",
               "transport": "Transport", "receptor": "Receptor",
               "phosphatase_cladeA": "Phosphatase clade A",
               "phosphatase_cladeB": "Phosphatase clade B", "kinase": "Kinase",
               "ligase": "Ligase", "targets": "Targets"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "pdf.fonttype": 42,
})

genes = pd.read_csv(GENES, sep="\t")
group_order = genes.group.drop_duplicates().tolist()
m = pd.read_csv(HERE / "Ath_ABA_orthologs_combined.tsv", sep="\t")
p = m[m[METHODS].any(axis=1)].copy()


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- summary
fig, (axv, axg) = plt.subplots(1, 2, figsize=(11, 4.6),
                               gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.45})

# A: Venn of Ath genes with ≥1 ortholog according to each method (equal circles; the
#    areas are not proportional). The union is the total of genes with an ortholog.
g_of, g_c = set(p[p.OrthoFinder].ath), set(p[p.Compara].ath)
only_of, only_c, both = len(g_of - g_c), len(g_c - g_of), len(g_of & g_c)
r, d = 1.0, 0.62
for cx, meth in [(-d, "OrthoFinder"), (d, "Compara")]:
    axv.add_patch(Circle((cx, 0), r, facecolor=COLOR[meth], alpha=0.22,
                         edgecolor=COLOR[meth], linewidth=2))
for x, n in [(-d - 0.42, only_of), (0, both), (d + 0.42, only_c)]:
    axv.text(x, 0, str(n), ha="center", va="center", fontsize=14, color=INK)
axv.text(-d - 0.2, r + 0.15, f"OrthoFinder ({only_of + both})", ha="right",
         va="bottom", fontsize=10, color=INK)
axv.text(d + 0.2, r + 0.15, f"Compara ({only_c + both})", ha="left", va="bottom",
         fontsize=10, color=INK)
axv.text(0, -r - 0.22, f"{only_of + both + only_c} of {len(genes)} Arabidopsis genes\n"
         "with ≥1 tomato ortholog", ha="center", va="top", fontsize=9, color=INK2,
         linespacing=1.4)
axv.set_xlim(-2.3, 2.3)
axv.set_ylim(-1.6, 1.6)
axv.set_aspect("equal")
axv.axis("off")

# B: Ath genes with ≥1 ortholog per method, by functional group
cov = []
for g in group_order:
    ids = genes[genes.group == g].gene_id
    row = {"group": g, "n": len(ids)}
    for meth in METHODS:
        row[meth] = ids.isin(p[p[meth]].ath).sum()
    row["any"] = ids.isin(p.ath).sum()
    cov.append(row)
cov = pd.DataFrame(cov)
UNION = "OrthoFinder or Compara"
bars = [(UNION, "any", INK2)] + [(mm, mm, COLOR[mm]) for mm in METHODS]
bh = 0.24
for k, (label, col, color) in enumerate(bars):
    ys = [i + (k - 1) * (bh + 0.03) for i in range(len(cov))]
    axg.barh(ys, cov[col], height=bh, color=color, label=label)
    if col == "any":
        for y_, v, n in zip(ys, cov["any"], cov.n):
            axg.text(v + 0.25, y_, f"{v}/{n}", va="center", ha="left", fontsize=7.5,
                     color=INK)
axg.set_yticks(range(len(cov)))
axg.set_yticklabels([f"{GROUP_LABEL[g]} ({n})" for g, n in zip(cov.group, cov.n)])
axg.invert_yaxis()
axg.set_xlim(0, cov.n.max() + 2.5)
axg.set_xlabel("Arabidopsis genes with ≥1 tomato ortholog")
axg.xaxis.grid(True, color=GRID, lw=0.6)
axg.set_axisbelow(True)
for s_ in ("top", "right", "left"):
    axg.spines[s_].set_visible(False)
axg.tick_params(length=0)
axg.legend(frameon=False, loc="lower right", fontsize=8)
# panel letters at the same height (the Venn has a 1:1 aspect and sits lower)
fig.canvas.draw()
top = axg.get_position().y1 + 0.04
for ax, letter in [(axv, "A"), (axg, "B")]:
    fig.text(ax.get_position().x0 - 0.02, top, letter, fontsize=12, fontweight="bold",
             color=INK, ha="left", va="bottom")
save(fig, "fig_orthologs_summary")

# ---------------------------------------------------------------- detail
rows = []
for g in group_order:
    for gg in genes[genes.group == g].itertuples():
        d_ = p[p.ath == gg.gene_id].copy()
        d_["score"] = d_[METHODS].sum(axis=1)
        d_ = d_.sort_values(["score", "OrthoFinder", "Compara"], ascending=False)
        if d_.empty:
            rows.append(dict(group=g, Symbol=gg.Symbol, sly=NONE, **{k: False for k in METHODS}))
        for r_ in d_.itertuples():
            rows.append(dict(group=g, Symbol=gg.Symbol, sly=r_.sly,
                             **{k: getattr(r_, k) for k in METHODS}))
rows = pd.DataFrame(rows)

# split into 3 columns without breaking functional groups
ncol = 3
target = len(rows) / ncol
cols, cur, acc = [[] for _ in range(ncol)], 0, 0
for g in group_order:
    n = (rows.group == g).sum()
    if acc and acc + n / 2 > target * (cur + 1) and cur < ncol - 1:
        cur += 1
    cols[cur].append(g)
    acc += n
max_rows = max(sum((rows.group == g).sum() + 1.5 for g in c) for c in cols)

X0, DX = 2.1, 0.45  # position of the first column of dots and spacing
fig, axes = plt.subplots(1, ncol, figsize=(14, max_rows * 0.165 + 1.2))
for ax, gl in zip(axes, cols):
    y = 0
    for g in gl:
        ax.text(-0.35, y, GROUP_LABEL[g], fontsize=9, fontweight="bold", color=INK,
                va="center")
        y += 1
        prev_sym = None
        for r_ in rows[rows.group == g].itertuples():
            if r_.Symbol != prev_sym:
                ax.text(-0.35, y, r_.Symbol, fontsize=7.5, color=INK, va="center",
                        style="italic")
                ax.axhline(y - 0.5, color=GRID, lw=0.5)
                prev_sym = r_.Symbol
            none = r_.sly == NONE
            ax.text(0.75, y, r_.sly, fontsize=7, color=MUTED if none else INK2,
                    va="center", ha="left", style="italic" if none else "normal")
            for j, meth in enumerate(METHODS):
                on = getattr(r_, meth)
                ax.scatter(X0 + j * DX, y, s=30 if on else 12,
                           color=COLOR[meth] if on else EMPTY, zorder=3)
            y += 1
        y += 0.5
    for j, meth in enumerate(METHODS):
        ax.text(X0 + j * DX, -1.3, meth, rotation=45, ha="left", va="bottom",
                fontsize=8, color=INK)
    ax.set_xlim(-0.4, X0 + (len(METHODS) - 1) * DX + 0.35)
    ax.set_ylim(max_rows, -2.5)
    ax.axis("off")
fig.subplots_adjust(top=0.99, bottom=0.01, left=0.02, right=0.99, wspace=0.08)
save(fig, "fig_orthologs_detail")
print(f"Ath genes — OF only: {only_of}  both: {both}  Compara only: {only_c}  pairs: {len(p)}")
print(cov.to_string(index=False))
