#!/usr/bin/env python3
"""
Figure 5 of the BNI3 paper: BNI3 vs BoolNet on DREAM4 size10.

Writes each panel as a separate file (SVG with editable text + PNG at 600 dpi),
without an overall title or panel letter:
    Figure5_A_BalancedAccuracy  Figure5_B_F1  Figure5_C_Coverage  Figure5_D_PrecisionRecall

Style: white background with a light grid, #BD759C / #6CAAC4 boxes with a dark edge
and a thick black median, gray points without jitter, and no separate outlier markers
because each series is already a point.

Each point is one series (n = 25 per method and condition). The tests are Wilcoxon
signed-rank, paired by (Size, series, condition), two-sided, zero_method='wilcox',
with continuity correction.

Usage: python3 make_figure5.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon

# ==== PATHS (relative to reproducibility/DREAM4) ====
DREAM4_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BNI3_TSV = os.path.join(DREAM4_DIR, "results", "BNI3_size10_complete_evaluation_results.tsv")
BOOLNET_TSV = os.path.join(DREAM4_DIR, "results", "BoolNet_size10_complete_evaluation_results.tsv")
OUTPUT_DIR = os.path.join(DREAM4_DIR, "figures")

# ==== STYLE ====
COLORS = {"BoolNet": "#BD759C", "BNI3": "#6CAAC4"}
POINT_COLOR = "#808080"
BOX_EDGE = "#333333"
METHODS = ["BoolNet", "BNI3"]
CONDITIONS = ["21 Time Points", "5 Time Points"]
FIGSIZE = (6.5, 4.8)
# Fixed margins, the same for the 4 panels (no bbox_inches="tight", which crops each
# panel differently): same image size and same axes area.
MARGINS = dict(left=0.13, right=0.97, top=0.95, bottom=0.14)
FS_TICKS, FS_LABEL, FS_LEGEND, FS_STARS = 12, 14, 11, 18
BOX_WIDTH = 0.8        # width of each group's pair of boxes (same dodge as stripplot)
BOX_GAP = 0.35         # gap between the two boxes of a group

sns.set_style("whitegrid", {"grid.color": "#E6E6E6"})
plt.rcParams.update({
    "svg.fonttype": "none",          # text as text in the SVG
    "font.size": FS_TICKS,
    "xtick.labelsize": FS_TICKS,
    "ytick.labelsize": FS_TICKS,
    "axes.labelsize": FS_LABEL,
    "axes.labelweight": "bold",
    "legend.fontsize": FS_LEGEND,
})


# ==== DATA ====
def load(path, method):
    df = pd.read_csv(path, sep="\t")
    df["Rep"] = df["Rep"].astype(str)
    df["Condition"] = np.where(df["Rep"].str.endswith("_5times"), "5 Time Points", "21 Time Points")
    df["SeriesId"] = df["Rep"].str.replace("_5times$", "", regex=True)
    df["Tool"] = method
    return df


bni3 = load(BNI3_TSV, "BNI3")
boolnet = load(BOOLNET_TSV, "BoolNet")
data = pd.concat([boolnet, bni3], ignore_index=True)
paired = boolnet.merge(bni3, on=["Size", "SeriesId", "Condition"], suffixes=("_BoolNet", "_BNI3"))


def paired_p(metric, condition):
    g = paired[paired["Condition"] == condition]
    return wilcoxon(g[f"{metric}_BoolNet"], g[f"{metric}_BNI3"],
                    zero_method="wilcox", correction=True, alternative="two-sided").pvalue


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


# ==== CHECK TABLE ====
def print_table():
    print(f"{'Panel':<20}{'Condition':<16}{'n':>4}  {'BoolNet mean/median':>22}  "
          f"{'BNI3 mean/median':>20}  {'p Wilcoxon':>11}")
    for panel, metric in [("A Balanced Acc.", "Balanced_Accuracy"), ("B F1-Score", "F1_score"),
                          ("C Coverage", "Coverage_Ratio"), ("D Precision", "Precision"),
                          ("D Recall", "Recall")]:
        for cond in CONDITIONS:
            g = paired[paired["Condition"] == cond]
            bn, b3 = g[f"{metric}_BoolNet"], g[f"{metric}_BNI3"]
            print(f"{panel:<20}{cond:<16}{len(g):>4}  {bn.mean():>11.3f} / {bn.median():<8.3f}  "
                  f"{b3.mean():>9.3f} / {b3.median():<8.3f}  {paired_p(metric, cond):>11.4g}")


# ==== COMPONENTS ====
def save(fig, name):
    fig.subplots_adjust(**MARGINS)
    for ext, kw in (("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"), **kw)
    plt.close(fig)
    print(f"  -> {name}.svg / .png")


def boxplot_panel(metric, ylabel, ylim, name, hline=None, hline_color="#555555",
                  significance=False, legend=False):
    fig, ax = plt.subplots(figsize=FIGSIZE)

    if hline is not None:
        # zorder 1: above the grid line at the same value, below the boxes
        ax.axhline(hline, color=hline_color, linestyle="--", linewidth=1.5, zorder=1)

    sns.boxplot(data=data, x="Condition", y=metric, hue="Tool", order=CONDITIONS,
                hue_order=METHODS, palette=COLORS, saturation=1, width=BOX_WIDTH,
                gap=BOX_GAP, showfliers=False,  # outliers are already shown as points
                boxprops=dict(edgecolor=BOX_EDGE, linewidth=1.5),
                medianprops=dict(color="black", linewidth=2),
                whiskerprops=dict(color="black", linewidth=1.5),
                capprops=dict(color="black", linewidth=1.5), ax=ax)
    sns.stripplot(data=data, x="Condition", y=metric, hue="Tool", order=CONDITIONS,
                  hue_order=METHODS, dodge=True, jitter=False,
                  palette=[POINT_COLOR, POINT_COLOR], alpha=0.6, size=6,
                  legend=False, clip_on=False, zorder=3, ax=ax)

    if significance:
        half = BOX_WIDTH / 4  # center of each box relative to the group center
        span = ylim[1] - ylim[0]
        for ci, cond in enumerate(CONDITIONS):
            p = paired_p(metric, cond)
            y = data[data["Condition"] == cond][metric].max() + 0.05 * span
            h = 0.025 * span
            ax.plot([ci - half, ci - half, ci + half, ci + half], [y, y + h, y + h, y],
                    "k-", lw=1.5)
            ax.text(ci, y + h, stars(p), ha="center", va="bottom", fontsize=FS_STARS,
                    fontweight="bold")

    ax.set_xlabel("")
    ax.xaxis.grid(True)  # light vertical line at the center of each group
    ax.set_ylabel(ylabel, fontsize=FS_LABEL, fontweight="bold")
    ax.set_ylim(*ylim)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:len(METHODS)], labels[:len(METHODS)], loc="lower right",
                  frameon=True)
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

    save(fig, name)


def precision_recall_panel(name):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    # Limits slightly below 0 so markers at the origin are not clipped
    xlim, ylim = (-0.015, 0.6), (-0.015, 0.7)

    # Iso-F1 curves: R = F·P / (2P − F)
    for f in (0.1, 0.2, 0.3, 0.4):
        p = np.linspace(f / 2 + 1e-3, xlim[1], 400)  # R ≥ 0 only for P > F/2
        r = f * p / (2 * p - f)
        keep = (r >= 0) & (r <= ylim[1])
        ax.plot(p[keep], r[keep], linestyle=":", color="#AAAAAA", linewidth=1.2, zorder=1)
        ax.text(xlim[1] - 0.005, r[keep][-1] + 0.012, f"F1 = {f}", ha="right", va="bottom",
                fontsize=9, color="#888888",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1))

    markers = {"21 Time Points": "o", "5 Time Points": "^"}
    short = {"21 Time Points": "21 TP", "5 Time Points": "5 TP"}
    for method in METHODS:
        for cond in CONDITIONS:
            g = data[(data["Tool"] == method) & (data["Condition"] == cond)]
            # clip_on=False: some series have precision and recall 0, right in the corner
            ax.scatter(g["Precision"], g["Recall"], s=50, marker=markers[cond],
                       color=COLORS[method], edgecolor="#333333", linewidth=0.6,
                       alpha=0.85, zorder=3, clip_on=False, label=f"{method} {short[cond]}")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.arange(0, 0.61, 0.1))
    ax.set_yticks(np.arange(0, 0.71, 0.1))
    ax.set_xlabel("Precision", fontsize=FS_LABEL, fontweight="bold")
    ax.set_ylabel("Recall", fontsize=FS_LABEL, fontweight="bold")
    # Upper right: the only area without data (above the F1 = 0.4 curve);
    # in the upper left the legend would cover BoolNet 21 TP series.
    ax.legend(loc="upper right", frameon=True)
    save(fig, name)


# ==== MAIN ====
if __name__ == "__main__":
    print(f"BNI3:    {len(bni3)} rows | BoolNet: {len(boolnet)} rows | pairs: {len(paired)}\n")
    print_table()
    origin = data[(data["Precision"] == 0) & (data["Recall"] == 0)]
    print("\nPanel D - points at (0, 0) per group:")
    print(origin.groupby(["Tool", "Condition"]).size().to_string())
    print("\nWriting panels to", OUTPUT_DIR)
    boxplot_panel("Balanced_Accuracy", "Balanced Accuracy", (0.0, 0.8),
                  "Figure5_A_BalancedAccuracy", hline=0.5, legend=True)
    boxplot_panel("F1_score", "F1-Score", (0.0, 0.5), "Figure5_B_F1")
    # 3.5 instead of 3.2: with 3.2 the 21 Time Points asterisks fall outside the axis
    boxplot_panel("Coverage_Ratio", "Coverage Ratio", (0.0, 3.5), "Figure5_C_Coverage",
                  hline=1.0, hline_color="#7B6FD0", significance=True)
    precision_recall_panel("Figure5_D_PrecisionRecall")
