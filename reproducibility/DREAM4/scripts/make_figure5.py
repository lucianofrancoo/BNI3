#!/usr/bin/env python3
"""
Figura 5 del paper BNI3: BNI3 vs BoolNet en DREAM4 size10 (resultados nuevos).

Genera cada panel como archivo separado (SVG con texto editable + PNG 600 dpi),
sin título general ni letra de panel:
    Figure5_A_BalancedAccuracy  Figure5_B_F1  Figure5_C_Coverage  Figure5_D_PrecisionRecall

Estética: la de la Figura 5 anterior (fondo blanco con grilla suave, cajas
#BD759C / #6CAAC4 con borde oscuro y mediana negra gruesa, puntos grises sin
jitter, sin marcar outliers aparte porque cada serie ya es un punto).

Cada punto es una serie (n = 25 por método y condición). Los tests son Wilcoxon
signed-rank pareados por (Size, serie, condición), two-sided, zero_method='wilcox',
con corrección de continuidad.

Uso: python3 make_figure5.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon

# ==== RUTAS (relativas a reproducibility/DREAM4) ====
DREAM4_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BNI3_TSV = os.path.join(DREAM4_DIR, "results", "BNI3_size10_complete_evaluation_results.tsv")
BOOLNET_TSV = os.path.join(DREAM4_DIR, "results", "BoolNet_size10_complete_evaluation_results.tsv")
OUTPUT_DIR = os.path.join(DREAM4_DIR, "figures")

# ==== ESTILO (igual a BNI3_vs_Boolnet_graph.py) ====
COLORS = {"BoolNet": "#BD759C", "BNI3": "#6CAAC4"}
POINT_COLOR = "#808080"
BOX_EDGE = "#333333"
METHODS = ["BoolNet", "BNI3"]
CONDITIONS = ["21 Time Points", "5 Time Points"]
FIGSIZE = (6.5, 4.8)
# Márgenes fijos e iguales en los 4 paneles (sin bbox_inches="tight", que recorta
# cada panel distinto): mismo tamaño de imagen y misma área de ejes.
MARGINS = dict(left=0.13, right=0.97, top=0.95, bottom=0.14)
FS_TICKS, FS_LABEL, FS_LEGEND, FS_STARS = 12, 14, 11, 18
BOX_WIDTH = 0.8        # ancho del par de cajas de cada grupo (mismo dodge que stripplot)
BOX_GAP = 0.35         # separación entre las dos cajas de un grupo

sns.set_style("whitegrid", {"grid.color": "#E6E6E6"})
plt.rcParams.update({
    "svg.fonttype": "none",          # texto como texto en el SVG
    "font.size": FS_TICKS,
    "xtick.labelsize": FS_TICKS,
    "ytick.labelsize": FS_TICKS,
    "axes.labelsize": FS_LABEL,
    "axes.labelweight": "bold",
    "legend.fontsize": FS_LEGEND,
})


# ==== DATOS ====
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


# ==== TABLA DE VERIFICACIÓN ====
def print_table():
    print(f"{'Panel':<20}{'Condición':<16}{'n':>4}  {'BoolNet media/mediana':>22}  "
          f"{'BNI3 media/mediana':>20}  {'p Wilcoxon':>11}")
    for panel, metric in [("A Balanced Acc.", "Balanced_Accuracy"), ("B F1-Score", "F1_score"),
                          ("C Coverage", "Coverage_Ratio"), ("D Precision", "Precision"),
                          ("D Recall", "Recall")]:
        for cond in CONDITIONS:
            g = paired[paired["Condition"] == cond]
            bn, b3 = g[f"{metric}_BoolNet"], g[f"{metric}_BNI3"]
            print(f"{panel:<20}{cond:<16}{len(g):>4}  {bn.mean():>11.3f} / {bn.median():<8.3f}  "
                  f"{b3.mean():>9.3f} / {b3.median():<8.3f}  {paired_p(metric, cond):>11.4g}")


# ==== COMPONENTES ====
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
        # zorder 1: sobre la línea de grilla del mismo valor, bajo las cajas
        ax.axhline(hline, color=hline_color, linestyle="--", linewidth=1.5, zorder=1)

    sns.boxplot(data=data, x="Condition", y=metric, hue="Tool", order=CONDITIONS,
                hue_order=METHODS, palette=COLORS, saturation=1, width=BOX_WIDTH,
                gap=BOX_GAP, showfliers=False,  # outliers ya están como puntos
                boxprops=dict(edgecolor=BOX_EDGE, linewidth=1.5),
                medianprops=dict(color="black", linewidth=2),
                whiskerprops=dict(color="black", linewidth=1.5),
                capprops=dict(color="black", linewidth=1.5), ax=ax)
    sns.stripplot(data=data, x="Condition", y=metric, hue="Tool", order=CONDITIONS,
                  hue_order=METHODS, dodge=True, jitter=False,
                  palette=[POINT_COLOR, POINT_COLOR], alpha=0.6, size=6,
                  legend=False, clip_on=False, zorder=3, ax=ax)

    if significance:
        half = BOX_WIDTH / 4  # centro de cada caja respecto del centro del grupo
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
    ax.xaxis.grid(True)  # línea vertical suave en el centro de cada grupo, como la referencia
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
    # Límites algo bajo 0 para que los marcadores en el origen no queden cortados
    xlim, ylim = (-0.015, 0.6), (-0.015, 0.7)

    # Curvas iso-F1: R = F·P / (2P − F)
    for f in (0.1, 0.2, 0.3, 0.4):
        p = np.linspace(f / 2 + 1e-3, xlim[1], 400)  # R ≥ 0 solo para P > F/2
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
            # clip_on=False: hay series con precisión y recall 0 justo en la esquina
            ax.scatter(g["Precision"], g["Recall"], s=50, marker=markers[cond],
                       color=COLORS[method], edgecolor="#333333", linewidth=0.6,
                       alpha=0.85, zorder=3, clip_on=False, label=f"{method} {short[cond]}")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks(np.arange(0, 0.61, 0.1))
    ax.set_yticks(np.arange(0, 0.71, 0.1))
    ax.set_xlabel("Precision", fontsize=FS_LABEL, fontweight="bold")
    ax.set_ylabel("Recall", fontsize=FS_LABEL, fontweight="bold")
    # Arriba a la derecha: es la única zona sin datos (sobre la curva F1 = 0.4);
    # arriba a la izquierda la leyenda tapa series de BoolNet 21 TP.
    ax.legend(loc="upper right", frameon=True)
    save(fig, name)


# ==== MAIN ====
if __name__ == "__main__":
    print(f"BNI3:    {len(bni3)} filas | BoolNet: {len(boolnet)} filas | pares: {len(paired)}\n")
    print_table()
    origin = data[(data["Precision"] == 0) & (data["Recall"] == 0)]
    print("\nPanel D - puntos en (0, 0) por grupo:")
    print(origin.groupby(["Tool", "Condition"]).size().to_string())
    print("\nGenerando paneles en", OUTPUT_DIR)
    boxplot_panel("Balanced_Accuracy", "Balanced Accuracy", (0.0, 0.8),
                  "Figure5_A_BalancedAccuracy", hline=0.5, legend=True)
    boxplot_panel("F1_score", "F1-Score", (0.0, 0.5), "Figure5_B_F1")
    # 3.5 en vez de 3.2: con 3.2 los asteriscos de 21 Time Points quedan fuera del eje
    boxplot_panel("Coverage_Ratio", "Coverage Ratio", (0.0, 3.5), "Figure5_C_Coverage",
                  hline=1.0, hline_color="#7B6FD0", significance=True)
    precision_recall_panel("Figure5_D_PrecisionRecall")
