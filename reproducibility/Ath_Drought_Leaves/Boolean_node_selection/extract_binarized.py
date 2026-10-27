#!/usr/bin/env python3
"""
Extrae de las matrices binarizadas (SSD) solo los 15 nodos de selected_nodes.tsv.

Entrada:  {Ath,Sly}_Drought_Leaves/Binarization/
          ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD.tsv
          (un gen por columna, un tiempo por fila, sin nombres de fila)
Salida:   mismo directorio, sufijo _selected_nodes.tsv, mismo formato; columnas en el
          orden de selected_nodes.tsv, nombradas "ID (símbolo)". En tomate el símbolo es
          el del gen de Arabidopsis del que es ortólogo. Se verifica que cada columna
          reproduzca el patrón de selected_nodes.tsv.

También extrae el conteo normalizado promedio (vst) de esos genes, que es lo que se
binarizó:
Entrada:  {Ath,Sly}_Drought_Leaves/Selected_genes/ConteoNormalizado_mean_{Control,Drought}_aba_genes.tsv
Salida:   mismo directorio, sufijo _selected_nodes.tsv; filas en el orden de
          selected_nodes.tsv, columna ID = "ID (símbolo)".

De cada salida se escribe además una versión *_symbol.tsv con solo el símbolo como nombre
(en tomate, el símbolo del ortólogo de Arabidopsis), para que los nodos se llamen igual en
las dos especies.
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
METHOD = "SSD"

nodes = pd.read_csv(HERE / "selected_nodes.tsv", sep="\t",
                    dtype={"ath_pattern": str, "sly_pattern": str})

for sp, col, pat in [("Ath", "arabidopsis_locus", "ath_pattern"),
                     ("Sly", "tomato_gene", "sly_pattern")]:
    genes = nodes[col].tolist()
    for i, cond in enumerate(("Control", "Drought")):
        f = (ROOT / f"{sp}_Drought_Leaves" / "Binarization" /
             f"ConteoNormalizado_mean_{cond}_aba_genes_binarized_{METHOD}.tsv")
        m = pd.read_csv(f, sep="\t", dtype=str)
        sub = m[genes]
        # verificación: la columna leída = mitad control / sequía del patrón de 8 bits
        for g, p in zip(genes, nodes[pat]):
            assert "".join(sub[g]) == p.split("|")[i], (sp, cond, g)
        out = f.with_name(f.stem + "_selected_nodes.tsv")
        sub.set_axis([f"{g} ({s})" for g, s in zip(genes, nodes.arabidopsis_gene)],
                     axis=1).to_csv(out, sep="\t", index=False)
        sub.set_axis(nodes.arabidopsis_gene.tolist(), axis=1).to_csv(
            out.with_name(out.stem + "_symbol.tsv"), sep="\t", index=False)
        print(f"{out.relative_to(ROOT)}: {sub.shape[1]} genes x {sub.shape[0]} tiempos")

    # conteo normalizado promedio de los mismos genes
    labels = {g: f"{g} ({s})" for g, s in zip(genes, nodes.arabidopsis_gene)}
    for cond in ("Control", "Drought"):
        f = (ROOT / f"{sp}_Drought_Leaves" / "Selected_genes" /
             f"ConteoNormalizado_mean_{cond}_aba_genes.tsv")
        m = pd.read_csv(f, sep="\t", dtype={"ID": str}).set_index("ID")
        sub = m.loc[genes]
        out = f.with_name(f.stem + "_selected_nodes.tsv")
        sub.rename(index=labels).reset_index().to_csv(out, sep="\t", index=False)
        sub.set_axis(nodes.arabidopsis_gene.tolist()).rename_axis("ID").reset_index().to_csv(
            out.with_name(out.stem + "_symbol.tsv"), sep="\t", index=False)
        print(f"{out.relative_to(ROOT)}: {len(sub)} genes x {sub.shape[1]} tiempos")
