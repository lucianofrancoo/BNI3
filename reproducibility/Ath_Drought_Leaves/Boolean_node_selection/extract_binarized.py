#!/usr/bin/env python3
"""
Extracts only the 15 nodes of selected_nodes.tsv from the binarized (SSD) matrices.

Input:    {Ath,Sly}_Drought_Leaves/Binarization/
          ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD.tsv
          (one gene per column, one time point per row, no row names)
Output:   same folder, suffix _selected_nodes.tsv, same format; columns in the order of
          selected_nodes.tsv, named "ID (symbol)". In tomato the symbol is that of the
          Arabidopsis gene it is an ortholog of. Each column is checked to reproduce the
          pattern in selected_nodes.tsv.

It also extracts the mean normalized count (vst) of those genes, which is what was
binarized:
Input:    {Ath,Sly}_Drought_Leaves/Selected_genes/ConteoNormalizado_mean_{Control,Drought}_aba_genes.tsv
Output:   same folder, suffix _selected_nodes.tsv; rows in the order of
          selected_nodes.tsv, ID column = "ID (symbol)".

Each output also gets a *_symbol.tsv version with only the symbol as name (in tomato, the
symbol of the Arabidopsis ortholog), so the nodes have the same names in both species.
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
        # check: the column read = control / drought half of the 8-bit pattern
        for g, p in zip(genes, nodes[pat]):
            assert "".join(sub[g]) == p.split("|")[i], (sp, cond, g)
        out = f.with_name(f.stem + "_selected_nodes.tsv")
        sub.set_axis([f"{g} ({s})" for g, s in zip(genes, nodes.arabidopsis_gene)],
                     axis=1).to_csv(out, sep="\t", index=False)
        sub.set_axis(nodes.arabidopsis_gene.tolist(), axis=1).to_csv(
            out.with_name(out.stem + "_symbol.tsv"), sep="\t", index=False)
        print(f"{out.relative_to(ROOT)}: {sub.shape[1]} genes x {sub.shape[0]} time points")

    # mean normalized count of the same genes
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
        print(f"{out.relative_to(ROOT)}: {len(sub)} genes x {sub.shape[1]} time points")
