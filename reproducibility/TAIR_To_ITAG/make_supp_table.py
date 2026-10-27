#!/usr/bin/env python3
"""
Supplementary table (TSV) with the exact data behind fig_orthologs_summary / _detail.

Input:  Ath_ABA_orthologs_final.tsv (from combine_orthologs.py) and the Ensembl Plants 62
        GFF (ITAG descriptions of the tomato genes).
Output: Table_S_ABA_orthologs.tsv — one row per pair (Ath gene, tomato gene) and one row
        per Ath gene without an ortholog. A pair is an ortholog if OrthoFinder or
        Ensembl Compara infers it.
"""
import gzip
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

HERE = Path(__file__).resolve().parent
GENES = HERE.parent / "NCBI_To_ITAG" / "Ath" / "Ath_ABA_genes.tsv"
ENS_GFF = HERE.parent / "NCBI_To_ITAG" / "Ath" / "Solanum_lycopersicum.SL3.0.62.gff3.gz"
OUT = HERE / "Table_S_ABA_orthologs.tsv"

GROUP_LABEL = {"TF": "TF", "metabolic": "Metabolic", "catabolic": "Catabolic",
               "transport": "Transport", "receptor": "Receptor",
               "phosphatase_cladeA": "Phosphatase clade A",
               "phosphatase_cladeB": "Phosphatase clade B", "kinase": "Kinase",
               "ligase": "Ligase", "targets": "Targets"}
REL = {"one2one": "one-to-one", "one2many": "one-to-many", "many2one": "many-to-one",
       "many2many": "many-to-many"}


def itag_descriptions(path):
    desc = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                gid = re.search(r"gene_id=([^;\n]+)", f[8])
                d = re.search(r"description=([^;\n]+)", f[8])
                if gid and d:
                    desc[gid.group(1)] = re.sub(r"\s*\[Source:.*\]$", "", unquote(d.group(1)))
    return desc


genes = pd.read_csv(GENES, sep="\t")
fin = pd.read_csv(HERE / "Ath_ABA_orthologs_final.tsv", sep="\t")
desc = itag_descriptions(ENS_GFF)

src = fin.sources.fillna("")
t = pd.DataFrame({
    "functional_group": fin.group.map(GROUP_LABEL),
    "arabidopsis_gene": fin.Symbol,
    "arabidopsis_locus": fin.ath,
    "tomato_gene": fin.sly.fillna("No ortholog"),
    "orthofinder": src.str.contains("OF").map({True: "Yes", False: "No"}),
    "ensembl_compara": src.str.contains("Compara").map({True: "Yes", False: "No"}),
    "orthofinder_relationship": fin.of_type.map(REL).fillna(""),
    "compara_relationship": fin.compara_type.str.replace("ortholog_", "").map(REL).fillna(""),
    "orthofinder_orthogroup": fin.of_orthogroup.fillna(""),
    "tomato_description": fin.sly.map(desc).fillna(""),
})
# order: that of the gene list; within each gene, pairs with both methods first
t["_g"] = t.arabidopsis_locus.map({a: i for i, a in enumerate(genes.gene_id)})
t["_n"] = (t.orthofinder == "Yes").astype(int) + (t.ensembl_compara == "Yes").astype(int)
t = t.sort_values(["_g", "_n", "orthofinder"], ascending=[True, False, False])
t.drop(columns=["_g", "_n"]).to_csv(OUT, sep="\t", index=False)

has = t[t.tomato_gene != "No ortholog"]
g_of = set(has[has.orthofinder == "Yes"].arabidopsis_locus)
g_c = set(has[has.ensembl_compara == "Yes"].arabidopsis_locus)
print(f"{OUT.name}: {len(t)} rows, {len(has)} pairs, {t.arabidopsis_locus.nunique()} Ath genes")
print(f"Venn genes — OF only {len(g_of - g_c)}, both {len(g_of & g_c)}, "
      f"Compara only {len(g_c - g_of)}, union {len(g_of | g_c)}")
