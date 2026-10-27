#!/usr/bin/env python3
"""
Selección de nodos para redes booleanas de la respuesta ABA a sequía en Arabidopsis y
tomate (genes de Song et al. 2016 y sus ortólogos en tomate).

Filtros y criterios (ver README.md):
  F1  gen Ath con datos de expresión y con >=1 ortólogo de tomate con datos
  F2  cambia significativamente (padj < 0.05, cualquier |log2FC| > 0) en algún tiempo
      en al menos una especie (en tomate: cualquiera de sus ortólogos)
  Representante de tomate por gen Ath: significativo > menor frecuencia de patrón >
      apoyado por OrthoFinder y Compara > mayor |log2FC| máximo
  Patrón: binarización SSD (más genes con patrón único que WCSS en 3 de 4 tablas y empate
      en la otra), control|sequía concatenados (8 bits). Frecuencia calculada sobre todos
      los genes binarizados de cada especie.
  Orden: max(freq_Ath, freq_Sly), luego freq_Ath + freq_Sly, luego significativo en ambas.
  Regla de no co-inclusión: en cada especie no entran dos nodos con el mismo patrón, ni
      dos genes Ath con el mismo ortólogo de tomate.
  Cobertura: como máximo MAX_PER_GROUP nodos por grupo funcional; total MAX_NODES.
  Curación manual (después de la selección automática): se suman MANUAL_INCLUDE y se sacan
      MANUAL_EXCLUDE (motivos en README.md); se verifica que se siga cumpliendo la regla de
      no co-inclusión.

Salidas:
  candidates.tsv      los 74 candidatos (F1+F2) con patrones, frecuencias y orden
  selected_nodes.tsv  nodos elegidos
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]   # BNI3_paper_backup/
OUT = Path(__file__).resolve().parent
METHOD = "SSD"
MAX_NODES = 15
MAX_PER_GROUP = 2   # cobertura de toda la vía: no más de 2 nodos por grupo funcional
# curación manual sobre la selección automática (ver README.md)
MANUAL_INCLUDE = ["HB7", "ABF3"]   # TF de interés; solo quedaban afuera por el tope de TF
MANUAL_EXCLUDE = ["DWA2", "ABA4"]  # patrón casi sin variación en ambas / en tomate


def patterns(sp):
    b = ROOT / f"{sp}_Drought_Leaves" / "Binarization"
    p = {}
    for cond in ("Control", "Drought"):
        f = b / f"ConteoNormalizado_mean_{cond}_aba_genes_binarized_{METHOD}_pattern.tsv"
        p[cond] = pd.read_csv(f, sep="\t", dtype={"Pattern": str}).set_index("Gene").Pattern
    df = pd.DataFrame({"ctrl": p["Control"], "wd": p["Drought"]})
    df["pattern"] = df.ctrl + "|" + df.wd
    df["freq"] = df.pattern.map(df.pattern.value_counts())
    return df


def dge(sp):
    t = pd.read_csv(ROOT / f"{sp}_Drought_Leaves" / "Selected_genes" /
                    "dotplot_aba_genes_table.tsv", sep="\t")
    sig = t[t.padj < 0.05]
    return pd.DataFrame({
        "sig": t.groupby("gene_id").padj.apply(lambda x: (x < 0.05).any()),
        "n_sig": sig.groupby("gene_id").size(),
        "max_abs_lfc": t.groupby("gene_id").log2FoldChange.apply(lambda x: x.abs().max()),
        "sig_times": sig.groupby("gene_id").comparison.apply(", ".join),
    }).fillna({"n_sig": 0, "sig_times": ""})


orth = pd.read_csv(ROOT / "TAIR_To_ITAG" / "Table_S_ABA_orthologs.tsv", sep="\t")
orth = orth[orth.tomato_gene != "No ortholog"]
pa, ps = patterns("Ath"), patterns("Sly")
da, ds = dge("Ath"), dge("Sly")

# F1: pares con datos en ambas especies
o = orth[orth.arabidopsis_locus.isin(pa.index) & orth.tomato_gene.isin(ps.index)].copy()
o["both_methods"] = (o.orthofinder == "Yes") & (o.ensembl_compara == "Yes")
o = o.join(ps.add_prefix("sly_"), on="tomato_gene").join(ds.add_prefix("sly_"), on="tomato_gene")
o = o.join(pa.add_prefix("ath_"), on="arabidopsis_locus").join(da.add_prefix("ath_"),
                                                             on="arabidopsis_locus")
o["sly_sig"] = o.sly_sig.fillna(False).astype(bool)
o["ath_sig"] = o.ath_sig.fillna(False).astype(bool)

# F2: significativo en al menos una especie (tomate: cualquier ortólogo)
any_sly = o.groupby("arabidopsis_locus").sly_sig.transform("any")
o = o[o.ath_sig | any_sly]

# representante de tomate por gen Ath
o = o.sort_values(["arabidopsis_locus", "sly_sig", "sly_freq", "both_methods",
                   "sly_max_abs_lfc"], ascending=[True, False, True, False, False])
o["n_tomato_orthologs"] = o.groupby("arabidopsis_locus").tomato_gene.transform("size")
c = o.drop_duplicates("arabidopsis_locus").copy()
c["sig_in"] = c.apply(lambda r: "both" if r.ath_sig and r.sly_sig
                      else ("Ath" if r.ath_sig else "Sly"), axis=1)
c["max_freq"] = c[["ath_freq", "sly_freq"]].max(axis=1)
c["sum_freq"] = c.ath_freq + c.sly_freq
c = c.sort_values(["max_freq", "sum_freq", "sig_in"],
                  key=lambda s: s.map({"both": 0, "Ath": 1, "Sly": 1}) if s.name == "sig_in"
                  else s).reset_index(drop=True)
c["rank"] = c.index + 1

# selección greedy con regla de no co-inclusión
used_a, used_s, used_t, keep, reason = set(), set(), set(), [], []
per_group = {}
for r in c.itertuples():
    why = []
    if r.ath_pattern in used_a:
        why.append("patrón Ath repetido")
    if r.sly_pattern in used_s:
        why.append("patrón Sly repetido")
    if r.tomato_gene in used_t:
        why.append("ortólogo de tomate ya usado")
    if not why and per_group.get(r.functional_group, 0) >= MAX_PER_GROUP:
        why.append(f"grupo {r.functional_group} completo ({MAX_PER_GROUP})")
    ok = not why and sum(keep) < MAX_NODES
    keep.append(ok)
    reason.append("; ".join(why) if why else ("" if ok else "fuera del cupo"))
    if ok:
        used_a.add(r.ath_pattern); used_s.add(r.sly_pattern); used_t.add(r.tomato_gene)
        per_group[r.functional_group] = per_group.get(r.functional_group, 0) + 1
c["selected"], c["excluded_reason"] = keep, reason
c["auto_selected"] = c.selected

# curación manual
inc = c.arabidopsis_gene.isin(MANUAL_INCLUDE)
exc = c.arabidopsis_gene.isin(MANUAL_EXCLUDE)
assert inc.sum() == len(MANUAL_INCLUDE) and exc.sum() == len(MANUAL_EXCLUDE)
c.loc[inc, ["selected", "excluded_reason"]] = [True, "incluido manualmente"]
c.loc[exc, ["selected", "excluded_reason"]] = [False, "excluido manualmente"]
fin = c[c.selected]
for col in ("ath_pattern", "sly_pattern", "tomato_gene"):
    assert fin[col].is_unique, f"curación manual rompe la regla de no co-inclusión ({col})"

cols = ["rank", "selected", "auto_selected", "excluded_reason", "functional_group", "arabidopsis_gene",
        "arabidopsis_locus", "tomato_gene", "n_tomato_orthologs", "both_methods", "sig_in",
        "ath_pattern", "ath_freq", "sly_pattern", "sly_freq", "max_freq",
        "ath_sig_times", "sly_sig_times", "ath_max_abs_lfc", "sly_max_abs_lfc"]
c[cols].to_csv(OUT / "candidates.tsv", sep="\t", index=False, float_format="%.2f")
c.loc[c.selected, cols[4:]].to_csv(OUT / "selected_nodes.tsv", sep="\t", index=False,
                                    float_format="%.2f")
print(f"F1: {orth[orth.arabidopsis_locus.isin(pa.index) & orth.tomato_gene.isin(ps.index)].arabidopsis_locus.nunique()}"
      f"  F1+F2: {len(c)}  automáticos: {c.auto_selected.sum()}  finales: {c.selected.sum()}")
