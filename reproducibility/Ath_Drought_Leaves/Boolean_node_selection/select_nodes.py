#!/usr/bin/env python3
"""
Node selection for the Boolean networks of the ABA response to drought in Arabidopsis and
tomato (genes from Song et al. 2016 and their tomato orthologs).

Filters and criteria (see README.md):
  F1  Ath gene with expression data and >=1 tomato ortholog with data
  F2  changes significantly (padj < 0.05, any |log2FC| > 0) at some time point in at
      least one species (in tomato: any of its orthologs)
  Tomato representative per Ath gene: significant > lower pattern frequency >
      supported by OrthoFinder and Compara > larger maximum |log2FC|
  Pattern: SSD binarization (more unique-pattern genes than WCSS in 3 of 4 tables and a
      tie in the other), control|drought concatenated (8 bits). Frequency computed over all
      the binarized genes of each species.
  Order: max(freq_Ath, freq_Sly), then freq_Ath + freq_Sly, then significant in both.
  No co-inclusion rule: within each species no two nodes with the same pattern, nor two
      Ath genes with the same tomato ortholog.
  Coverage: at most MAX_PER_GROUP nodes per functional group; MAX_NODES in total.
  Manual curation (after the automatic selection): MANUAL_INCLUDE are added and
      MANUAL_EXCLUDE removed (reasons in README.md); the no co-inclusion rule is checked
      again.

Outputs:
  candidates.tsv      the 74 candidates (F1+F2) with patterns, frequencies and rank
  selected_nodes.tsv  chosen nodes
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]   # reproducibility/
OUT = Path(__file__).resolve().parent
METHOD = "SSD"
MAX_NODES = 15
MAX_PER_GROUP = 2   # cover the whole pathway: no more than 2 nodes per functional group
# manual curation on top of the automatic selection (see README.md)
MANUAL_INCLUDE = ["HB7", "ABF3"]   # TFs of interest; only left out by the TF cap
MANUAL_EXCLUDE = ["DWA2", "ABA4"]  # pattern with almost no variation in both / in tomato


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

# F1: pairs with data in both species
o = orth[orth.arabidopsis_locus.isin(pa.index) & orth.tomato_gene.isin(ps.index)].copy()
o["both_methods"] = (o.orthofinder == "Yes") & (o.ensembl_compara == "Yes")
o = o.join(ps.add_prefix("sly_"), on="tomato_gene").join(ds.add_prefix("sly_"), on="tomato_gene")
o = o.join(pa.add_prefix("ath_"), on="arabidopsis_locus").join(da.add_prefix("ath_"),
                                                             on="arabidopsis_locus")
o["sly_sig"] = o.sly_sig.fillna(False).astype(bool)
o["ath_sig"] = o.ath_sig.fillna(False).astype(bool)

# F2: significant in at least one species (tomato: any ortholog)
any_sly = o.groupby("arabidopsis_locus").sly_sig.transform("any")
o = o[o.ath_sig | any_sly]

# tomato representative per Ath gene
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

# greedy selection with the no co-inclusion rule
used_a, used_s, used_t, keep, reason = set(), set(), set(), [], []
per_group = {}
for r in c.itertuples():
    why = []
    if r.ath_pattern in used_a:
        why.append("repeated Ath pattern")
    if r.sly_pattern in used_s:
        why.append("repeated Sly pattern")
    if r.tomato_gene in used_t:
        why.append("tomato ortholog already used")
    if not why and per_group.get(r.functional_group, 0) >= MAX_PER_GROUP:
        why.append(f"group {r.functional_group} full ({MAX_PER_GROUP})")
    ok = not why and sum(keep) < MAX_NODES
    keep.append(ok)
    reason.append("; ".join(why) if why else ("" if ok else "over the node limit"))
    if ok:
        used_a.add(r.ath_pattern); used_s.add(r.sly_pattern); used_t.add(r.tomato_gene)
        per_group[r.functional_group] = per_group.get(r.functional_group, 0) + 1
c["selected"], c["excluded_reason"] = keep, reason
c["auto_selected"] = c.selected

# manual curation
inc = c.arabidopsis_gene.isin(MANUAL_INCLUDE)
exc = c.arabidopsis_gene.isin(MANUAL_EXCLUDE)
assert inc.sum() == len(MANUAL_INCLUDE) and exc.sum() == len(MANUAL_EXCLUDE)
c.loc[inc, ["selected", "excluded_reason"]] = [True, "included manually"]
c.loc[exc, ["selected", "excluded_reason"]] = [False, "excluded manually"]
fin = c[c.selected]
for col in ("ath_pattern", "sly_pattern", "tomato_gene"):
    assert fin[col].is_unique, f"manual curation breaks the no co-inclusion rule ({col})"

cols = ["rank", "selected", "auto_selected", "excluded_reason", "functional_group", "arabidopsis_gene",
        "arabidopsis_locus", "tomato_gene", "n_tomato_orthologs", "both_methods", "sig_in",
        "ath_pattern", "ath_freq", "sly_pattern", "sly_freq", "max_freq",
        "ath_sig_times", "sly_sig_times", "ath_max_abs_lfc", "sly_max_abs_lfc"]
c[cols].to_csv(OUT / "candidates.tsv", sep="\t", index=False, float_format="%.2f")
c.loc[c.selected, cols[4:]].to_csv(OUT / "selected_nodes.tsv", sep="\t", index=False,
                                    float_format="%.2f")
print(f"F1: {orth[orth.arabidopsis_locus.isin(pa.index) & orth.tomato_gene.isin(ps.index)].arabidopsis_locus.nunique()}"
      f"  F1+F2: {len(c)}  automatic: {c.auto_selected.sum()}  final: {c.selected.sum()}")
