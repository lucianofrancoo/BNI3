#!/usr/bin/env python3
"""
Combina tres fuentes de ortólogos Arabidopsis -> tomate (SL3.0) para los genes ABA.

  BLAST   : blast_fwd_ath_aba_vs_sly.tsv / blast_rev_sly_vs_tair12.tsv
            (de blast_ath_aba_to_sly.py; tblastx, e<=1e-10)
  Compara : compara62_ath_sly_homologies.tsv.gz  (Ensembl Plants release 62,
            Compara.115.protein_default, filas arabidopsis_thaliana -> solanum_lycopersicum)
            https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/tsv/ensembl-compara/homologies/arabidopsis_thaliana/Compara.115.protein_default.homologies.tsv.gz
  OrthoFinder : .../Results_sptree/Orthologues/.../
            Arabidopsis_thaliana__v__Solanum_lycopersicum.tsv (OrthoFinder 3.1.5, 10 especies,
            proteomas Ensembl Plants 62 + TAIR12, una isoforma por gen; -M msa, árbol de
            especies de referencia con -s)

Salida: Ath_ABA_orthologs_combined.tsv, una fila por par (gen Ath, gen Sly) propuesto por
alguna fuente. Pares BLAST incluidos: RBH, o top BLAST_TOP por gen Ath con bitscore
>= BLAST_FRAC * el mejor bitscore de ese gen Ath.
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
GENES = HERE.parent / "NCBI_To_ITAG" / "Ath" / "Ath_ABA_genes.tsv"
OUT = HERE / "Ath_ABA_orthologs_combined.tsv"
# corrida con árbol de especies de referencia (-s Orthofinder/species_tree_rooted.nwk)
OF_RESULTS = [HERE / "Orthofinder" / "results" / "Results_Sep24" / "WorkingDirectory"
              / "OrthoFinder" / "Results_sptree"]
BLAST_TOP = 5
BLAST_FRAC = 0.5
COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".split()


def sly_gene(s):
    return s.str.replace(r"\.\d+$", "", regex=True)


def ath_gene(s):
    return s.str.split(".").str[0].str.upper()


genes = pd.read_csv(GENES, sep="\t")
genes["gene_id"] = genes.gene_id.str.upper()

# --- BLAST
f = pd.read_csv(HERE / "blast_fwd_ath_aba_vs_sly.tsv", sep="\t", names=COLS)
f = f.assign(ath=ath_gene(f.qseqid), sly=sly_gene(f.sseqid))
f = f.sort_values("bitscore", ascending=False).drop_duplicates(["ath", "sly"])
f["blast_rank"] = f.groupby("ath").cumcount() + 1
f["blast_rel"] = f.bitscore / f.groupby("ath").bitscore.transform("max")
r = pd.read_csv(HERE / "blast_rev_sly_vs_tair12.tsv", sep="\t", names=COLS)
r = r.assign(sly=sly_gene(r.qseqid), ath=ath_gene(r.sseqid))
back = r.sort_values("bitscore", ascending=False).drop_duplicates("sly").set_index("sly").ath
f["blast_back_to"] = f.sly.map(back)
# RBH estricto: el mejor gen de tomate de A tiene como mejor gen Ath a A
f["blast_RBH"] = (f.blast_rank == 1) & (f.blast_back_to == f.ath)
blast = f[f.blast_RBH | ((f.blast_rank <= BLAST_TOP) & (f.blast_rel >= BLAST_FRAC))]
blast = blast[["ath", "sly", "blast_rank", "bitscore", "pident", "blast_back_to", "blast_RBH"]]
blast = blast.rename(columns={"bitscore": "blast_bitscore", "pident": "blast_pident"})

# --- Compara
c = pd.read_csv(HERE / "compara62_ath_sly_homologies.tsv.gz", sep="\t")
c = c[c.gene_stable_id.isin(genes.gene_id)]
c = c.rename(columns={"gene_stable_id": "ath", "homology_gene_stable_id": "sly",
                      "homology_type": "compara_type", "identity": "compara_identity",
                      "is_high_confidence": "compara_high_conf"})
c = c[["ath", "sly", "compara_type", "compara_identity", "compara_high_conf"]]

# --- OrthoFinder (el último Results_*)
of_file = next(OF_RESULTS[-1].rglob("Arabidopsis_thaliana__v__Solanum_lycopersicum.tsv"))
o = pd.read_csv(of_file, sep="\t")
strip = lambda x: x.str.strip().str.replace(r"^(Arabidopsis_thaliana|Solanum_lycopersicum)_", "", regex=True)
o["ath"] = o.Arabidopsis_thaliana.str.split(",")
o["sly"] = o.Solanum_lycopersicum.str.split(",")
o["n_ath"] = o.ath.str.len()
o["n_sly"] = o.sly.str.len()
o["of_type"] = o.n_ath.map(lambda n: "one" if n == 1 else "many") + "2" + \
    o.n_sly.map(lambda n: "one" if n == 1 else "many")
o = o.explode("ath").explode("sly")
o["ath"], o["sly"] = strip(o.ath), strip(o.sly)
o = o[o.ath.isin(genes.gene_id)][["ath", "sly", "Orthogroup", "of_type"]]
o = o.rename(columns={"Orthogroup": "of_orthogroup"})

m = (blast.merge(c, on=["ath", "sly"], how="outer")
          .merge(o, on=["ath", "sly"], how="outer"))
m["BLAST"] = m.blast_rank.notna()
m["Compara"] = m.compara_type.notna()
m["OrthoFinder"] = m.of_type.notna()
SOURCES = ["OrthoFinder", "Compara", "BLAST"]
m["n_sources"] = m[SOURCES].sum(axis=1)
# posición del par en el BLAST aunque no haya pasado el filtro
rank_all = f.set_index(["ath", "sly"]).blast_rank
m["blast_rank"] = [rank_all.get((a, s)) for a, s in zip(m.ath, m.sly)]

m = m.merge(genes.rename(columns={"gene_id": "ath"})[["ath", "Symbol", "group"]], on="ath")
m = m.sort_values(["group", "Symbol", "n_sources", "blast_rank"],
                  ascending=[True, True, False, True])
cols = ["group", "Symbol", "ath", "sly", "n_sources"] + SOURCES + ["of_type", "of_orthogroup",
        "blast_rank", "blast_bitscore", "blast_pident", "blast_RBH", "blast_back_to",
        "compara_type", "compara_identity", "compara_high_conf"]
m[cols].to_csv(OUT, sep="\t", index=False, float_format="%.1f")

# --- resumen
cov = pd.DataFrame({s: genes.gene_id.isin(m[m[s]].ath) for s in SOURCES})
cov.index = genes.Symbol
print("genes Ath con al menos un ortólogo por fuente:", cov.sum().to_dict())
print("sin ninguno:", cov.index[~cov.any(axis=1)].tolist())
print("pares:", len(m), " por n_sources:", m.n_sources.value_counts().sort_index().to_dict())
print("OrthoFinder vs Compara (pares OrthoFinder):", m[m.OrthoFinder].Compara.value_counts().to_dict())
print("Compara vs OrthoFinder (pares Compara):", m[m.Compara].OrthoFinder.value_counts().to_dict())
print("Compara vs BLAST (pares Compara con blast_rank):",
      m[m.Compara].blast_rank.fillna(-1).astype(int).clip(upper=6).value_counts().sort_index().to_dict())

# --- tabla final: ortólogos = pares propuestos por OrthoFinder o Compara (métodos
#     basados en árboles). El RBH no define ortólogos; solo se reporta como validación.
OUT_FINAL = HERE / "Ath_ABA_orthologs_final.tsv"
tree = m[m.OrthoFinder | m.Compara].copy()
tree["sources"] = tree.apply(lambda r: "+".join(
    lab for lab, col in [("OF", "OrthoFinder"), ("Compara", "Compara")] if r[col]), axis=1)
fin = genes[["group", "Symbol", "gene_id"]].rename(columns={"gene_id": "ath"}).merge(
    tree.drop(columns=["group", "Symbol"]), on="ath", how="left")
fin["has_ortholog"] = fin.sly.notna()
fin = fin.sort_values(["group", "Symbol", "OrthoFinder", "Compara"],
                      ascending=[True, True, False, False])
fin[["group", "Symbol", "ath", "has_ortholog", "sly", "sources", "of_type",
     "compara_type", "of_orthogroup"]].to_csv(OUT_FINAL, sep="\t", index=False)
rbh = m[m.blast_RBH.eq(True)]
print("final:", fin.drop_duplicates("ath").has_ortholog.value_counts().to_dict(), "genes;",
      fin.sly.notna().sum(), "pares;", fin.sly.nunique(), "genes de tomate distintos")
print(f"validación: {(rbh.OrthoFinder | rbh.Compara).sum()} de {len(rbh)} pares RBH "
      "también los dan OrthoFinder o Compara")
