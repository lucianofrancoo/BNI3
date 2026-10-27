#!/usr/bin/env python3
"""
Ortólogos en tomate (SL3.0) de genes ABA de Arabidopsis por mejor hit recíproco (RBH).

Entradas (en este directorio salvo la lista de genes):
  Ath_ABA_genes_cdna.fa   cDNA de los genes ABA (TAIR12, todas las isoformas)
  Sly3.0_cdna.fa          cDNA de tomate SL3.0 (Ensembl Plants)
  TAIR12_cdna.fa          cDNA completo de Arabidopsis (para la búsqueda recíproca)
  ../NCBI_To_ITAG/Ath/Ath_ABA_genes.tsv  Symbol, gene_id, hormone, group

Método:
  1. tblastx (traducción en los 6 marcos de query y subject; más sensible que blastn
     entre especies tan divergentes) Ath ABA -> Sly, evalue <= EVALUE.
  2. Los hits se colapsan a nivel de gen (Ath: AT1G49720.1 -> AT1G49720;
     Sly: Solyc01g067380.3.1 -> Solyc01g067380.3) quedándose con el HSP de mayor
     bitscore de cada par de genes.
  3. Recíproco: los cDNA de los TOP_N mejores genes de tomato por gen Ath se buscan
     con tblastx contra todo TAIR12; mejor gen Ath de cada gen de tomate.
  4. RBH = el mejor gen de tomate de A tiene como mejor gen de Arabidopsis a A.

Salidas:
  blast_fwd_ath_aba_vs_sly.tsv   hits crudos Ath->Sly (outfmt 6 + qlen slen)
  blast_rev_sly_vs_tair12.tsv    hits crudos Sly->TAIR12
  Ath_ABA_orthologs_sly.tsv      una fila por gen Ath (ver columnas en README/abajo)
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
QUERY = HERE / "Ath_ABA_genes_cdna.fa"
SLY = HERE / "Sly3.0_cdna.fa"
TAIR = HERE / "TAIR12_cdna.fa"
GENES = HERE.parent / "NCBI_To_ITAG" / "Ath" / "Ath_ABA_genes.tsv"
DB = HERE / "blastdb"
FWD = HERE / "blast_fwd_ath_aba_vs_sly.tsv"
REV = HERE / "blast_rev_sly_vs_tair12.tsv"
OUT = HERE / "Ath_ABA_orthologs_sly.tsv"

EVALUE = 1e-10
TOP_N = 5          # genes de tomate por gen Ath que se prueban en la vuelta
THREADS = 32
FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
COLS = FMT.split()[1:]


def run(cmd):
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    subprocess.run(list(map(str, cmd)), check=True)


def makedb(fa):
    db = DB / fa.stem
    if not Path(str(db) + ".nsq").exists() and not Path(str(db) + ".00.nsq").exists():
        DB.mkdir(exist_ok=True)
        run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db, "-parse_seqids"])
    return db


def tblastx(query, db, out):
    if not out.exists():
        run(["tblastx", "-query", query, "-db", db, "-out", out, "-outfmt", FMT,
             "-evalue", EVALUE, "-max_target_seqs", 50, "-num_threads", THREADS])
    return pd.read_csv(out, sep="\t", names=COLS)


def ath_gene(s):
    return s.str.split(".").str[0].str.upper()


def sly_gene(s):
    return s.str.replace(r"\.\d+$", "", regex=True)


def best_per_gene_pair(h, qgene, sgene):
    h = h.assign(qgene=qgene(h.qseqid), sgene=sgene(h.sseqid))
    return (h.sort_values("bitscore", ascending=False)
             .drop_duplicates(["qgene", "sgene"]).reset_index(drop=True))


def main():
    genes = pd.read_csv(GENES, sep="\t")
    genes["gene_id"] = genes.gene_id.str.upper()

    fwd = best_per_gene_pair(tblastx(QUERY, makedb(SLY), FWD), ath_gene, sly_gene)
    fwd["rank"] = fwd.groupby("qgene").cumcount() + 1

    # vuelta: cDNA de los TOP_N genes de tomato de cada gen Ath contra todo TAIR12
    top = fwd[fwd["rank"] <= TOP_N]
    rev_ids = set(top.sseqid)  # transcrito con el mejor HSP de cada par
    rev_q = HERE / "blast_rev_query.fa"
    run(["blastdbcmd", "-db", makedb(SLY), "-entry", ",".join(sorted(rev_ids)),
         "-out", rev_q])
    rev = best_per_gene_pair(tblastx(rev_q, makedb(TAIR), REV), sly_gene, ath_gene)
    rev_q.unlink()
    rev_best = rev.drop_duplicates("qgene").set_index("qgene")

    rows = []
    for g in genes.itertuples():
        h = fwd[fwd.qgene == g.gene_id]
        row = dict(Symbol=g.Symbol, gene_id=g.gene_id, group=g.group)
        if h.empty:
            rows.append(row | dict(n_sly_hits=0))
            continue
        b = h.iloc[0]
        back = rev_best.sgene.get(b.sgene)
        rbh_hits = [s for s in h[h["rank"] <= TOP_N].sgene
                    if rev_best.sgene.get(s) == g.gene_id]
        others = h[(h["rank"] > 1) & (h["rank"] <= TOP_N)]
        row |= dict(
            best_sly=b.sgene, pident=b.pident, aln_len=b.length, evalue=b.evalue,
            bitscore=b.bitscore, best_sly_back_to=back,
            best_back_symbol=genes.set_index("gene_id").Symbol.get(back, ""),
            RBH=back == g.gene_id,
            other_rbh_sly=";".join(s for s in rbh_hits if s != b.sgene),
            n_sly_hits=len(h),
            next_hits=";".join(f"{r.sgene}:{r.bitscore:.0f}" for r in others.itertuples()))
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"{len(res)} genes Ath; con hit: {(res.n_sly_hits > 0).sum()}; "
          f"RBH: {res.RBH.fillna(False).astype(bool).sum()}", file=sys.stderr)


if __name__ == "__main__":
    main()
