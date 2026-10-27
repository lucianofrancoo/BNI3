#!/usr/bin/env python3
"""
Tomato (SL3.0) orthologs of Arabidopsis ABA genes by reciprocal best hit (RBH).

Inputs (in this folder except the gene list; not included in the repository):
  Ath_ABA_genes_cdna.fa   cDNA of the ABA genes (TAIR12, all isoforms)
  Sly3.0_cdna.fa          tomato SL3.0 cDNA (Ensembl Plants)
  TAIR12_cdna.fa          full Arabidopsis cDNA (for the reciprocal search)
  ../NCBI_To_ITAG/Ath/Ath_ABA_genes.tsv  Symbol, gene_id, hormone, group

Method:
  1. tblastx (translation in the 6 frames of query and subject; more sensitive than blastn
     between species this divergent) Ath ABA -> Sly, evalue <= EVALUE.
  2. Hits are collapsed to gene level (Ath: AT1G49720.1 -> AT1G49720;
     Sly: Solyc01g067380.3.1 -> Solyc01g067380.3), keeping the HSP with the highest
     bitscore of each gene pair.
  3. Reciprocal: the cDNAs of the TOP_N best tomato genes per Ath gene are searched with
     tblastx against all of TAIR12; best Ath gene of each tomato gene.
  4. RBH = the best tomato gene of A has A as its best Arabidopsis gene.

Outputs:
  blast_fwd_ath_aba_vs_sly.tsv   raw Ath->Sly hits (outfmt 6 + qlen slen)
  blast_rev_sly_vs_tair12.tsv    raw Sly->TAIR12 hits
  Ath_ABA_orthologs_sly.tsv      one row per Ath gene (columns in README/below)
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
TOP_N = 5          # tomato genes per Ath gene tested in the reverse search
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

    # reverse: cDNA of the TOP_N tomato genes of each Ath gene against all of TAIR12
    top = fwd[fwd["rank"] <= TOP_N]
    rev_ids = set(top.sseqid)  # transcript with the best HSP of each pair
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
    print(f"{len(res)} Ath genes; with a hit: {(res.n_sly_hits > 0).sum()}; "
          f"RBH: {res.RBH.fillna(False).astype(bool).sum()}", file=sys.stderr)


if __name__ == "__main__":
    main()
