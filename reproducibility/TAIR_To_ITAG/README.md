# Tomato orthologs of Arabidopsis ABA genes

Analysis run on 23–24 September 2026.

## Goal

For a list of 95 ABA-related Arabidopsis genes (synthesis, degradation, transport,
perception, signaling, transcription factors and target genes), find **all** their tomato
orthologs, as Solyc IDs of the ITAG annotation on the **SL3.0** assembly.

All orthologs are kept, not just one per gene, because most of these families expanded
independently in each lineage after Arabidopsis and tomato split (about 110 million years
ago). The expected relationship is therefore often one-to-many or many-to-many, not
one-to-one.

## Summary

- **86 of the 95 Arabidopsis genes have at least one tomato ortholog.** That is 153
  Arabidopsis–tomato pairs, with 116 distinct tomato genes.
- Genes covered by each method:
  - OrthoFinder: 72.
  - Ensembl Compara: 71.
  - Both methods: 57.
  - OrthoFinder only: 15.
  - Compara only: 14.
- **No ortholog (9):**
  - NFYC2, BG1, CYP707A2, ABCG30, TMAC2 and COR15.
  - DREB2A, PYL10 and RCAR3: supported only by BLAST, which is not used as a criterion.
- The supplementary table is `Table_S_ABA_orthologs.tsv`.
- The figures are `fig_orthologs_summary.pdf` and `fig_orthologs_detail.pdf`.

**Definition used:** a tomato gene is an ortholog of an Arabidopsis gene if **OrthoFinder or
Ensembl Compara** infers it. Both methods are based on gene trees.

**Included in this repository:** the scripts, `compara62_ath_sly_homologies.tsv.gz`, the
combined and final ortholog tables, the supplementary table and the figures. The sequence
files, the raw BLAST hits and the OrthoFinder run are not included, because of their size.
Section 1 says where each input comes from.

---

## 1. Input data

### Provided

| File | Content | md5 |
|---|---|---|
| `Ath_ABA_genes.tsv` | Curated list of the 95 ABA genes: `Symbol`, `gene_id` (AGI), `hormone`, `group`. The same genes are in `../Ath_Drought_Leaves/Selected_genes/aba_genes.tsv` | — |
| `TAIR12_cdna.fa` | Arabidopsis cDNA, TAIR12 (54,906 transcripts) | 8d6cc1c8ab263af0833353eca3228b5c |
| `Orthofinder/TAIR12_pep.fa.gz` | Arabidopsis proteins, TAIR12 (36,693) | 4c9889f5e046b53ad9d7e4698de0500a |
| `Sly3.0_cdna.fa` | Tomato SL3.0 cDNA in Ensembl format (35,825), decompressed from the original `Sly3.0_cdna.bgz` | — |
| `Orthofinder/SL3.0.fa.bgz` | Tomato SL3.0 proteins in Ensembl format (34,429, one per gene) | 45a34c22821cc311edd69d2a704f04ba |

The tomato file headers follow the Ensembl Plants format, with Solyc IDs of the ITAG3
annotation on SL3.0.

### Downloaded (24 September 2026)

**Ensembl Compara**, Ensembl Plants homologies, **release 62**:
`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/tsv/ensembl-compara/homologies/arabidopsis_thaliana/Compara.115.protein_default.homologies.tsv.gz`

- Only the Arabidopsis → *Solanum lycopersicum* rows were kept, 28,220 pairs, in
  `compara62_ath_sly_homologies.tsv.gz` (md5 8756d4c39f115d87c1199cc708338837).
- Release 62 is used because it is **the last Ensembl Plants release with tomato on SL3.0**.
  From release 63 on, tomato is on SL4.0 (`solanum_lycopersicum_gca000188115v5cm`), with
  different IDs.

**Proteomes of 8 more species for OrthoFinder**, from Ensembl Plants release 62:
`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/<species>/pep/`.
All checksums were verified against the Ensembl `CHECKSUMS` files.

| Species | File | Genes (one isoform) |
|---|---|---|
| *Brassica rapa* | Brassica_rapa.Brapa_1.0.pep.all.fa.gz | 41,018 |
| *Solanum tuberosum* | Solanum_tuberosum.SolTub_3.0.pep.all.fa.gz | 39,021 |
| *Capsicum annuum* | Capsicum_annuum.ASM51225v2.pep.all.fa.gz | 35,845 |
| *Nicotiana attenuata* | Nicotiana_attenuata.NIATTr2.pep.all.fa.gz | 33,320 |
| *Coffea canephora* | Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa.gz | 25,574 |
| *Vitis vinifera* | Vitis_vinifera.ASM3070453v1.pep.all.fa.gz | 37,480 |
| *Populus trichocarpa* | Populus_trichocarpa.Pop_tri_v4.pep.all.fa.gz | 34,699 |
| *Oryza sativa* | Oryza_sativa.IRGSP-1.0.pep.all.fa.gz | 35,804 |

The tomato gene descriptions come from the Ensembl Plants 62 GFF
(`Solanum_lycopersicum.SL3.0.62.gff3.gz`).

---

## 2. What was done and why

### 2.1 Concepts

- **Homologs:** genes with a common ancestor.
- **Orthologs:** homologs that diverged through **speciation**.
- **Paralogs:** homologs that diverged through **duplication**.

A good alignment (low e-value) shows that two genes are homologs, but does not tell
orthologs from paralogs. That requires reconstructing the family's history: a gene tree
reconciled with the species tree.

### 2.2 First attempt: BLAST reciprocal best hits (`blast_ath_aba_to_sly.py`)

- **What:**
  - The cDNAs of the 95 genes were extracted (`Ath_ABA_genes_cdna.fa`, 120 transcripts,
    all isoforms, with the symbol in the header).
  - **tblastx** search (BLAST+ 2.12.0) against tomato cDNA, with e-value ≤ 1e-10.
  - Reverse search: the 5 best tomato genes of each gene were searched against all of TAIR12.
  - Everything was summarized at gene level, keeping the best HSP of each pair.
- **Why tblastx and not blastn:** between species this divergent, DNA is saturated with
  third-codon-position changes. Comparing translations keeps much more signal.
- **Reciprocal best hit (RBH):** the best tomato gene of gene A has A as its best
  Arabidopsis gene. This gives 38 pairs.
- **Limitation:** RBH picks a single gene on each side. When a family duplicated in one
  lineage (SnRK2, PYL, clade A PP2C, ABF), the other copies, which are orthologs too, are
  missed. If a copy was lost, RBH can take a paralog for an ortholog.
- **Final use:** **validation only.** 31 of the 38 RBH pairs (82%) are also found by
  OrthoFinder or Compara. The other 7 are one gene with no other support (DREB2A) and six
  pairs that the tree-based methods reject, probably paralogs.

### 2.3 DAVID (discarded)

- The *DAVID Ortholog* tool does not infer orthology. It converts IDs with **OMA** or
  **Ensembl Compara** tables.
- Its output had 53 genes with **a single** tomato gene each, and missed core genes such as
  ABI1, ABI2, RCAR1 and PYL1.
- It was discarded because it is a reduced subset of Compara, which is used directly.

### 2.4 Ensembl Compara (release 62)

- Ensembl's comparative genomics system.
- It builds gene trees for about 100 plant species and reconciles them with the species
  tree, classifying each pair as ortholog (`one2one`, `one2many`, `many2many`) or paralog.
- All Arabidopsis → tomato orthologs of the 95 genes were taken: 71 genes with an
  ortholog, 106 pairs.

### 2.5 OrthoFinder

OrthoFinder was run because it allows control over the annotation (SL3.0), the species and
the versions, and because it is the most widely used phylogenetic orthology method today.

1. **Proteomes:**
   - One protein per gene, the longest isoform, as OrthoFinder recommends, so isoforms do
     not show up as false duplications.
   - Headers are the gene ID.
2. **Species:** 10 in total. With only two species the trees cannot be rooted and
   duplications cannot be placed on a lineage.
   - *Arabidopsis* + *Brassica rapa*: separate duplications specific to the Brassicaceae.
   - Tomato + potato, pepper and *Nicotiana*: separate the Solanaceae duplications,
     including their whole-genome triplication.
   - Coffee: another asterid.
   - Grapevine and poplar: rosids and relatives.
   - Rice: outgroup.
3. **Run:**
   - OrthoFinder 3.1.5 with DIAMOND 2.2.8 (all against all).
   - MCL 22.282 to build the orthogroups.
   - Multiple alignments with FAMSA 2.4.1 and trees with FastTree 2.2.0 (`-M msa`).
   - 32 threads for the searches and 16 for the analysis.
4. **Species tree:**
   - The tree inferred by OrthoFinder grouped grapevine with the asterids and left the
     Brassicaceae outside. This is a known artifact (long-branch attraction), because the
     Brassicaceae evolve fast and grapevine slowly.
   - It was rerun with the accepted tree (option `-s`), reusing the DIAMOND searches (`-b`):
     ```
     (Oryza,((Vitis,(Populus,(Arabidopsis,Brassica))),(Coffea,(Nicotiana,(Capsicum,(S.tuberosum,S.lycopersicum))))))
     ```
   - For these 95 genes, both runs give **the same 125 pairs**. The run with the correct
     tree is reported.
5. **Result:** 72 genes with an ortholog and 125 pairs, from
   `Orthologues/Orthologues_Arabidopsis_thaliana/Arabidopsis_thaliana__v__Solanum_lycopersicum.tsv`.

### 2.6 Combination (`combine_orthologs.py`)

- Merges the BLAST, Compara and OrthoFinder pairs into `Ath_ABA_orthologs_combined.tsv`,
  with the information from each source. The BLAST columns are kept for reference only.
- The final table (`Ath_ABA_orthologs_final.tsv`) takes as orthologs the pairs from
  **OrthoFinder or Compara**.
- **Why the union and not the intersection:** each method fails on different families.
  - Compara assigns no orthologs to ABI1, ABI2, AHG1, AHG3, RCAR1, PYL3, PYL5 or PYL7,
    which OrthoFinder resolves.
  - OrthoFinder leaves NCED3, NCED9, AAO3, ABF1 and ABF4 without an ortholog, which
    Compara assigns.
  - The intersection would lose core genes of the pathway. The union covers 86 of 95.

### 2.7 Figures and supplementary table

- `plot_orthologs.py` draws both figures:
  - `fig_orthologs_summary`:
    - Panel A: Venn diagram of the Arabidopsis genes with an ortholog, by method.
    - Panel B: coverage by functional group, for the union and for each method.
  - `fig_orthologs_detail`: every pair with the method that supports it.
- `make_supp_table.py` writes `Table_S_ABA_orthologs.tsv`, with one row per pair and one
  row per gene without an ortholog. Columns:
  - `functional_group`, `arabidopsis_gene`, `arabidopsis_locus`, `tomato_gene`.
  - `orthofinder`, `ensembl_compara` (Yes/No).
  - `orthofinder_relationship`: Arabidopsis:tomato direction.
  - `compara_relationship`: type according to Compara.
  - `orthofinder_orthogroup`.
  - `tomato_description`: the ITAG description exists for only 39 of the 153 pairs.

---

## 3. Cases to check before interpreting

- **ABI1 and ABI2 → `Solyc07g040990.3`, OrthoFinder only.** The best BLAST hit was
  `Solyc08g062650.3`. The orthogroup tree is worth a look.
- **NCED3 and NCED9 → `Solyc07g056570.1` (SlNCED1, *notabilis*), Compara only.**
  OrthoFinder assigns that gene to NCED2 and NCED5, and infers that the tomato lineage lost
  the NCED3-type copy. That part of the (FastTree) tree is poorly resolved.
- **ABCG40 (9 tomato genes) and ABCG31 (6):** these families are heavily expanded and the
  list probably includes copies with different functions.
- **RD29A → `Solyc01g009660.2`, OrthoFinder only.** RD29A is considered nearly exclusive to
  the Brassicaceae.

The family trees are in OrthoFinder's `Resolved_Gene_Trees/Resolved_Gene_Trees.txt`
(Newick format; they can be viewed, for example, with iTOL).

---

## 4. Files

| File | What it is |
|---|---|
| `blast_ath_aba_to_sly.py` | Forward and reverse tblastx search (RBH) |
| `combine_orthologs.py` | Combines the sources and writes the final table |
| `plot_orthologs.py` | Figures |
| `make_supp_table.py` | Supplementary table |
| `compara62_ath_sly_homologies.tsv.gz` | Arabidopsis–tomato homologies from Compara 62 |
| `Ath_ABA_orthologs_combined.tsv` | All sources, pair by pair |
| `Ath_ABA_orthologs_final.tsv` | Final orthologs (OrthoFinder or Compara) |
| `Table_S_ABA_orthologs.tsv` | Supplementary table |
| `fig_orthologs_summary.pdf/.png`, `fig_orthologs_detail.pdf/.png` | Figures |

## 5. Reproduce

The scripts expect the inputs of section 1 (sequences, gene list, tomato GFF) and the
OrthoFinder results, which are not in this repository.

```bash
# 1. BLAST (needs BLAST+ on the PATH)
python3 blast_ath_aba_to_sly.py
# 2. OrthoFinder 3.1.5 on the 10 proteomes, with the species tree above
#    (about 35 min with 32 threads)
# 3. Combination, figures and supplementary table
python3 combine_orthologs.py
python3 plot_orthologs.py
python3 make_supp_table.py
```

---

## 6. Downstream use of these orthologs

`Table_S_ABA_orthologs.tsv` is the input of two analyses:

1. **Expression of the tomato orthologs under drought**
   (`../Sly_Drought_Leaves/Selected_genes/selected_genes.R`).
   - The same analysis as in Arabidopsis (`../Ath_Drought_Leaves/Selected_genes/`),
     applied to the 116 tomato genes with an ortholog.
   - Of those 116, 97 have DGE data. The other 19 are barely expressed in leaves and were
     removed by the DESeq2 filter.
2. **Node selection for the Boolean networks of both species**
   (`../Ath_Drought_Leaves/Boolean_node_selection/`).
   - Keeps the Arabidopsis genes that have a tomato ortholog, change in at least one
     species and have rare binary patterns. One representative ortholog is chosen per gene.
   - The final list has 15 pairs: MYB44, HB6, HB7, ABF3, NCED5, CYP707A3, CYP707A4, AIT1,
     ABCG31, RCAR1, PYL1, AHG1, HAI3, SnRK2.8 and CPK1.

## 7. Methods summary

> Tomato orthologs of 95 Arabidopsis ABA-related genes were identified with two gene-tree
> based methods. First, OrthoFinder v3.1.5 (Emms & Kelly, 2019) was run on the proteomes
> of ten species (longest isoform per gene): *Arabidopsis thaliana* (TAIR12), *Solanum
> lycopersicum* (SL3.0, ITAG3), and, from Ensembl Plants release 62, *S. tuberosum*,
> *Capsicum annuum*, *Nicotiana attenuata*, *Coffea canephora*, *Vitis vinifera*,
> *Populus trichocarpa*, *Brassica rapa* and *Oryza sativa*. We used DIAMOND searches,
> MSA-based gene trees (FAMSA, FastTree) and a user-supplied rooted species tree. Second,
> Arabidopsis–tomato orthologs were retrieved from Ensembl Compara (Ensembl Plants
> release 62, the last release based on SL3.0). A tomato gene was considered an ortholog
> when inferred by either method, retaining one-to-many and many-to-many relationships.
> As an independent check, 31 of 38 reciprocal best tBLASTX hits were also recovered by
> at least one of these methods.

References:
- OrthoFinder: Emms & Kelly (2019) *Genome Biology* 20:238; Emms et al. (2025) bioRxiv
  doi:10.1101/2025.07.15.664860.
- Ensembl Plants / Compara: cite the current Ensembl Genomes publication.
