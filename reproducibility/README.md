# Reproducing the BNI3 paper

This folder holds the data, scripts and results behind the analyses in the BNI3 paper.
Every result here was produced with the code in this repository at release `v1.3.0`.
Only curated files are included. Raw reads, genome references and intermediate
per-experiment metrics are left out, and each section says where they come from.

| Analysis | Folder | Main result |
|---|---|---|
| DREAM4 size10 benchmark, BNI3 vs BoolNet (Figure 5) | [`DREAM4/`](DREAM4/) | `DREAM4/figures/Figure5_*` |
| Synthetic cell-cycle network: attractors and perturbations | [`synthetic_cellcycle/`](synthetic_cellcycle/) | `synthetic_cellcycle/results/` |
| ABA response to drought, *Arabidopsis thaliana* rosette leaves | [`Ath_Drought_Leaves/`](Ath_Drought_Leaves/) | `Boolean_Rules_Inference/{Control,Drought}/` |
| ABA response to drought, *Solanum lycopersicum* mature leaves | [`Sly_Drought_Leaves/`](Sly_Drought_Leaves/) | `Boolean_Rules_Inference/{Control,Drought}/` |
| Arabidopsis–tomato orthologs of the ABA pathway genes | [`TAIR_To_ITAG/`](TAIR_To_ITAG/) | `Table_S_ABA_orthologs.tsv` |

The species folders and `TAIR_To_ITAG/` carry their own READMEs (in Spanish), with the
details of each step.

## Environment

Install BNI3 as described in the [main README](../README.md). The paper was run with
Python 3.12.2, pandas 2.2.3, numpy 2.4.6, scipy 1.15.2, scikit-learn 1.6.1, deap 1.4.2,
geppy 0.1.3, numba 0.65.1, networkx 3.4.2, matplotlib 3.10.1 and seaborn 0.13.2.

The R steps (gold-standard scoring, BoolNet, DESeq2) need R with `dplyr`, `BoolNet`,
`DESeq2`, `ggplot2`, `BiocParallel` and `patchwork`.

Scripts that end in `.sh` resolve their paths from their own location and can be run
from any directory. The R scripts in the species folders and in
`synthetic_cellcycle/data/` use paths relative to their folder: run them from there.

The rule inference and the rule-combination evaluation are stochastic (GEP evolution;
hill climbing when the combination space is too large to enumerate). A rerun gives
networks of the same quality, but not necessarily the same rules.

---

## DREAM4 size10 (Figure 5)

**Data.** DREAM4 *in silico* network challenge, size 10 (Marbach et al., 2010).
`data/time_series/size10_X/` has, for each of the 5 networks, the 5 time series of
21 points (`size10_X_Y`) and the same series reduced to 5 points (`size10_X_Y_5times`).
`data/gold_standard/` has the corresponding gold standards.

**Pipeline** (`scripts/run_bni3_dream4.sh`), for each of the 50 series:

1. Binarization with SSD and WCSS (`BNI3_Binarizator.py`), with the behavior reviewer.
2. The method with the higher `Unique_Ratio` is kept (WCSS on ties).
3. Rule inference with `-p 50 -n_rep 10 -pop 50 -gen 50`, with the chained analysis
   and rule-combination evaluation.
4. Structural scoring against the gold standard (`scripts/eval_gold_standard.r`).

Self-regulation is excluded from the scoring on both sides, because the gold standard has
none: a gene is dropped from its own regulators and the diagonal is skipped, so 90
ordered pairs are scored per network. BoolNet is scored the same way.

**BoolNet.** `scripts/run_boolnet_dream4.r` reconstructs the same 50 series with
BoolNet (`reconstructNetwork`, method `bestfit`), and
`scripts/eval_gold_standard_boolnet.r` scores them.

**Results.**

| Path | Content |
|---|---|
| `results/BNI3/binarized/size10_X/<series>/` | SSD and WCSS matrices and the reviewer summary |
| `results/BNI3/rules/size10_X/rules_<METHOD>_<series>/` | candidate rules (`rules_by_gene.tsv`), selected network (`rules_by_gene_evaluated.tsv`), ranked combinations (`evaluation_results.tsv.gz`) and the evaluation log |
| `results/BNI3_size10_complete_evaluation_results.tsv` | BNI3 metrics per series |
| `results/BoolNet/size10_X/` | BoolNet binarized matrices and rules |
| `results/BoolNet_size10_complete_evaluation_results.tsv` | BoolNet metrics per series |
| `figures/Figure5_{A,B,C,D}_*.{svg,png}` | Figure 5 panels |

**Reproduce.**

```bash
# Figure 5 and its statistics, from the included results (seconds)
python3 reproducibility/DREAM4/scripts/make_figure5.py

# Gold-standard scoring of the included networks (reproduces the tables exactly)
Rscript reproducibility/DREAM4/scripts/eval_gold_standard.r
Rscript reproducibility/DREAM4/scripts/eval_gold_standard_boolnet.r

# Full rerun into DREAM4/rerun/ (~2.2 h per 21-point series and ~2 min per 5-point
# series with 50 processes; about two days in total)
nohup bash reproducibility/DREAM4/scripts/run_bni3_dream4.sh > run_bni3_dream4.log 2>&1 &
Rscript reproducibility/DREAM4/scripts/run_boolnet_dream4.r
```

`make_figure5.py` prints, for each panel and condition, n, the mean and median of each
method and the paired Wilcoxon signed-rank p-value (two-sided, continuity correction).

---

## Synthetic cell-cycle network

**Data.** `data/synthetic_expression_data.r` simulates expression from the mammalian
cell-cycle Boolean model distributed with BoolNet (`data(cellcycle)`; Fauré et al., 2006).
`data/reglas_gold_standard.tsv` lists the model's rules, and `rules/rules_by_gene.tsv`
holds the rules inferred by BNI3 (defaults: `-pop 100 -gen 100 -n_rep 10`).

**Reproduce.**

```bash
bash reproducibility/synthetic_cellcycle/scripts/run_attractors.sh
```

This computes the attractors for the wild type, the CycD knockout (`CycD:0`) and the CycD
overexpression (`CycD:1`), with their heatmaps and attractor diagrams
(`results/attractors/`), and draws the regulatory network (`results/network/`).

---

## ABA response to drought in Arabidopsis and tomato

Both species follow the same steps, one folder each:

| Folder | Step |
|---|---|
| `Counts/` | gene-level count matrix and sample table |
| `DGE/` | `DESeq2.R`: normalization (vst), mean per condition and WD vs CTRL contrasts |
| `Selected_genes/` | `selected_genes.R`: the ABA pathway genes (Song et al., 2016) and their expression |
| `Boolean_node_selection/` (Arabidopsis folder, both species) | `select_nodes.py`: the 15 network nodes, each an Arabidopsis gene paired with a tomato ortholog |
| `Binarization/` | SSD and WCSS matrices; SSD is the one used, restricted to the 15 nodes |
| `Boolean_Rules_Inference/{Control,Drought}/` | rules, ranked combinations, network figure and attractors |
| `BoolNet/` | `boolnet.r`: BoolNet (`bestfit`) on the same binarized matrices |

`TAIR_To_ITAG/` documents how the tomato orthologs were chosen.

**Data sources.**

- Tomato: BioProject PRJNA1125861 (21 paired-end runs). `Counts/run_kallisto.sh`
  quantified them against Ensembl Plants SL3.0 cDNA with kallisto 0.48.0, and
  `Counts/extraerCounts.r` built the gene matrix.
- Arabidopsis: <!-- TODO: GEO/SRA accession of the Arabidopsis drought time course -->

**Not included,** to keep the repository small. All of it is regenerated from the files
that are here:

- the full normalized matrices and the per-contrast DESeq2 tables (`DGE/`): run
  `DESeq2.R` from `DGE/` before `selected_genes.R`;
- reads, the kallisto index and TPM tables;
- per-experiment GEP metrics from the rule inference.

**Reproduce**, for each species:

```bash
cd reproducibility/Ath_Drought_Leaves          # or Sly_Drought_Leaves
(cd DGE && Rscript DESeq2.R)
(cd Selected_genes && Rscript selected_genes.R)
# binarization, node selection and inference: see the README of each folder
(cd BoolNet && Rscript boolnet.r)
```

Rule inference used the BNI3 defaults (`-pop 100 -gen 100 -n_rep 10`, regulator penalty
with K = 2) on the SSD matrix of the 15 nodes, `Binarization/*_binarized_SSD_selected_nodes_symbol.tsv`,
with the chained evaluation. Its settings (search strategy, cascade, number of
combinations) are at the top of each `evaluation_log.tsv`. Attractors were computed with
`3.Attractors/BNI3_Attractorator.py` on `rules_by_gene_evaluated.tsv` and the same
binarized matrix.
