# Node selection for the Boolean networks: ABA response to drought in Arabidopsis and tomato

Analysis of 24 September 2026. Script: `select_nodes.py`.

## Goal

Choose a short list (≤ 15) of ABA pathway genes to build Boolean networks in Arabidopsis
and tomato, and compare the drought response between the two species. Each node is a pair
made of an Arabidopsis gene and a tomato ortholog.

The Arabidopsis genes are the 95 ABA pathway genes taken from Song et al. (2016). The
tomato orthologs are described in `../../TAIR_To_ITAG/README.md`.

## Data

| What | File |
|---|---|
| Arabidopsis–tomato orthologs | `../../TAIR_To_ITAG/Table_S_ABA_orthologs.tsv` (OrthoFinder or Ensembl Compara) |
| Arabidopsis DGE (rosette; WD vs CTRL at 7200, 15840 and 20160) | `../Selected_genes/dotplot_aba_genes_table.tsv` |
| Tomato DGE (mature leaf; WD vs CTRL at d4, d5, d8 and d10) | `../../Sly_Drought_Leaves/Selected_genes/dotplot_aba_genes_table.tsv` |
| Binary patterns | `../../{Ath,Sly}_Drought_Leaves/Binarization/ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_pattern.tsv` |

## Criteria, in order

1. **F1, orthology:** the Arabidopsis gene has expression data and at least one tomato
   ortholog with data. **76 of 95** genes remain.
2. **F2, change:** padj < 0.05 (any |log2FC| > 0) at one time point or more, in at least
   one species. In tomato any of the orthologs counts. **74** remain: ABI1 and AFP2 drop.
3. **Binarization method: SSD.** The method with more unique-pattern genes was chosen.

   | Table | SSD | WCSS |
   |---|---|---|
   | Arabidopsis, control | 3 | 1 |
   | Arabidopsis, drought | 3 | 3 |
   | Tomato, control | 1 | 1 |
   | Tomato, drought | 2 | 1 |

   SSD is used on ties.
4. **8-bit pattern: control|drought.** Each condition was binarized separately (4 time
   points each). They are concatenated because the network models both conditions and
   because 4 bits give very little resolution: in Arabidopsis drought, 51 of 89 genes share
   the pattern `0001`. With 8 bits there are 35 distinct patterns in Arabidopsis (21 unique)
   and 51 in tomato (36 unique).
   - The **frequency** is the number of genes with that pattern among all the binarized
     genes of the species: 89 in Arabidopsis and 97 in tomato.
5. **Representative ortholog:** one per Arabidopsis gene. Preferred, in this order:
   - the one that changes significantly;
   - the one with the less frequent pattern;
   - the one found by both OrthoFinder and Compara;
   - the one with the largest |log2FC|.
6. **Priority order:**
   - First, the higher frequency of the two species, so the pattern that is rare in
     **both** is preferred.
   - Then, the sum of the two frequencies.
   - Last, whether it changes significantly in both species.
7. **No co-inclusion rule:** within each species, no two nodes with the same pattern are
   included. If two genes share a pattern, the Boolean network cannot tell which one causes
   a downstream change. Nor are two Arabidopsis genes with the same tomato ortholog included.
8. **Pathway coverage:** at most 2 nodes per functional group, 15 in total.
   - Without this cap, the frequency order filled the list with receptors and synthesis
     genes and left out the clade A PP2Cs, a central link of the PYL → PP2C → SnRK2 axis.

## Result: 15 nodes

### Automatic selection and manual curation

The automatic selection (criteria 1–8) gave 15 nodes. Two manual changes were then made,
recorded in `MANUAL_INCLUDE` and `MANUAL_EXCLUDE` in `select_nodes.py`. The script checks
that the final selection still satisfies the no co-inclusion rule.

- **Included: HB7 and ABF3.** They are TFs of interest for the discussion; ABF3 is the
  direct transcriptional output of SnRK2. They had been left out only by the cap of 2 TFs:
  neither shares a pattern or a tomato ortholog with the other nodes.
  - **ABF3** changes in both species, but its pattern is less rare in Arabidopsis
    (frequency 7).
  - **HB7** does not change significantly in Arabidopsis (log2FC 0.34) and its tomato
    pattern is shared by 6 genes.
- **Excluded: DWA2 and ABA4.**
  - **DWA2** had the smallest change of the 15 (log2FC 0.54 in Arabidopsis and 0.19, not
    significant, in tomato) and barely varies in either species: the vst range is 0.5 in
    Arabidopsis and 0.4 in tomato, so its binary patterns were practically noise.
  - **ABA4** is flat in tomato (range 0.3, no significant change), and synthesis is still
    represented by NCED5, which changes in both species.

In `candidates.tsv`, the `auto_selected` column marks the automatic selection and
`excluded_reason` records the manual changes.

### Final nodes

| # | Group | Arabidopsis | AGI | Tomato | Ath pattern (freq) | Sly pattern (freq) | Changes in |
|---|---|---|---|---|---|---|---|
| 1 | TF | MYB44 | AT5G67300 | Solyc02g092930.1 | 0010\|1011 (1) | 1000\|0111 (1) | both |
| 2 | TF | HB6 | AT2G22430 | Solyc05g051460.3 | 0111\|0111 (2) | 1001\|0001 (1) | both |
| 3 | TF | HB7 | AT5G46880 | Solyc08g076370.3 | 1100\|1110 (1) | 1000\|1110 (6) | Sly |
| 4 | TF | ABF3 | AT4G34000 | Solyc11g044560.2 | 0001\|0011 (7) | 0010\|1100 (2) | both |
| 5 | Synthesis | NCED5 | AT1G30100 | Solyc07g056570.1 | 0110\|0111 (1) | 0100\|1000 (1) | both |
| 6 | Degradation | CYP707A3 | AT5G45340 | Solyc08g005610.3 | 1010\|1010 (1) | 1011\|0101 (1) | both |
| 7 | Degradation | CYP707A4 | AT3G19270 | Solyc01g108210.3 | 1000\|1011 (1) | 0111\|1101 (1) | Sly |
| 8 | Transport | AIT1 | AT1G69850 | Solyc04g005790.3 | 1001\|1110 (1) | 0011\|1001 (1) | both |
| 9 | Transport | ABCG31 | AT2G29940 | Solyc02g055530.3 | 1000\|1000 (1) | 0010\|1101 (1) | Sly |
| 10 | Receptor | RCAR1 | AT1G01360 | Solyc01g095700.3 | 0001\|0111 (1) | 1100\|1011 (1) | both |
| 11 | Receptor | PYL1 | AT5G46790 | Solyc06g061180.1 | 0110\|1110 (1) | 1101\|1110 (1) | both |
| 12 | Clade A PP2C | AHG1 | AT5G51760 | Solyc03g006960.3 | 0010\|0001 (2) | 0111\|0101 (1) | both |
| 13 | Clade A PP2C | HAI3 | AT2G29380 | Solyc06g076400.3 | 0100\|0001 (2) | 0111\|0111 (2) | both |
| 14 | Kinase | SnRK2.8 | AT1G78290 | Solyc04g074500.3 | 0001\|0110 (1) | 1011\|1011 (1) | Ath |
| 15 | Kinase | CPK1 | AT5G04870 | Solyc04g009800.3 | 1100\|1100 (1) | 1100\|1000 (2) | both |

The 15 patterns are all distinct within each species. All have frequency 1 or 2, except
ABF3 in Arabidopsis (7) and HB7 in tomato (6).

Pattern = control|drought.
- Arabidopsis: CTRL_0, 7200, 15840, 20160 | CTRL_0, WD 7200, 15840, 20160.
- Tomato: d4, d5, d8, d10 | d4, d5, d8, d10.

The details of each node (significant time points, maximum |log2FC|, number of orthologs)
are in `selected_nodes.tsv`. The 74 candidates, each with its exclusion reason, are in
`candidates.tsv`.

## What is left out and why

- **ABF1, ABF2, ABF4, OST1, ABI2, HAB1 and HAB2** are not included because their patterns
  are very frequent, especially in Arabidopsis (freq 31) or in tomato (OST1, freq 14). They
  are canonical nodes, but in these data their binary dynamics cannot be told apart from
  those of many other genes.
- **PYR1** shares a tomato ortholog with PYL1 (`Solyc06g061180.1`), and **HB5** with HB6
  (`Solyc05g051460.3`).
- **SnRK2.2 and SnRK2.3** are left out by the cap of 2 kinases; SnRK2.8 and CPK1 have rarer
  patterns.
- **DWA2 and ABA4** were picked by the automatic selection, but removed by hand (see above).

## Caveats about network nodes

- **AHG1 and HAI3:** in Arabidopsis they are barely expressed in control (vst between 4.2
  and 5.4). Their control pattern comes from small differences, close to noise. Their
  drought pattern does reflect a real increase at the end (log2FC 9 and 10).
- **ABF3:** its Arabidopsis pattern is shared by 7 genes of the full list, though by none of
  the other 14 nodes. If it shows up as a regulator in the model, those genes could play the
  same role.
- **HB7:** does not change significantly in Arabidopsis, and its tomato pattern is shared by
  6 genes.
- **Nodes that change significantly in one species only:** SnRK2.8 (Arabidopsis only),
  CYP707A4 and ABCG31 (tomato only). In the other species, their binary pattern reflects
  non-significant variation.

## Figure

`fig_boolean_nodes.pdf/.png` (script `plot_nodes.py`) shows, for each node:
- the functional group;
- the Arabidopsis and tomato genes;
- the control|drought SSD binary pattern of each species (dark cell = 1);
- the maximum log2FC (WD vs CTRL) of each species.

For the maximum log2FC, the time point with the largest |log2FC| among the significant ones
is taken. If none is significant, the largest |log2FC| is used and the bar is drawn as an
outline only. The axis is truncated at ±6 and out-of-range values are annotated (MYB44 in
tomato: 24.6; AHG1 and HAI3 in Arabidopsis: 9.1 and 10.7). The values and the time point of
each are in `fig_boolean_nodes_lfc.tsv`.

## Files for the Boolean network

`extract_binarized.py` saves, for the 15 nodes, the SSD binarized matrix and the mean
normalized count (vst, which is what was binarized), for control and drought in each
species:

| Type | Folder | File |
|---|---|---|
| Binary | `{Ath,Sly}_Drought_Leaves/Binarization/` | `ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_selected_nodes.tsv` |
| Continuous | `{Ath,Sly}_Drought_Leaves/Selected_genes/` | `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv` |

- **Format:** the same as the source tables.
  - Binary: one node per column and one time point per row (4 rows, no row names).
  - Continuous: one node per row (`ID` column) and the 4 time points in columns.
- **Names:** `ID (symbol)`, for example `AT1G01360 (RCAR1)` or `Solyc01g095700.3 (RCAR1)`.
  In tomato the symbol is that of the orthologous Arabidopsis gene.
- **`*_symbol.tsv` version:** each file also comes with only the symbol (`RCAR1`,
  `NCED5`, …), so the nodes have the same names in both species.
- **Node order:** the same in the 16 files and in both species, so position *i* is the same
  ortholog pair:
  RCAR1, NCED5, AIT1, CYP707A3, PYL1, MYB44, SnRK2.8, ABCG31, CYP707A4, HB6, CPK1, AHG1, HAI3, HB7, ABF3.
  This is the order of `selected_nodes.tsv`, not that of the table above, which is grouped
  by function.
- **Check:** the script verifies that each binary column reproduces the pattern used in the
  selection.

## Reproduce

```bash
python3 select_nodes.py        # candidates.tsv and selected_nodes.tsv
python3 plot_nodes.py          # figure and fig_boolean_nodes_lfc.tsv
python3 extract_binarized.py   # binary and continuous matrices of the 15 nodes
```

- The parameters are at the top of `select_nodes.py`: `MAX_NODES`, `MAX_PER_GROUP`,
  `MANUAL_INCLUDE` and `MANUAL_EXCLUDE`.
- The folder is in `Ath_Drought_Leaves/`; the scripts locate the other project folders from
  their own position.
