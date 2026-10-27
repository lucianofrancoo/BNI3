# Tomato orthologs of the ABA genes: mature leaves, drought (PRJNA1125861)

The same analysis as `../../Ath_Drought_Leaves/Selected_genes/`, applied to the tomato
(SL3.0) orthologs of the Arabidopsis ABA genes.

## Genes

- `selected_genes.R` builds the list from `../../TAIR_To_ITAG/Table_S_ABA_orthologs.tsv`,
  with the orthologs from OrthoFinder or Ensembl Compara.
- `aba_orthologs_sly.tsv` has one tomato gene per row, with the Arabidopsis genes it is an
  ortholog of (`ath_genes`, `ath_loci`) and its functional group. That is 116 genes.
- 97 pass the DESeq2 filter. The other 19 are barely expressed in leaves: they have between
  0 and 46 reads summed over the 21 samples. The script prints the list.

## `selected_genes.R`

1. **Normalized count** (vst) averaged per condition, for these genes only, from `../DGE/`:
   - `ConteoNormalizado_mean_aba_genes.tsv`: all conditions, with `ath_genes` and `group`.
   - `ConteoNormalizado_mean_Control_aba_genes.tsv`: `ID` and `CTRL_d4` to `CTRL_d10`.
   - `ConteoNormalizado_mean_Drought_aba_genes.tsv`: `ID` and `WD_d4` to `WD_d10`.
   - The control and drought tables have only the ID and the values; they are the ones
     that were binarized (`../Binarization/`).
   - This experiment has no time 0.
2. **Dot plot** of the WD vs CTRL log2FC at d4, d5, d8 and d10, with significance
   (FDR < 0.05): `dotplot_aba_genes.png/.svg`. The values are in
   `dotplot_aba_genes_table.tsv`.
   - Each row is labeled `Solyc… (Arabidopsis ortholog)`.
   - Within each group, the rows follow the order of the Arabidopsis list, so the copies of
     a family stay together.
   - Point size scales with the maximum |log2FC|, which is 24.5 (the MYB44 ortholog at
     d10). That is why most points look small.

## Boolean network nodes

These files are written by
`../../Ath_Drought_Leaves/Boolean_node_selection/extract_binarized.py` (see the README in
that folder):
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv`: the normalized
  count of the 15 nodes only, with `ID` as `Solyc (symbol)`. The symbol is that of the
  orthologous Arabidopsis gene.
- `..._selected_nodes_symbol.tsv`: the same table, with only the symbol as `ID`, as in
  Arabidopsis.

The rows are in the same order as in the Arabidopsis files.
