# Selected ABA genes: Arabidopsis, rosette leaves, drought

## Genes

`aba_genes.tsv` has the 95 ABA pathway genes taken from Song et al. (2016), with these
columns:
- `geneName`: the symbol.
- `gene_id`: the AGI.
- `hormone`.
- `group`: the functional group.

89 of them pass the DESeq2 filter. The other 6 have no data because of low counts: NCED6,
PYL10, PYL11, PYL12, PYL13 and PYL2.

## `selected_genes.R`

1. **Normalized count** (vst) averaged per condition, for the listed genes only, from
   `../DGE/`:
   - `ConteoNormalizado_mean_aba_genes.tsv`: all conditions, with `geneName` and `group`.
   - `ConteoNormalizado_mean_Control_aba_genes.tsv`: `ID`, `CTRL_0`, `CTRL_7200`,
     `CTRL_15840` and `CTRL_20160`.
   - `ConteoNormalizado_mean_Drought_aba_genes.tsv`: `ID`, `CTRL_0`, `WD_7200`,
     `WD_15840` and `WD_20160`.
   - The control and drought tables have only the ID and the values; they are the ones
     that were binarized (`../Binarization/`).
2. **Dot plot** of the WD vs CTRL log2FC at each time point, with significance
   (FDR < 0.05): `dotplot_aba_genes.png/.svg`. The values are in
   `dotplot_aba_genes_table.tsv`.

## Boolean network nodes

These files are written by `../Boolean_node_selection/extract_binarized.py` (see the README
in that folder):
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv`: the same
  normalized count, for the 15 chosen nodes only, with `ID` as `AGI (symbol)`.
- `..._selected_nodes_symbol.tsv`: the same table, with only the symbol as `ID`.

The rows follow the order of `../Boolean_node_selection/selected_nodes.tsv`, which is the
same as in the binary matrices and in the tomato files.
