# Expression binarization: ABA genes, tomato (mature leaves, drought; PRJNA1125861)

## Input and output files

- **What was binarized:** the mean normalized count (vst) in
  `../Selected_genes/ConteoNormalizado_mean_{Control,Drought}_aba_genes.tsv`. That is
  97 genes, with 4 time points per condition: d4, d5, d8 and d10.
- **Methods:** two, SSD and WCSS (BNI3 binarization, with the behavior reviewer). Each
  condition was binarized separately.
  - `*_binarized_{SSD,WCSS}.tsv`: the binary matrix, one gene per column and one time point
    per row (no row names).
  - `*_binarized_{SSD,WCSS}_pattern.tsv`: each gene's 4-bit pattern and its frequency,
    that is, how many genes share that pattern.
  - `*_reviewer_summary.tsv`: the summary per method. It includes how many genes have a
    unique pattern, how many distinct patterns there are and the maximum frequency.

## Boolean network nodes

**SSD** was chosen for the network: the method with more unique-pattern genes, with SSD on
ties. The node selection is documented in
`../../Ath_Drought_Leaves/Boolean_node_selection/README.md`.

`../../Ath_Drought_Leaves/Boolean_node_selection/extract_binarized.py` writes:
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_selected_nodes.tsv`:
  the SSD binary matrix of the 15 nodes only, with `Solyc (symbol)` columns.
- `..._selected_nodes_symbol.tsv`: the same matrix, with only the symbol as column name.
  The symbol is that of the orthologous Arabidopsis gene, so the nodes have the same names
  in both species.

The columns follow the order of `selected_nodes.tsv`, which is the same in both species and
the same as in the continuous tables `../Selected_genes/*_selected_nodes*.tsv`.
