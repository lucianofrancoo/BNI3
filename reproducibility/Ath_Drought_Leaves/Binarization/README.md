# Expression binarization: ABA genes, Arabidopsis (rosette leaves, drought)

## Input and output files

- **What was binarized:** the mean normalized count (vst) in
  `../Selected_genes/ConteoNormalizado_mean_{Control,Drought}_aba_genes.tsv`. That is
  89 genes, with 4 time points per condition: 0, 7200, 15840 and 20160. In drought the first
  time point is CTRL_0, the shared time 0.
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
ties. The node selection is documented in `../Boolean_node_selection/README.md`.

`../Boolean_node_selection/extract_binarized.py` writes:
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_selected_nodes.tsv`:
  the SSD binary matrix of the 15 nodes only, with `AGI (symbol)` columns.
- `..._selected_nodes_symbol.tsv`: the same matrix, with only the symbol as column name.

The columns follow the order of `selected_nodes.tsv`, which is the same in both species and
the same as in the continuous tables `../Selected_genes/*_selected_nodes*.tsv`.
