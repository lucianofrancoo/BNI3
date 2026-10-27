# Genes ABA seleccionados: Arabidopsis, hoja de roseta, sequía

## Genes

`aba_genes.tsv` tiene los 95 genes de la vía ABA tomados de Song et al. (2016), con estas
columnas:
- `geneName`: el símbolo.
- `gene_id`: el AGI.
- `hormone`.
- `group`: el grupo funcional.

89 de ellos pasan el filtro de DESeq2. Los otros 6 no tienen datos por bajo conteo: NCED6,
PYL10, PYL11, PYL12, PYL13 y PYL2.

## `selected_genes.R`

1. **Conteo normalizado** (vst) promediado por condición, solo de los genes de la lista, a
   partir de `../DGE/`:
   - `ConteoNormalizado_mean_aba_genes.tsv`: todas las condiciones, con `geneName` y
     `group`.
   - `ConteoNormalizado_mean_Control_aba_genes.tsv`: `ID`, `CTRL_0`, `CTRL_7200`,
     `CTRL_15840` y `CTRL_20160`.
   - `ConteoNormalizado_mean_Drought_aba_genes.tsv`: `ID`, `CTRL_0`, `WD_7200`,
     `WD_15840` y `WD_20160`.
   - Las tablas de control y sequía tienen solo el ID y los valores; son las que se
     binarizaron (`../Binarization/`).
2. **Dotplot** del log2FC WD vs CTRL en cada tiempo, con significancia (FDR < 0,05):
   `dotplot_aba_genes.png/.svg`. Los valores están en `dotplot_aba_genes_table.tsv`.

## Nodos de la red booleana

Estos archivos los genera `../Boolean_node_selection/extract_binarized.py` (ver el README
de esa carpeta):
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv`: el mismo conteo
  normalizado, solo para los 15 nodos elegidos, con `ID` como `AGI (símbolo)`.
- `..._selected_nodes_symbol.tsv`: la misma tabla, con solo el símbolo como `ID`.

Las filas están en el orden de `../Boolean_node_selection/selected_nodes.tsv`, que es el
mismo que el de las matrices binarias y el de los archivos de tomate.
