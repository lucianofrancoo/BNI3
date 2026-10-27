# Ortólogos en tomate de los genes ABA: hoja madura, sequía (PRJNA1125861)

Es el mismo análisis que `../../Ath_Drought_Leaves/Selected_genes/`, aplicado a los
ortólogos de tomate (SL3.0) de los genes ABA de Arabidopsis.

## Genes

- `selected_genes.R` arma la lista a partir de
  `../../TAIR_To_ITAG/Table_S_ABA_orthologs.tsv`, con los ortólogos por OrthoFinder o
  Ensembl Compara.
- `aba_orthologs_sly.tsv` tiene un gen de tomate por fila, con los genes de Arabidopsis de
  los que es ortólogo (`ath_genes`, `ath_loci`) y su grupo funcional. Son 116 genes.
- 97 pasan el filtro de DESeq2. Los otros 19 casi no se expresan en hoja: tienen entre 0 y
  46 lecturas sumando las 21 muestras. El script imprime la lista.

## `selected_genes.R`

1. **Conteo normalizado** (vst) promediado por condición, solo de estos genes, a partir de
   `../DGE/`:
   - `ConteoNormalizado_mean_aba_genes.tsv`: todas las condiciones, con `ath_genes` y
     `group`.
   - `ConteoNormalizado_mean_Control_aba_genes.tsv`: `ID` y `CTRL_d4` a `CTRL_d10`.
   - `ConteoNormalizado_mean_Drought_aba_genes.tsv`: `ID` y `WD_d4` a `WD_d10`.
   - Las tablas de control y sequía tienen solo el ID y los valores; son las que se
     binarizaron (`../Binarization/`).
   - Este experimento no tiene tiempo 0.
2. **Dotplot** del log2FC WD vs CTRL en d4, d5, d8 y d10, con significancia (FDR < 0,05):
   `dotplot_aba_genes.png/.svg`. Los valores están en `dotplot_aba_genes_table.tsv`.
   - Cada fila se rotula `Solyc… (ortólogo de Arabidopsis)`.
   - Dentro de cada grupo, las filas siguen el orden de la lista de Arabidopsis, para que
     las copias de una familia queden juntas.
   - El tamaño de los puntos se escala con el |log2FC| máximo, que es 24,5 (el ortólogo de
     MYB44 en d10). Por eso la mayoría de los puntos se ven chicos.

## Nodos de la red booleana

Estos archivos los genera
`../../Ath_Drought_Leaves/Boolean_node_selection/extract_binarized.py` (ver el README de esa
carpeta):
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv`: el conteo
  normalizado solo de los 15 nodos, con `ID` como `Solyc (símbolo)`. El símbolo es el del
  gen de Arabidopsis ortólogo.
- `..._selected_nodes_symbol.tsv`: la misma tabla, con solo el símbolo como `ID`, igual que
  en Arabidopsis.

Las filas están en el mismo orden que en los archivos de Arabidopsis.
