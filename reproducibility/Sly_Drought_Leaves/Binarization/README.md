# Binarización de la expresión: genes ABA, tomate (hoja madura, sequía; PRJNA1125861)

## Archivos de entrada

Provistos por el usuario; el script de binarización no está en esta carpeta.

- **Qué se binarizó:** el conteo normalizado promedio (vst) de
  `../Selected_genes/ConteoNormalizado_mean_{Control,Drought}_aba_genes.tsv`. Son
  97 genes, con 4 tiempos por condición: d4, d5, d8 y d10.
- **Métodos:** hay dos, SSD y WCSS. Cada condición se binarizó por separado.
  - `*_binarized_{SSD,WCSS}.tsv`: la matriz binaria, con un gen por columna y un tiempo
    por fila (sin nombres de fila).
  - `*_binarized_{SSD,WCSS}_pattern.tsv`: el patrón de 4 bits de cada gen y su frecuencia,
    es decir, cuántos genes comparten ese patrón.
  - `*_reviewer_summary.tsv`: el resumen por método. Incluye cuántos genes tienen un patrón
    único, cuántos patrones distintos hay y la frecuencia máxima.

## Nodos de la red booleana

Para la red se eligió **SSD**, el método con más genes de patrón único, usando SSD en los
empates. La selección de nodos está documentada en
`../../Ath_Drought_Leaves/Boolean_node_selection/README.md`.

`../../Ath_Drought_Leaves/Boolean_node_selection/extract_binarized.py` genera:
- `ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_selected_nodes.tsv`:
  la matriz binaria SSD solo de los 15 nodos, con columnas `Solyc (símbolo)`.
- `..._selected_nodes_symbol.tsv`: la misma matriz, con solo el símbolo como nombre de
  columna. El símbolo es el del gen de Arabidopsis ortólogo, así que los nodos se llaman
  igual en las dos especies.

Las columnas están en el orden de `selected_nodes.tsv`, que es el mismo en las dos especies
y el mismo que el de las tablas continuas `../Selected_genes/*_selected_nodes*.tsv`.
