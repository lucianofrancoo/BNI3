# Selección de nodos para redes booleanas: respuesta ABA a sequía en Arabidopsis y tomate

Análisis del 24 de septiembre de 2026. Script: `select_nodes.py`.

## Objetivo

Elegir una lista corta (≤ 15) de genes de la vía ABA para construir redes booleanas en
Arabidopsis y en tomate, y comparar la respuesta a sequía entre las dos especies. Cada nodo
es un par formado por un gen de Arabidopsis y un ortólogo de tomate.

Los genes de Arabidopsis son los 95 de la vía ABA tomados de Song et al. (2016). Los
ortólogos de tomate están descritos en `../../TAIR_To_ITAG/README.md`.

## Datos

| Qué | Archivo |
|---|---|
| Ortólogos Arabidopsis–tomate | `../../TAIR_To_ITAG/Table_S_ABA_orthologs.tsv` (OrthoFinder o Ensembl Compara) |
| DGE de Arabidopsis (roseta; WD vs CTRL a 7200, 15840 y 20160) | `../Selected_genes/dotplot_aba_genes_table.tsv` |
| DGE de tomate (hoja madura; WD vs CTRL en d4, d5, d8 y d10) | `../../Sly_Drought_Leaves/Selected_genes/dotplot_aba_genes_table.tsv` |
| Patrones binarios | `../../{Ath,Sly}_Drought_Leaves/Binarization/ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_pattern.tsv` |

## Criterios, en orden

1. **F1, ortología:** el gen de Arabidopsis tiene datos de expresión y al menos un ortólogo
   de tomate con datos. Quedan **76 de 95** genes.
2. **F2, cambio:** padj < 0,05 (cualquier |log2FC| > 0) en al menos un tiempo, en al menos
   una especie. En tomate cuenta cualquiera de los ortólogos. Quedan **74**: se van ABI1 y
   AFP2.
3. **Método de binarización: SSD.** Se eligió el método con más genes de patrón único.

   | Tabla | SSD | WCSS |
   |---|---|---|
   | Arabidopsis, control | 3 | 1 |
   | Arabidopsis, sequía | 3 | 3 |
   | Tomate, control | 1 | 1 |
   | Tomate, sequía | 2 | 1 |

   En los empates se usa SSD.
4. **Patrón de 8 bits: control|sequía.** Cada condición se binarizó por separado (4 tiempos
   cada una). Se concatenan porque la red modela las dos condiciones y porque 4 bits dan
   muy poca resolución: en la sequía de Arabidopsis, 51 de 89 genes comparten el patrón
   `0001`. Con 8 bits hay 35 patrones distintos en Arabidopsis (21 únicos) y 51 en tomate
   (36 únicos).
   - La **frecuencia** es el número de genes con ese patrón entre todos los binarizados de
     la especie: 89 en Arabidopsis y 97 en tomate.
5. **Ortólogo representante:** uno por gen de Arabidopsis. Se prefiere, en este orden:
   - el que cambia significativamente;
   - el de patrón menos frecuente;
   - el que dan OrthoFinder y Compara a la vez;
   - el de mayor |log2FC|.
6. **Orden de prioridad:**
   - Primero, la frecuencia más alta de las dos especies, así se prefiere el patrón raro
     en **ambas**.
   - Después, la suma de las dos frecuencias.
   - Por último, que cambie significativamente en las dos especies.
7. **Regla de no co-inclusión:** dentro de cada especie no entran dos nodos con el mismo
   patrón. Si dos genes comparten patrón, la red booleana no puede distinguir cuál causa un
   cambio aguas abajo. Tampoco entran dos genes de Arabidopsis con el mismo ortólogo de
   tomate.
8. **Cobertura de la vía:** como máximo 2 nodos por grupo funcional, en total 15.
   - Sin este tope, el orden por frecuencia llenaba la lista con receptores y genes de
     síntesis y dejaba afuera a las PP2C del clado A, un eslabón central del eje
     PYL → PP2C → SnRK2.

## Resultado: 15 nodos

### Selección automática y curación manual

La selección automática (criterios 1–8) dio 15 nodos. Después se hicieron dos cambios a
mano, registrados en `MANUAL_INCLUDE` y `MANUAL_EXCLUDE` de `select_nodes.py`. El script
verifica que la selección final siga cumpliendo la regla de no co-inclusión.

- **Incluidos: HB7 y ABF3.** Son TF de interés para la discusión; ABF3 es la salida
  transcripcional directa de SnRK2. Solo habían quedado afuera por el tope de 2 TF: ninguno
  comparte patrón con los demás nodos ni ortólogo de tomate.
  - **ABF3** cambia en las dos especies, pero su patrón es menos raro en Arabidopsis
    (frecuencia 7).
  - **HB7** no cambia significativamente en Arabidopsis (log2FC 0,34) y su patrón de
    tomate lo comparten 6 genes.
- **Excluidos: DWA2 y ABA4.**
  - **DWA2** tenía el cambio más chico de los 15 (log2FC 0,54 en Arabidopsis y 0,19, no
    significativo, en tomate) y casi no varía en ninguna especie: el rango de vst es de 0,5
    en Arabidopsis y 0,4 en tomate, así que sus patrones binarios eran prácticamente ruido.
  - **ABA4** está plano en tomate (rango 0,3, sin cambio significativo) y la síntesis sigue
    representada por NCED5, que cambia en las dos especies.

En `candidates.tsv`, la columna `auto_selected` indica la selección automática y
`excluded_reason` marca los cambios manuales.

### Nodos finales

| # | Grupo | Arabidopsis | AGI | Tomate | Patrón Ath (freq) | Patrón Sly (freq) | Cambia en |
|---|---|---|---|---|---|---|---|
| 1 | TF | MYB44 | AT5G67300 | Solyc02g092930.1 | 0010\|1011 (1) | 1000\|0111 (1) | ambas |
| 2 | TF | HB6 | AT2G22430 | Solyc05g051460.3 | 0111\|0111 (2) | 1001\|0001 (1) | ambas |
| 3 | TF | HB7 | AT5G46880 | Solyc08g076370.3 | 1100\|1110 (1) | 1000\|1110 (6) | Sly |
| 4 | TF | ABF3 | AT4G34000 | Solyc11g044560.2 | 0001\|0011 (7) | 0010\|1100 (2) | ambas |
| 5 | Síntesis | NCED5 | AT1G30100 | Solyc07g056570.1 | 0110\|0111 (1) | 0100\|1000 (1) | ambas |
| 6 | Degradación | CYP707A3 | AT5G45340 | Solyc08g005610.3 | 1010\|1010 (1) | 1011\|0101 (1) | ambas |
| 7 | Degradación | CYP707A4 | AT3G19270 | Solyc01g108210.3 | 1000\|1011 (1) | 0111\|1101 (1) | Sly |
| 8 | Transporte | AIT1 | AT1G69850 | Solyc04g005790.3 | 1001\|1110 (1) | 0011\|1001 (1) | ambas |
| 9 | Transporte | ABCG31 | AT2G29940 | Solyc02g055530.3 | 1000\|1000 (1) | 0010\|1101 (1) | Sly |
| 10 | Receptor | RCAR1 | AT1G01360 | Solyc01g095700.3 | 0001\|0111 (1) | 1100\|1011 (1) | ambas |
| 11 | Receptor | PYL1 | AT5G46790 | Solyc06g061180.1 | 0110\|1110 (1) | 1101\|1110 (1) | ambas |
| 12 | PP2C clado A | AHG1 | AT5G51760 | Solyc03g006960.3 | 0010\|0001 (2) | 0111\|0101 (1) | ambas |
| 13 | PP2C clado A | HAI3 | AT2G29380 | Solyc06g076400.3 | 0100\|0001 (2) | 0111\|0111 (2) | ambas |
| 14 | Quinasa | SnRK2.8 | AT1G78290 | Solyc04g074500.3 | 0001\|0110 (1) | 1011\|1011 (1) | Ath |
| 15 | Quinasa | CPK1 | AT5G04870 | Solyc04g009800.3 | 1100\|1100 (1) | 1100\|1000 (2) | ambas |

Los 15 patrones son distintos entre sí en cada especie. Todos tienen frecuencia 1 o 2,
salvo ABF3 en Arabidopsis (7) y HB7 en tomate (6).

Patrón = control|sequía.
- Arabidopsis: CTRL_0, 7200, 15840, 20160 | CTRL_0, WD 7200, 15840, 20160.
- Tomate: d4, d5, d8, d10 | d4, d5, d8, d10.

El detalle de cada nodo (tiempos significativos, |log2FC| máximo, número de ortólogos) está
en `selected_nodes.tsv`. Los 74 candidatos, con el motivo de exclusión de cada uno, están en
`candidates.tsv`.

## Qué queda afuera y por qué

- **ABF1, ABF2, ABF4, OST1, ABI2, HAB1 y HAB2** no entran porque sus patrones son muy
  frecuentes, sobre todo en Arabidopsis (freq 31) o en tomate (OST1, freq 14). Son nodos
  canónicos, pero en estos datos su dinámica binaria es indistinguible de la de muchos
  otros genes.
- **PYR1** comparte ortólogo de tomate con PYL1 (`Solyc06g061180.1`), y **HB5** con HB6
  (`Solyc05g051460.3`).
- **SnRK2.2 y SnRK2.3** quedan fuera por el tope de 2 quinasas; SnRK2.8 y CPK1 tienen
  patrones más raros.
- **DWA2 y ABA4** los eligió la selección automática, pero se sacaron a mano (ver arriba).

## Advertencias sobre nodos de la red

- **AHG1 y HAI3:** en Arabidopsis casi no se expresan en control (vst entre 4,2 y 5,4). Su
  patrón de control sale de diferencias chicas, cercanas al ruido. Su patrón de sequía sí
  refleja un aumento real al final (log2FC 9 y 10).
- **ABF3:** su patrón de Arabidopsis lo comparten 7 genes de la lista completa, aunque
  ninguno de los otros 14 nodos. Si aparece como regulador en el modelo, esos genes podrían
  cumplir el mismo papel.
- **HB7:** no cambia significativamente en Arabidopsis, y su patrón de tomate lo comparten
  6 genes.
- **Nodos que cambian significativamente en una sola especie:** SnRK2.8 (solo en
  Arabidopsis), CYP707A4 y ABCG31 (solo en tomate). En la otra especie, su patrón binario
  refleja variaciones no significativas.

## Figura

`fig_boolean_nodes.pdf/.png` (script `plot_nodes.py`) muestra, para cada nodo:
- el grupo funcional;
- los genes de Arabidopsis y de tomate;
- el patrón binario SSD control|sequía de cada especie (celda oscura = 1);
- el log2FC máximo (WD vs CTRL) de cada especie.

Para el log2FC máximo se toma el tiempo con mayor |log2FC| entre los significativos. Si
ninguno es significativo, se usa el mayor |log2FC| y la barra se dibuja solo con contorno.
El eje está truncado en ±6 y los valores fuera de rango van anotados (MYB44 en tomate: 24,6;
AHG1 y HAI3 en Arabidopsis: 9,1 y 10,7). Los valores y el tiempo de cada uno están en
`fig_boolean_nodes_lfc.tsv`.

## Archivos para la red booleana

`extract_binarized.py` guarda, para los 15 nodos, la matriz binarizada SSD y el conteo
normalizado promedio (vst, que es lo que se binarizó), para control y sequía en cada
especie:

| Tipo | Carpeta | Archivo |
|---|---|---|
| Binaria | `{Ath,Sly}_Drought_Leaves/Binarization/` | `ConteoNormalizado_mean_{Control,Drought}_aba_genes_binarized_SSD_selected_nodes.tsv` |
| Continua | `{Ath,Sly}_Drought_Leaves/Selected_genes/` | `ConteoNormalizado_mean_{Control,Drought}_aba_genes_selected_nodes.tsv` |

- **Formato:** el mismo que las tablas de origen.
  - En la binaria, un nodo por columna y un tiempo por fila (4 filas, sin nombre de fila).
  - En la continua, un nodo por fila (columna `ID`) y los 4 tiempos en columnas.
- **Nombres:** `ID (símbolo)`, por ejemplo `AT1G01360 (RCAR1)` o `Solyc01g095700.3 (RCAR1)`.
  En tomate el símbolo es el del gen de Arabidopsis ortólogo.
- **Versión `*_symbol.tsv`:** de cada archivo hay una con solo el símbolo (`RCAR1`,
  `NCED5`, …). Así los nodos se llaman igual en las dos especies.
- **Orden de los nodos:** es el mismo en los 16 archivos y en las dos especies, así que la
  posición *i* es el mismo par de ortólogos:
  RCAR1, NCED5, AIT1, CYP707A3, PYL1, MYB44, SnRK2.8, ABCG31, CYP707A4, HB6, CPK1, AHG1, HAI3, HB7, ABF3.
  Es el orden de `selected_nodes.tsv`, no el de la tabla de arriba, que está agrupada por
  función.
- **Verificación:** el script comprueba que cada columna binaria reproduzca el patrón usado
  en la selección.

## Cómo reproducir

```bash
python3 select_nodes.py        # candidates.tsv y selected_nodes.tsv
python3 plot_nodes.py          # figura y fig_boolean_nodes_lfc.tsv
python3 extract_binarized.py   # matrices binarias y continuas de los 15 nodos
```

- Los parámetros están al principio de `select_nodes.py`: `MAX_NODES`, `MAX_PER_GROUP`,
  `MANUAL_INCLUDE` y `MANUAL_EXCLUDE`.
- La carpeta está en `Ath_Drought_Leaves/`; los scripts ubican las demás carpetas del
  proyecto a partir de su propia posición.
