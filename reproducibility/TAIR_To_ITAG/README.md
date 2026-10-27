# Ortólogos en tomate de genes ABA de Arabidopsis

Análisis realizado el 23 y 24 de septiembre de 2026.

## Objetivo

Encontrar, para una lista de 95 genes de Arabidopsis relacionados con ABA (síntesis,
degradación, transporte, percepción, señalización, factores de transcripción y genes
blanco), **todos** sus ortólogos en tomate, con IDs Solyc de la anotación ITAG sobre el
ensamblaje **SL3.0**.

Se buscan todos los ortólogos, no solo uno por gen, porque la mayoría de estas familias
se expandieron por separado en cada linaje después de que Arabidopsis y tomate se
separaron (hace unos 110 millones de años). Por eso la relación esperable suele ser de
uno a muchos o de muchos a muchos, no de uno a uno.

## Resultado en breve

- **86 de 95 genes de Arabidopsis tienen al menos un ortólogo en tomate.** Son 153 pares
  Arabidopsis–tomate, con 116 genes de tomate distintos.
- Cuántos genes cubre cada método:
  - OrthoFinder: 72.
  - Ensembl Compara: 71.
  - Ambos métodos: 57.
  - Solo OrthoFinder: 15.
  - Solo Compara: 14.
- **Sin ortólogo (9):**
  - NFYC2, BG1, CYP707A2, ABCG30, TMAC2 y COR15.
  - DREB2A, PYL10 y RCAR3: solo tenían apoyo por BLAST, que no se usa como criterio.
- La tabla suplementaria es `Table_S_ABA_orthologs.tsv`.
- Las figuras son `fig_orthologs_summary.pdf` y `fig_orthologs_detail.pdf`.

**Definición usada:** un gen de tomate es ortólogo de un gen de Arabidopsis si lo infiere
**OrthoFinder o Ensembl Compara**. Los dos métodos se basan en árboles génicos.

---

## 1. Datos de entrada

### Provistos por el usuario (origen a completar)

| Archivo | Contenido | md5 |
|---|---|---|
| `../NCBI_To_ITAG/Ath/Ath_ABA_genes.tsv` | Lista curada de 95 genes ABA: `Symbol`, `gene_id` (AGI), `hormone`, `group` | — |
| `TAIR12_cdna.fa` | cDNA de Arabidopsis, TAIR12 (54.906 transcritos) | 8d6cc1c8ab263af0833353eca3228b5c |
| `Orthofinder/TAIR12_pep.fa.gz` | Proteínas de Arabidopsis, TAIR12 (36.693) | 4c9889f5e046b53ad9d7e4698de0500a |
| `Sly3.0_cdna.fa` | cDNA de tomate SL3.0 en formato Ensembl (35.825), descomprimido del `Sly3.0_cdna.bgz` original | — |
| `Orthofinder/SL3.0.fa.bgz` | Proteínas de tomate SL3.0 en formato Ensembl (34.429, una por gen) | 45a34c22821cc311edd69d2a704f04ba |

Los encabezados de los archivos de tomate tienen el formato de Ensembl Plants, con IDs
Solyc de la anotación ITAG3 sobre SL3.0. Conviene anotar de qué release de Ensembl
salieron.

### Descargados (24 de septiembre de 2026)

**Ensembl Compara**, homologías de Ensembl Plants **release 62**:
`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/tsv/ensembl-compara/homologies/arabidopsis_thaliana/Compara.115.protein_default.homologies.tsv.gz`

- Se guardaron solo las filas Arabidopsis → *Solanum lycopersicum*, 28.220 pares, en
  `compara62_ath_sly_homologies.tsv.gz` (md5 8756d4c39f115d87c1199cc708338837).
- Se usa la release 62 porque es **la última de Ensembl Plants con tomate en SL3.0**.
  Desde la 63 el tomate está en SL4.0 (`solanum_lycopersicum_gca000188115v5cm`), con
  otros IDs.

**Proteomas de 8 especies más para OrthoFinder**, de Ensembl Plants release 62:
`https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/<especie>/pep/`.
Están en `Orthofinder/ensembl62_pep/`, cada uno con su archivo `CHECKSUMS` de Ensembl, y
todas las sumas se verificaron.

| Especie | Archivo | Genes (una isoforma) |
|---|---|---|
| *Brassica rapa* | Brassica_rapa.Brapa_1.0.pep.all.fa.gz | 41.018 |
| *Solanum tuberosum* | Solanum_tuberosum.SolTub_3.0.pep.all.fa.gz | 39.021 |
| *Capsicum annuum* | Capsicum_annuum.ASM51225v2.pep.all.fa.gz | 35.845 |
| *Nicotiana attenuata* | Nicotiana_attenuata.NIATTr2.pep.all.fa.gz | 33.320 |
| *Coffea canephora* | Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa.gz | 25.574 |
| *Vitis vinifera* | Vitis_vinifera.ASM3070453v1.pep.all.fa.gz | 37.480 |
| *Populus trichocarpa* | Populus_trichocarpa.Pop_tri_v4.pep.all.fa.gz | 34.699 |
| *Oryza sativa* | Oryza_sativa.IRGSP-1.0.pep.all.fa.gz | 35.804 |

Las descripciones de los genes de tomate salen del GFF de Ensembl Plants 62
(`../NCBI_To_ITAG/Ath/Solanum_lycopersicum.SL3.0.62.gff3.gz`, documentado en el README
de esa carpeta).

---

## 2. Qué se hizo y por qué

### 2.1 Conceptos

- **Homólogos:** genes con un ancestro común.
- **Ortólogos:** homólogos que se separaron por **especiación**.
- **Parálogos:** homólogos que se separaron por **duplicación**.

Tener un buen alineamiento (e-value bajo) prueba que dos genes son homólogos, pero no
distingue ortólogos de parálogos. Para eso hay que reconstruir la historia de la familia:
un árbol génico comparado con el árbol de especies.

### 2.2 Primer intento: BLAST con mejor hit recíproco (`blast_ath_aba_to_sly.py`)

- **Qué:**
  - Se extrajeron los cDNA de los 95 genes (`Ath_ABA_genes_cdna.fa`, 120 transcritos,
    todas las isoformas, con el símbolo en el encabezado).
  - Búsqueda **tblastx** (BLAST+ 2.12.0) contra el cDNA de tomate, con e-value ≤ 1e-10.
  - Vuelta: los 5 mejores genes de tomate de cada gen se buscaron contra todo TAIR12.
  - Todo se resumió a nivel de gen, quedándose con el mejor HSP de cada par.
- **Por qué tblastx y no blastn:** entre especies tan divergentes el ADN está saturado de
  cambios en las terceras posiciones de los codones. Comparar las traducciones conserva
  mucha más señal.
- **Mejor hit recíproco (RBH):** el mejor gen de tomate del gen A tiene como mejor gen de
  Arabidopsis a A. Da 38 pares.
- **Limitación:** el RBH elige un solo gen por cada lado. Cuando una familia se duplicó en
  un linaje (SnRK2, PYL, PP2C clado A, ABF), las demás copias, que son ortólogas igual,
  quedan afuera. Si una copia se perdió, el RBH puede tomar un parálogo por ortólogo.
- **Uso final:** **solo como validación.** 31 de los 38 pares RBH (82 %) también los dan
  OrthoFinder o Compara. Los otros 7 son un gen sin otro apoyo (DREB2A) y seis pares que
  los métodos basados en árboles rechazan, probablemente parálogos.

### 2.3 DAVID (descartado)

- La herramienta *DAVID Ortholog* no calcula ortología. Convierte los IDs con tablas de
  **OMA** o **Ensembl Compara**.
- El resultado (un archivo ya eliminado de la carpeta) tenía 53 genes con **un solo** gen de tomate
  cada uno, y omitía genes centrales como ABI1, ABI2, RCAR1 y PYL1.
- Se descartó porque es un subconjunto reducido de Compara, que se usa directamente.

### 2.4 Ensembl Compara (release 62)

- Es el sistema de genómica comparativa de Ensembl.
- Arma árboles génicos para unas 100 especies de plantas y los compara con el árbol de
  especies, para clasificar cada par como ortólogo (`one2one`, `one2many`,
  `many2many`) o parálogo.
- Se tomaron todos los ortólogos Arabidopsis → tomate de los 95 genes: 71 genes con
  ortólogo, 106 pares.

### 2.5 OrthoFinder (`Orthofinder/`)

Se corrió OrthoFinder porque permite controlar la anotación (SL3.0), las especies y las
versiones, y porque es el método de ortología filogenética más usado hoy.

1. **Proteomas** (`prepare_proteomes.py` → `primary/`):
   - Una proteína por gen, la isoforma más larga, como recomienda OrthoFinder, para que
     las isoformas no aparezcan como falsas duplicaciones.
   - Los encabezados son el ID de gen.
2. **Especies:** 10 en total. Con solo dos especies no se pueden enraizar los árboles ni
   ubicar las duplicaciones en un linaje.
   - *Arabidopsis* + *Brassica rapa*: separan duplicaciones propias de las Brassicaceae.
   - Tomate + papa, pimiento y *Nicotiana*: separan las duplicaciones de las solanáceas,
     incluida su triplicación genómica.
   - Café: otra astérida.
   - Vid y álamo: rósidas y afines.
   - Arroz: grupo externo.
3. **Corrida** (`run_orthofinder.sh`):
   - OrthoFinder 3.1.5 con DIAMOND 2.2.8 (todos contra todos).
   - MCL 22.282 para formar los orthogroups.
   - Alineamientos múltiples con FAMSA 2.4.1 y árboles con FastTree 2.2.0 (`-M msa`).
   - 32 hilos para las búsquedas y 16 para el análisis.
4. **Árbol de especies:**
   - El árbol que infirió OrthoFinder agrupaba la vid con las astéridas y dejaba afuera a
     las Brassicaceae. Es un artefacto conocido (atracción de ramas largas), porque las
     Brassicaceae evolucionan rápido y la vid lento.
   - Se relanzó con el árbol aceptado (`species_tree_rooted.nwk`, opción `-s`),
     reutilizando las búsquedas de DIAMOND (`-b`):
     ```
     (Oryza,((Vitis,(Populus,(Arabidopsis,Brassica))),(Coffea,(Nicotiana,(Capsicum,(S.tuberosum,S.lycopersicum))))))
     ```
   - Para estos 95 genes, las dos corridas dan **los mismos 125 pares**. Se reporta la
     corrida con el árbol correcto.
5. **Resultado:** 72 genes con ortólogo y 125 pares, de
   `Results_sptree/Orthologues/Orthologues_Arabidopsis_thaliana/Arabidopsis_thaliana__v__Solanum_lycopersicum.tsv`.

### 2.6 Combinación (`combine_orthologs.py`)

- Une en `Ath_ABA_orthologs_combined.tsv` los pares de BLAST, Compara y OrthoFinder, con
  la información de cada fuente. Las columnas de BLAST quedan solo como referencia.
- La tabla final (`Ath_ABA_orthologs_final.tsv`) toma como ortólogos los pares de
  **OrthoFinder o Compara**.
- **Por qué la unión y no la intersección:** cada método falla en familias distintas.
  - Compara no asigna ortólogos a ABI1, ABI2, AHG1, AHG3, RCAR1, PYL3, PYL5 ni PYL7, que
    OrthoFinder sí resuelve.
  - OrthoFinder deja sin ortólogo a NCED3, NCED9, AAO3, ABF1 y ABF4, que Compara sí
    asigna.
  - La intersección perdería genes centrales de la vía. La unión cubre 86 de 95.

### 2.7 Figuras y tabla suplementaria

- `plot_orthologs.py` genera las dos figuras:
  - `fig_orthologs_summary`:
    - Panel A: Venn de genes de Arabidopsis con ortólogo según cada método.
    - Panel B: cobertura por grupo funcional, con la unión y cada método por separado.
  - `fig_orthologs_detail`: cada par con el método que lo apoya.
- `make_supp_table.py` genera `Table_S_ABA_orthologs.tsv`, con una fila por par y una
  fila por cada gen sin ortólogo. Columnas:
  - `functional_group`, `arabidopsis_gene`, `arabidopsis_locus`, `tomato_gene`.
  - `orthofinder`, `ensembl_compara` (Yes/No).
  - `orthofinder_relationship`: dirección Arabidopsis:tomato.
  - `compara_relationship`: tipo según Compara.
  - `orthofinder_orthogroup`.
  - `tomato_description`: la descripción ITAG existe solo para 39 de los 153 pares.

---

## 3. Casos a revisar antes de interpretar

- **ABI1 y ABI2 → `Solyc07g040990.3`, solo según OrthoFinder.** El mejor hit de BLAST era
  `Solyc08g062650.3`. Conviene mirar el árbol del orthogroup.
- **NCED3 y NCED9 → `Solyc07g056570.1` (SlNCED1, *notabilis*), solo según Compara.**
  OrthoFinder asigna ese gen a NCED2 y NCED5, e infiere que el linaje de tomate perdió la
  copia de tipo NCED3. Esa parte del árbol (FastTree) está poco resuelta.
- **ABCG40 (9 genes de tomate) y ABCG31 (6):** son familias muy expandidas y la lista
  seguramente incluye copias con funciones distintas.
- **RD29A → `Solyc01g009660.2`, solo según OrthoFinder.** RD29A se considera casi
  exclusivo de las Brassicaceae.

Los árboles de cada familia están en
`Results_sptree/Resolved_Gene_Trees/Resolved_Gene_Trees.txt` (formato Newick; se pueden
visualizar, por ejemplo, con iTOL).

---

## 4. Archivos

| Archivo | Qué es |
|---|---|
| `blast_ath_aba_to_sly.py` | Búsqueda tblastx de ida y vuelta (RBH) |
| `combine_orthologs.py` | Combina las fuentes y genera la tabla final |
| `plot_orthologs.py` | Figuras |
| `make_supp_table.py` | Tabla suplementaria |
| `Orthofinder/prepare_proteomes.py` | Proteomas con una isoforma por gen |
| `Orthofinder/run_orthofinder.sh` | Comandos de OrthoFinder |
| `Orthofinder/species_tree_rooted.nwk` | Árbol de especies de referencia |
| `Ath_ABA_genes_cdna.fa` | cDNA de los 95 genes, con el símbolo en el encabezado |
| `blast_fwd_ath_aba_vs_sly.tsv`, `blast_rev_sly_vs_tair12.tsv` | Hits tblastx crudos |
| `compara62_ath_sly_homologies.tsv.gz` | Homologías Arabidopsis–tomate de Compara 62 |
| `Ath_ABA_orthologs_combined.tsv` | Todas las fuentes, par por par |
| `Ath_ABA_orthologs_final.tsv` | Ortólogos finales (OrthoFinder o Compara) |
| `Table_S_ABA_orthologs.tsv` | Tabla suplementaria |
| `fig_orthologs_summary.pdf/.png`, `fig_orthologs_detail.pdf/.png` | Figuras |
| `Orthofinder/results/Results_Sep24/WorkingDirectory/OrthoFinder/Results_sptree/` | Resultados finales de OrthoFinder |

## 5. Cómo reproducir

```bash
# 1. BLAST (requiere BLAST+ en el PATH)
python3 blast_ath_aba_to_sly.py
# 2. OrthoFinder (entorno conda "orthofinder"; unos 35 min con 32 hilos)
cd Orthofinder && bash run_orthofinder.sh && cd ..
# 3. Combinación, figuras y tabla suplementaria
python3 combine_orthologs.py
python3 plot_orthologs.py
python3 make_supp_table.py
```

Los scripts leen la lista de genes (`Ath_ABA_genes.tsv`) desde `../NCBI_To_ITAG/Ath/`, y
`make_supp_table.py` también lee de ahí el GFF de tomate.

---

## 6. Uso posterior de estos ortólogos

`Table_S_ABA_orthologs.tsv` es la entrada de dos análisis:

1. **Expresión de los ortólogos de tomate en sequía**
   (`../Sly_Drought_Leaves/Selected_genes/selected_genes.R`).
   - Es el mismo análisis que se hizo en Arabidopsis
     (`../Ath_Drought_Leaves/Selected_genes/`), aplicado a los 116 genes de tomate con
     ortólogo.
   - De esos 116, 97 tienen datos de DGE. Los otros 19 casi no se expresan en hoja y el
     filtro de DESeq2 los descartó.
2. **Selección de nodos para redes booleanas en las dos especies**
   (`../Ath_Drought_Leaves/Boolean_node_selection/`).
   - Se quedan los genes de Arabidopsis con ortólogo de tomate, que cambian en al menos una
     especie y que tienen patrones binarios poco frecuentes. Se elige un ortólogo
     representante por gen.
   - La lista final tiene 15 pares: MYB44, HB6, HB7, ABF3, NCED5, CYP707A3, CYP707A4, AIT1,
     ABCG31, RCAR1, PYL1, AHG1, HAI3, SnRK2.8 y CPK1.

## 7. Texto sugerido para Métodos (inglés)

> Tomato orthologs of 95 Arabidopsis ABA-related genes were identified with two gene-tree
> based methods. First, OrthoFinder v3.1.5 (Emms & Kelly, 2019) was run on the proteomes
> of ten species (longest isoform per gene): *Arabidopsis thaliana* (TAIR12), *Solanum
> lycopersicum* (SL3.0, ITAG3), and, from Ensembl Plants release 62, *S. tuberosum*,
> *Capsicum annuum*, *Nicotiana attenuata*, *Coffea canephora*, *Vitis vinifera*,
> *Populus trichocarpa*, *Brassica rapa* and *Oryza sativa*. We used DIAMOND searches,
> MSA-based gene trees (FAMSA, FastTree) and a user-supplied rooted species tree. Second,
> Arabidopsis–tomato orthologs were retrieved from Ensembl Compara (Ensembl Plants
> release 62, the last release based on SL3.0). A tomato gene was considered an ortholog
> when inferred by either method, retaining one-to-many and many-to-many relationships.
> As an independent check, 31 of 38 reciprocal best tBLASTX hits were also recovered by
> at least one of these methods.

Citas:
- OrthoFinder: Emms & Kelly (2019) *Genome Biology* 20:238; Emms et al. (2025) bioRxiv
  doi:10.1101/2025.07.15.664860.
- Ensembl Plants / Compara: citar la publicación de Ensembl Genomes vigente.
