#!/usr/bin/env Rscript

# Agrega los abundance.h5 de kallisto en una matriz de conteos.
# Uso: Rscript extraerCounts.r [directorio_Counts]

library(tximport)

args <- commandArgs(trailingOnly = TRUE)
conteo_root <- if (length(args) >= 1) args[1] else getwd()

setwd(conteo_root)

meta_file   <- file.path(conteo_root, "..", "raw_data", "samples_metadata.tsv")
output_file <- file.path(conteo_root, "Conteo_genes_Sly_Drought_Leaves.tsv")
tpm_file    <- file.path(conteo_root, "TPM_genes_Sly_Drought_Leaves.tsv")
coldata_out <- file.path(conteo_root, "coldata.tsv")

if (!file.exists(meta_file)) stop("No se encuentra ", meta_file)
meta <- read.delim(meta_file, stringsAsFactors = FALSE)

# Ordenar por dia y tratamiento (control primero), no alfabeticamente por SRR
meta$Day_num <- as.integer(sub("^d", "", meta$Day))
meta <- meta[order(meta$Day_num, meta$Treatment != "Control", meta$Run), ]

# Quedarse con las muestras que efectivamente tienen abundance.h5
meta$h5 <- file.path(meta$Run, "abundance.h5")
falta <- !file.exists(meta$h5)
if (any(falta)) {
  cat("AVISO: sin abundance.h5 (se omiten):",
      paste(meta$Run[falta], collapse = ", "), "\n")
}
meta <- meta[!falta, ]

if (nrow(meta) == 0) stop("No se encontraron archivos abundance.h5 en ", conteo_root)

# Nombre de columna: SRR + grupo, para que la matriz se lea sola
sample_names <- paste0(meta$Run, "_", meta$Group)

h5_rutas <- meta$h5
names(h5_rutas) <- sample_names

cat("Encontrados", length(h5_rutas), "archivos de cuantificacion.\n")
cat("Muestras:", paste(sample_names, collapse = ", "), "\n\n")

# El indice SL3.0 tiene 34658 targets = 34658 genes (1 transcrito por gen),
# asi que txOut = TRUE ya entrega nivel de gen y no hace falta tx2gene.
txi.kallisto <- tximport(h5_rutas, type = "kallisto", txOut = TRUE)

# Los targets del FASTA son transcritos (Solyc06g043170.3.1) y el gen es
# Solyc06g043170.3. Se quita el sufijo de transcrito para que los IDs empaten
# con las matrices previas (ej. librerias_Luciano_Tomate_Sequia_Raices/DGE).
# Verificado: la conversion es exacta en los 34658 targets.
rownames(txi.kallisto$counts)    <- sub("\\.[0-9]+$", "", rownames(txi.kallisto$counts))
rownames(txi.kallisto$abundance) <- sub("\\.[0-9]+$", "", rownames(txi.kallisto$abundance))

escribir <- function(mat, archivo, sufijo) {
  IDs <- rownames(mat)
  out <- cbind(gene_id = IDs, mat)
  colnames(out) <- c("gene_id", paste0(colnames(mat), sufijo))
  write.table(out, file = archivo, sep = "\t", quote = FALSE, row.names = FALSE)
}

escribir(round(txi.kallisto$counts), output_file, "_count")
escribir(txi.kallisto$abundance, tpm_file, "_tpm")

# Tabla de diseno lista para DESeq2 / edgeR
coldata <- data.frame(
  sample    = sample_names,
  run       = meta$Run,
  day       = factor(meta$Day, levels = paste0("d", sort(unique(meta$Day_num)))),
  treatment = factor(ifelse(meta$Treatment == "Control", "CTRL", "WD"),
                     levels = c("CTRL", "WD")),
  group     = factor(meta$Group),
  stringsAsFactors = FALSE
)
write.table(coldata, file = coldata_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Conteos guardados en:", output_file, "\n")
cat("TPM guardados en:", tpm_file, "\n")
cat("Tabla de diseno en:", coldata_out, "\n")
cat("Dimensiones:", nrow(txi.kallisto$counts), "genes x", ncol(txi.kallisto$counts), "muestras\n")
