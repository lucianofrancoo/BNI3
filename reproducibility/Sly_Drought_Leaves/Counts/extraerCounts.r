#!/usr/bin/env Rscript

# Aggregates the kallisto abundance.h5 files into a count matrix.
# Usage: Rscript extraerCounts.r [Counts_directory]

library(tximport)

args <- commandArgs(trailingOnly = TRUE)
conteo_root <- if (length(args) >= 1) args[1] else getwd()

setwd(conteo_root)

meta_file   <- file.path(conteo_root, "..", "raw_data", "samples_metadata.tsv")
output_file <- file.path(conteo_root, "Conteo_genes_Sly_Drought_Leaves.tsv")
tpm_file    <- file.path(conteo_root, "TPM_genes_Sly_Drought_Leaves.tsv")
coldata_out <- file.path(conteo_root, "coldata.tsv")

if (!file.exists(meta_file)) stop("Not found: ", meta_file)
meta <- read.delim(meta_file, stringsAsFactors = FALSE)

# Order by day and treatment (control first), not alphabetically by SRR
meta$Day_num <- as.integer(sub("^d", "", meta$Day))
meta <- meta[order(meta$Day_num, meta$Treatment != "Control", meta$Run), ]

# Keep only the samples that actually have abundance.h5
meta$h5 <- file.path(meta$Run, "abundance.h5")
falta <- !file.exists(meta$h5)
if (any(falta)) {
  cat("WARNING: no abundance.h5 (skipped):",
      paste(meta$Run[falta], collapse = ", "), "\n")
}
meta <- meta[!falta, ]

if (nrow(meta) == 0) stop("No abundance.h5 files found in ", conteo_root)

# Column name: SRR + group, so the matrix is self-explanatory
sample_names <- paste0(meta$Run, "_", meta$Group)

h5_rutas <- meta$h5
names(h5_rutas) <- sample_names

cat("Found", length(h5_rutas), "quantification files.\n")
cat("Samples:", paste(sample_names, collapse = ", "), "\n\n")

# The SL3.0 index has 34658 targets = 34658 genes (1 transcript per gene),
# so txOut = TRUE already gives gene level and tx2gene is not needed.
txi.kallisto <- tximport(h5_rutas, type = "kallisto", txOut = TRUE)

# The FASTA targets are transcripts (Solyc06g043170.3.1) and the gene is
# Solyc06g043170.3. The transcript suffix is removed so the IDs match gene-level
# matrices. Checked: the conversion is exact for the 34658 targets.
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

# Design table ready for DESeq2 / edgeR
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

cat("Counts saved in:", output_file, "\n")
cat("TPM saved in:", tpm_file, "\n")
cat("Design table in:", coldata_out, "\n")
cat("Dimensions:", nrow(txi.kallisto$counts), "genes x", ncol(txi.kallisto$counts), "samples\n")
