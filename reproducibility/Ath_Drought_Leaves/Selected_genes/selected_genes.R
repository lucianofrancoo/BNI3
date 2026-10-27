#!/usr/bin/env Rscript
# ============================================================
# Selected genes (aba_genes.tsv) – Arabidopsis, rosette leaves, drought
# 1) Normalized count (vst) averaged per condition, listed genes only
# 2) Dot plot of the WD vs CTRL log2FC per time point, with significance
# ============================================================

library(ggplot2)
library(patchwork)

# Run from this folder (Selected_genes/); needs ../DGE/ from DGE/DESeq2.R
base_dir <- ".."
dge_dir  <- file.path(base_dir, "DGE")
out_dir  <- file.path(base_dir, "Selected_genes")

genes <- read.table(file.path(out_dir, "aba_genes.tsv"), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
genes$display_name <- make.unique(genes$geneName, sep = ".")

# -------------------------------------------------------
# Averaged normalized count, selected genes only
# -------------------------------------------------------
escribir_subset <- function(archivo_in, archivo_out, anotacion = TRUE) {
  mat <- read.table(file.path(dge_dir, archivo_in), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  sub <- merge(genes[, c("gene_id", "geneName", "group")], mat,
               by.x = "gene_id", by.y = "ID", sort = FALSE)
  # Keep the order of aba_genes.tsv
  sub <- sub[order(match(sub$gene_id, genes$gene_id)), ]
  colnames(sub)[1] <- "ID"
  # anotacion = FALSE: only ID and expression values (no geneName or group)
  if (!anotacion) sub <- sub[, setdiff(colnames(sub), c("geneName", "group"))]
  write.table(sub, file = file.path(out_dir, archivo_out),
              quote = FALSE, sep = "\t", row.names = FALSE)
  sub
}

sub_mean <- escribir_subset("ConteoNormalizado_mean.tsv", "ConteoNormalizado_mean_aba_genes.tsv")
invisible(escribir_subset("ConteoNormalizado_mean_Control.tsv", "ConteoNormalizado_mean_Control_aba_genes.tsv",
                          anotacion = FALSE))
invisible(escribir_subset("ConteoNormalizado_mean_Drought.tsv", "ConteoNormalizado_mean_Drought_aba_genes.tsv",
                          anotacion = FALSE))

faltan <- genes[!genes$gene_id %in% sub_mean$ID, ]
cat("Genes in the list:", nrow(genes), "\n")
cat("Genes with data (passed the DESeq2 filter):", nrow(sub_mean), "\n")
if (nrow(faltan) > 0) {
  cat("No data (filtered out for low counts):",
      paste0(faltan$geneName, " (", faltan$gene_id, ")", collapse = ", "), "\n")
}

# -------------------------------------------------------
# Load DGE
# -------------------------------------------------------
comp_levels <- c("WD 7200", "WD 15840", "WD 20160")

load_dge <- function(file, label) {
  df <- read.table(file.path(dge_dir, file), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  df$gene_id    <- df$ID
  df$comparison <- label
  df[, c("gene_id", "comparison", "log2FoldChange", "padj")]
}

dge <- rbind(
  load_dge("DGE_WD_7200_vs_CTRL_7200_All.txt",   "WD 7200"),
  load_dge("DGE_WD_15840_vs_CTRL_15840_All.txt", "WD 15840"),
  load_dge("DGE_WD_20160_vs_CTRL_20160_All.txt", "WD 20160")
)

d <- merge(genes[, c("gene_id", "display_name", "group")], dge, by = "gene_id")
d$sig_status <- ifelse(is.na(d$padj), "No data",
                       ifelse(d$padj < 0.05, "FDR < 0.05", "NS"))
d$log2FC_cap <- pmax(pmin(d$log2FoldChange, 2), -2)
d$comparison <- factor(d$comparison, levels = comp_levels)
# Groups in the order they appear in aba_genes.tsv
d$Description <- factor(d$group, levels = unique(genes$group))

# Table with the dot plot values
tabla <- d[order(d$Description, d$display_name, d$comparison),
           c("gene_id", "display_name", "group", "comparison",
             "log2FoldChange", "padj", "sig_status")]
write.table(tabla, file = file.path(out_dir, "dotplot_aba_genes_table.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# -------------------------------------------------------
# Shared scales
# -------------------------------------------------------
fc_color_scale <- scale_color_gradient2(
  low = "#4575B4", mid = "white", high = "#D73027",
  midpoint = 0, limits = c(-2, 2), oob = scales::squish,
  na.value = "grey70", name = "log2FC\n(cap ±2)"
)

sig_shape_scale <- scale_shape_manual(
  values = c("FDR < 0.05" = 16, "NS" = 1),
  limits = c("FDR < 0.05", "NS"),
  name   = "Significance"
)

# Fixed size limits so points are comparable across panels
size_limits <- c(0, max(abs(d$log2FoldChange[d$sig_status != "No data"])))

# -------------------------------------------------------
# Function: dot plot
# -------------------------------------------------------
make_dotplot <- function(data, title, subtitle = "",
                         facet_by_group = TRUE, width = 11) {

  data <- data[data$sig_status != "No data", ]
  data <- data[order(data$Description, data$display_name), ]
  sym_levels <- rev(unique(data$display_name))
  data$gene_f <- factor(data$display_name, levels = sym_levels)

  p <- ggplot(data, aes(
    x     = comparison,
    y     = gene_f,
    size  = abs(log2FoldChange),
    color = log2FC_cap,
    shape = sig_status
  )) +
    geom_point(stroke = 1.2) +
    fc_color_scale +
    sig_shape_scale +
    scale_size_continuous(range = c(2, 9), name = "|log2FC|",
                          breaks = c(1, 2, 3, 4), limits = size_limits) +
    theme_linedraw() +
    theme(
      axis.text.x        = element_text(angle = 40, hjust = 1, size = 11,
                                        face = "bold"),
      axis.text.y        = element_text(size = 8.5),
      strip.text.y       = element_text(size = 9.5, face = "bold",
                                        angle = 0, color = "black"),
      strip.background   = element_rect(fill = "white"),
      legend.title       = element_text(size = 10),
      plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "grey40"),
      panel.grid.major.y = element_line(color = "grey92", linetype = "dashed"),
      panel.grid.major.x = element_blank()
    ) +
    labs(x = "", y = "", title = title, subtitle = subtitle)

  if (facet_by_group) {
    p <- p + facet_grid(Description ~ ., scales = "free_y", space = "free_y")
  }

  n_genes <- length(unique(data$display_name))
  height  <- max(4, n_genes * 0.28 + 2.5)
  list(plot = p, height = height, width = width)
}

save_dotplot <- function(pl, filename_base) {
  svg(file.path(out_dir, paste0(filename_base, ".svg")),
      width = pl$width, height = pl$height)
  print(pl$plot); dev.off()

  png(file.path(out_dir, paste0(filename_base, ".png")),
      width = pl$width, height = pl$height, units = "in", res = 200)
  print(pl$plot); dev.off()

  cat("  → dotplot:", paste0(filename_base, "\n"))
}

# Three panels side by side in the same figure
paneles <- list(
  c("TF", "metabolic"),
  c("catabolic", "transport", "receptor", "phosphatase_cladeA", "phosphatase_cladeB"),
  c("kinase", "ligase", "targets")
)
pls <- lapply(paneles, function(grupos) {
  make_dotplot(d[d$group %in% grupos, ], title = "")
})
fig <- wrap_plots(lapply(pls, `[[`, "plot"), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "ABA Signaling – Drought vs Control (rosette leaf)",
    subtitle = "log2FC WD vs CTRL at each time point",
    theme    = theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"))
  )
save_dotplot(
  list(plot   = fig,
       height = max(sapply(pls, `[[`, "height")),
       width  = 18),
  "dotplot_aba_genes")
