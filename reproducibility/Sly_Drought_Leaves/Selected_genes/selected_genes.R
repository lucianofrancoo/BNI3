#!/usr/bin/env Rscript
# ============================================================
# Tomato orthologs of the Arabidopsis ABA genes – mature leaves, drought
# (PRJNA1125861). Same analysis as
# Ath_Drought_Leaves/Selected_genes/selected_genes.R:
# 1) Normalized count (vst) averaged per condition, selected genes only
# 2) Dot plot of the WD vs CTRL log2FC per day, with significance
# The tomato genes come from TAIR_To_ITAG/Table_S_ABA_orthologs.tsv
# (orthologs from OrthoFinder or Ensembl Compara; see TAIR_To_ITAG/README.md).
# ============================================================

library(ggplot2)
library(patchwork)

# Run from this folder (Selected_genes/); needs ../DGE/ from DGE/DESeq2.R
base_dir  <- ".."
dge_dir   <- file.path(base_dir, "DGE")
out_dir   <- file.path(base_dir, "Selected_genes")
orth_file <- "../../TAIR_To_ITAG/Table_S_ABA_orthologs.tsv"

# -------------------------------------------------------
# Gene list: one tomato gene per row, with its Arabidopsis genes
# -------------------------------------------------------
orth <- read.table(orth_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
orth <- orth[orth$tomato_gene != "No ortholog", ]
# Groups with the original names in aba_genes.tsv (those of the Arabidopsis dot plot)
group_key <- c("TF" = "TF", "Metabolic" = "metabolic", "Catabolic" = "catabolic",
               "Transport" = "transport", "Receptor" = "receptor",
               "Phosphatase clade A" = "phosphatase_cladeA",
               "Phosphatase clade B" = "phosphatase_cladeB", "Kinase" = "kinase",
               "Ligase" = "ligase", "Targets" = "targets")
orth$group <- unname(group_key[orth$functional_group])

sly_ids <- unique(orth$tomato_gene)   # keeps the table order (order of the Ath list)
genes <- data.frame(
  gene_id = sly_ids,
  ath_genes = sapply(sly_ids, function(s)
    paste(unique(orth$arabidopsis_gene[orth$tomato_gene == s]), collapse = "/")),
  ath_loci = sapply(sly_ids, function(s)
    paste(unique(orth$arabidopsis_locus[orth$tomato_gene == s]), collapse = "/")),
  group = sapply(sly_ids, function(s) orth$group[orth$tomato_gene == s][1]),
  stringsAsFactors = FALSE, row.names = NULL
)
genes$geneName <- paste0(genes$gene_id, " (", genes$ath_genes, ")")
genes$display_name <- genes$geneName
write.table(genes[, c("gene_id", "ath_genes", "ath_loci", "group")],
            file = file.path(out_dir, "aba_orthologs_sly.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# -------------------------------------------------------
# Averaged normalized count, selected genes only
# -------------------------------------------------------
escribir_subset <- function(archivo_in, archivo_out, anotacion = TRUE) {
  mat <- read.table(file.path(dge_dir, archivo_in), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  sub <- merge(genes[, c("gene_id", "ath_genes", "group")], mat,
               by.x = "gene_id", by.y = "ID", sort = FALSE)
  # Keep the order of the gene list
  sub <- sub[order(match(sub$gene_id, genes$gene_id)), ]
  colnames(sub)[1] <- "ID"
  # anotacion = FALSE: only ID and expression values (no ath_genes or group)
  if (!anotacion) sub <- sub[, setdiff(colnames(sub), c("ath_genes", "group"))]
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
cat("Tomato genes in the list:", nrow(genes), "\n")
cat("Genes with data (passed the DESeq2 filter):", nrow(sub_mean), "\n")
if (nrow(faltan) > 0) {
  cat("No data (filtered out for low counts):",
      paste0(faltan$gene_id, " (", faltan$ath_genes, ")", collapse = ", "), "\n")
}

# -------------------------------------------------------
# Load DGE
# -------------------------------------------------------
comp_levels <- c("WD d4", "WD d5", "WD d8", "WD d10")

load_dge <- function(file, label) {
  df <- read.table(file.path(dge_dir, file), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  df$gene_id    <- df$ID
  df$comparison <- label
  df[, c("gene_id", "comparison", "log2FoldChange", "padj")]
}

dge <- rbind(
  load_dge("DGE_WD_d4_vs_CTRL_d4_All.txt",   "WD d4"),
  load_dge("DGE_WD_d5_vs_CTRL_d5_All.txt",   "WD d5"),
  load_dge("DGE_WD_d8_vs_CTRL_d8_All.txt",   "WD d8"),
  load_dge("DGE_WD_d10_vs_CTRL_d10_All.txt", "WD d10")
)

d <- merge(genes[, c("gene_id", "display_name", "ath_genes", "group")], dge, by = "gene_id")
d$sig_status <- ifelse(is.na(d$padj), "No data",
                       ifelse(d$padj < 0.05, "FDR < 0.05", "NS"))
d$log2FC_cap <- pmax(pmin(d$log2FoldChange, 2), -2)
d$comparison <- factor(d$comparison, levels = comp_levels)
# Groups in the order of the list
d$Description <- factor(d$group, levels = unique(genes$group))
# Rows in the order of the list (not alphabetical), so paralogs stay together
d$display_name <- factor(d$display_name, levels = genes$display_name)

# Table with the dot plot values
tabla <- d[order(d$Description, d$display_name, d$comparison),
           c("gene_id", "ath_genes", "group", "comparison",
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
  sym_levels <- rev(unique(as.character(data$display_name)))
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

# Three panels side by side in the same figure (same groups as in Arabidopsis)
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
    title    = "ABA Signaling (tomato orthologs) – Drought vs Control (mature leaf)",
    subtitle = "log2FC WD vs CTRL at each day; labels: tomato gene (Arabidopsis ortholog)",
    theme    = theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"))
  )
save_dotplot(
  list(plot   = fig,
       height = max(sapply(pls, `[[`, "height")),
       width  = 22),
  "dotplot_aba_genes")
