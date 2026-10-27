########################################## Load libraries
library(DESeq2)
library(ggplot2)
library(BiocParallel)

register(MulticoreParam(8))

# Arabidopsis thaliana, rosette leaves, drought.
# Contrasts: Drought vs Control within each time point (7200, 15840, 20160).
# Control_0 (Treatment = None) enters the model and the PCA, but not the contrasts.
# Run from this folder (DGE/); reads ../Counts/
###################################### Read count data
countdata <- read.table("../Counts/Counts.tsv",
                        sep = "\t", header = TRUE, row.names = "ID",
                        check.names = FALSE)
countdata <- round(as.matrix(countdata))
head(countdata)
################################ Read sample information
coldata <- read.table("../Counts/Counts_Info.tsv",
                      sep = "\t",
                      header = TRUE, row.names = "Name")
coldata
# Order the countdata columns by the row order of coldata
order_coldata <- match(rownames(coldata), colnames(countdata))
countdata <- countdata[, order_coldata]
head(countdata)
############################# Convert each column to factor
tiempos <- c("7200", "15840", "20160")
coldata$Time <- factor(coldata$Time, levels = c("0", tiempos))
# "None" (time 0) is the initial control; "Drougth" is spelled this way in Counts_Info.tsv
coldata$Treatment <- factor(ifelse(coldata$Treatment == "Drougth", "WD", "CTRL"),
                            levels = c("CTRL", "WD"))
coldata$Condition <- factor(paste(coldata$Treatment, coldata$Time, sep = "_"),
                            levels = c("CTRL_0",
                                       paste(rep(c("CTRL", "WD"), times = length(tiempos)),
                                             rep(tiempos, each = 2), sep = "_")))
coldata$Name <- rownames(coldata)
### Make sure the row names of coldata
### are among the column names of countData
all(colnames(countdata) %in% rownames(coldata))
################################# Make sure they are in the same order
all(colnames(countdata) == row.names(coldata))
####################################### Build the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ Condition)
############################ Remove genes with fewer than X reads in total
nrow(dds)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
######################################## Run DESeq2
de <- DESeq(dds, parallel = TRUE)
################################################# Results
contraste <- function(de, contraste_b, contraste_a) {
  resultados <- results(de, contrast = c("Condition",
                                         contraste_b,
                                         contraste_a))
  summary(resultados, alpha = 0.05)
  ############################################### Export
  resultados <- resultados[order(resultados$padj), ]
  resultados <- cbind(ID = row.names(resultados),
                      as.data.frame(resultados))
  write.table(resultados,
              file = paste0("DGE_", contraste_b, "_vs_", contraste_a, "_All.txt"), # nolint: line_length_linter.
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
  sig <- subset(resultados, !is.na(padj) & padj < 0.05) # nolint
  data.frame(Comparison = paste0(contraste_b, " vs ", contraste_a),
             Up = sum(sig$log2FoldChange > 0),
             Down = sum(sig$log2FoldChange < 0))
}
coldata
# WD at each time point vs the control at the same time point
resumen <- do.call(rbind, lapply(tiempos, function(t) {
  contraste(de, paste0("WD_", t), paste0("CTRL_", t))
}))
resumen
write.table(resumen, file = "DGE_summary.tsv",
            quote = FALSE, row.names = FALSE, sep = "\t")
################################################## Summary barplot
resumen$Label <- sub("^WD_([0-9]+) vs CTRL_[0-9]+$", "WD \\1\nvs CTRL \\1",
                     resumen$Comparison)
resumen$Label <- factor(resumen$Label, levels = resumen$Label)
df_bar <- rbind(data.frame(Label = resumen$Label, n = resumen$Up,
                           Direction = "Up-regulated"),
                data.frame(Label = resumen$Label, n = -resumen$Down,
                           Direction = "Down-regulated"))
df_bar$Direction <- factor(df_bar$Direction,
                           levels = c("Up-regulated", "Down-regulated"))
lim <- max(abs(df_bar$n)) * 1.18
p_bar <- ggplot(df_bar, aes(x = Label, y = n, fill = Direction)) +
  geom_col(width = 0.55) +
  geom_text(data = subset(df_bar, n != 0),
            aes(label = abs(n), color = Direction,
                vjust = ifelse(n > 0, -0.4, 1.4)),
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Up-regulated" = "#E05252",
                               "Down-regulated" = "#5278E0")) +
  scale_color_manual(values = c("Up-regulated" = "#C03030",
                                "Down-regulated" = "#2040B0")) +
  scale_y_continuous(limits = c(-lim, lim), labels = function(x) abs(x)) +
  labs(x = NULL, y = "Number of DEGs", fill = NULL,
       title = "Differentially Expressed Genes per Comparison\n(padj < 0.05)") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = c(0.92, 0.95),
        panel.grid.major.y = element_line(linetype = "dashed",
                                          linewidth = 0.3, color = "grey80"))
svg("DGE_summary_barplot.svg", width = 10, height = 6)
print(p_bar)
dev.off()
png("DGE_summary_barplot.png", width = 10, height = 6, units = "in", res = 150)
print(p_bar)
dev.off()
################################################## PCA plot
vsd <- vst(de)
plotPCA(vsd, intgroup = c("Condition"))
df_pca <- plotPCA(vsd, intgroup = c("Condition"), returnData = TRUE)
percentVar <- round(100 * attr(df_pca, "percentVar"))
df_pca$Treatment <- coldata$Treatment
df_pca$Name <- coldata$Name
df_pca$Time <- coldata$Time
svg("PCA.svg", width = 12, height = 10)
ggplot(data = df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 6) +
  geom_text(label = df_pca$Condition, size = 4, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  expand_limits(x = c(0, max(df_pca$PC1) * 1.1),
                y = c(0, max(df_pca$PC2) * 1.1)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
dev.off()
################################################## Normalized counts
# vst (log2, variance-stabilized): for PCA, heatmaps, networks
conteo_normalizado <- assay(vsd)
conteo_normalizado <- cbind(ID = row.names(conteo_normalizado),
                            conteo_normalizado)
write.table(conteo_normalizado,
            file = "ConteoNormalizado.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
# vst averaged per condition (mean of the replicates)
vst_mat <- assay(vsd)
conteo_mean <- sapply(levels(coldata$Condition), function(cond) {
  rowMeans(vst_mat[, coldata$Condition == cond, drop = FALSE])
})
colnames(conteo_mean) <- paste0(colnames(conteo_mean), "_mean")
conteo_mean <- cbind(ID = row.names(conteo_mean), conteo_mean)
write.table(conteo_mean,
            file = "ConteoNormalizado_mean.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
# Same average split into control and water deficit. CTRL_0 is the starting point
# of both series, so it is the first column in both files.
for (grupo in c("CTRL", "WD")) {
  columnas <- unique(c("ID", "CTRL_0_mean",
                       grep(paste0("^", grupo, "_"), colnames(conteo_mean), value = TRUE)))
  sufijo <- if (grupo == "CTRL") "Control" else "Drought"
  write.table(conteo_mean[, columnas],
              file = paste0("ConteoNormalizado_mean_", sufijo, ".tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
}
# Normalized by the DESeq2 size factors (count scale, no log)
conteo_sf <- counts(de, normalized = TRUE)
conteo_sf <- cbind(ID = row.names(conteo_sf), conteo_sf)
write.table(conteo_sf,
            file = "ConteoNormalizado_sizeFactors.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
