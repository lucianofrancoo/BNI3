########################################## Leer Librerías
library(DESeq2)
library(ggplot2)
library(BiocParallel)

register(MulticoreParam(8))

# PRJNA1125861: Solanum lycopersicum, hoja madura, deficit hidrico progresivo.
# Contrastes: WD vs CTRL dentro de cada dia (d4, d5, d8, d10).
# Run from this folder (DGE/); reads ../Counts/
###################################### Leer Datos de conteo
# Columnas del estilo SRR29464690_CTRL_d4_Mature_count (ver Counts/extraerCounts.r)
countdata <- read.table("../Counts/Conteo_genes_Sly_Drought_Leaves.tsv",
                        sep = "\t", header = TRUE, row.names = "gene_id",
                        check.names = FALSE)
# Dejar solo el SRR como nombre de columna para empatar con la metadata
colnames(countdata) <- sub("_.*$", "", colnames(countdata))
countdata <- round(as.matrix(countdata))
head(countdata)
################################ Leer Información de Datos
coldata <- read.table("../raw_data/samples_metadata.tsv",
                      sep = "\t",
                      header = TRUE, row.names = "Run")
coldata
# Ordenar las columnas de countdata según el orden de filas de coldata
order_coldata <- match(rownames(coldata), colnames(countdata))
countdata <- countdata[, order_coldata]
head(countdata)
############################# Convirtiendo la clase de cada columna a factor
dias <- c("d4", "d5", "d8", "d10")
coldata$Day <- factor(coldata$Day, levels = dias)
coldata$Treatment <- factor(ifelse(coldata$Treatment == "Control", "CTRL", "WD"),
                            levels = c("CTRL", "WD"))
coldata$Condition <- factor(paste(coldata$Treatment, coldata$Day, sep = "_"),
                            levels = paste(rep(c("CTRL", "WD"), times = length(dias)),
                                           rep(dias, each = 2), sep = "_"))
coldata$Name <- rownames(coldata)
### Asegurarse que los nombres de filas en coldata
### estén en los nombres de columnas en countData
all(colnames(countdata) %in% rownames(coldata))
################################# Asegurarse que estén en el mismo orden
all(colnames(countdata) == row.names(coldata))
####################################### Construcción del objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ Condition)
############################ Remover conteos con menos de X reads en total
nrow(dds)
# d8 y d10 CTRL tienen solo 2 replicas
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
######################################## Ejecutar DESeq2
de <- DESeq(dds, parallel = TRUE)
################################################# Resultados
contraste <- function(de, contraste_b, contraste_a) {
  resultados <- results(de, contrast = c("Condition",
                                         contraste_b,
                                         contraste_a))
  summary(resultados, alpha = 0.05)
  ############################################### Exportar
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
# WD de cada dia vs su control del mismo dia
resumen <- do.call(rbind, lapply(dias, function(d) {
  contraste(de, paste0("WD_", d), paste0("CTRL_", d))
}))
resumen
write.table(resumen, file = "DGE_summary.tsv",
            quote = FALSE, row.names = FALSE, sep = "\t")
################################################## Barplot resumen
resumen$Label <- sub("^WD_(d[0-9]+) vs CTRL_d[0-9]+$", "WD \\1\nvs CTRL \\1",
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
df_pca$Day <- coldata$Day
svg("PCA.svg", width = 12, height = 10)
ggplot(data = df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 6) +
  geom_text(label = df_pca$Condition, size = 4, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +  # Puedes cambiar "minimal" a otros temas disponibles
  expand_limits(x = c(0, max(df_pca$PC1) * 1.1),
                y = c(0, max(df_pca$PC2) * 1.1)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
dev.off()
################################################## Conteos normalizados
# vst (log2, estabilizado en varianza): para PCA, heatmaps, redes
conteo_normalizado <- assay(vsd)
conteo_normalizado <- cbind(ID = row.names(conteo_normalizado),
                            conteo_normalizado)
write.table(conteo_normalizado,
            file = "ConteoNormalizado.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
# vst promediado por condicion (media de las replicas)
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
# Mismo promedio separado en control y deficit hidrico
for (grupo in c("CTRL", "WD")) {
  columnas <- c("ID", grep(paste0("^", grupo, "_"), colnames(conteo_mean), value = TRUE))
  sufijo <- if (grupo == "CTRL") "Control" else "Drought"
  write.table(conteo_mean[, columnas],
              file = paste0("ConteoNormalizado_mean_", sufijo, ".tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
}
# Normalizados por size factors de DESeq2 (escala de conteos, sin log)
conteo_sf <- counts(de, normalized = TRUE)
conteo_sf <- cbind(ID = row.names(conteo_sf), conteo_sf)
write.table(conteo_sf,
            file = "ConteoNormalizado_sizeFactors.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
