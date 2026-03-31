library("BoolNet")
library("dplyr")
library("igraph")

# ==== PARÁMETROS GLOBALES ====
ruta_base_raw <- "/home/lahumada/disco1/Valdivia2025/counts/Counts_RosetteLeaf_Control_reduced_final.tsv"
ruta_base_output <- "/home/lahumada/disco1/BNI3/Rules_Inference/BoolNet"

cat(sprintf("Salida: %s\n\n", ruta_base_output))

cat("Leyendo datos de time-series...\n")
timeseries <- read.table(ruta_base_raw,
                         sep = "\t",
                         header = TRUE,
                         row.names = "ID")
cat(sprintf("  Dimensiones: %d genes × %d timepoints\n", nrow(timeseries), ncol(timeseries)))

cat("  Binarizando con k-means...\n")
bin_result <- binarizeTimeSeries(timeseries, method = "kmeans")

nombre_archivo <- basename(ruta_base_raw)
output_binary <- file.path(ruta_base_output, paste0("bin_", nombre_archivo))
bin_df <- as.data.frame(bin_result$binarizedMeasurements)
bin_df <- cbind(Gen = rownames(bin_df), bin_df)
rownames(bin_df) <- NULL
write.table(bin_df, output_binary, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Matriz binarizada guardada: %s\n", basename(output_binary)))

cat("  Reconstruyendo red con BoolNet...\n")
red <- reconstructNetwork(bin_result$binarizedMeasurements,
                          method = "bestfit",
                          readableFunctions = TRUE)

# EXTRAER REGLAS
cat("  Extrayendo reglas...\n")
reglas_list <- list()

for (gene in names(red$interactions)) {
  gene_rules <- red$interactions[[gene]]
  
  if (length(gene_rules) > 0) {
    regla_expression <- gene_rules[[1]]$expression
    # Remover paréntesis externos
    regla_expression <- gsub("^\\((.*)\\)$", "\\1", regla_expression)
    
    # Reemplazar 1 por TRUE y 0 por FALSE
    if (regla_expression == "1") {
      regla_expression <- "TRUE"
    } else if (regla_expression == "0") {
      regla_expression <- "FALSE"
    }
    
    n_inputs <- length(gene_rules[[1]]$input)
    
    reglas_list[[gene]] <- data.frame(
      Gene = gene,
      Position = 1,
      Rule = regla_expression,
      Correct = "N/A",
      N_Regulators = n_inputs,
      MSE = "N/A",
      Score = "N/A",
      stringsAsFactors = FALSE
    )
  }
}

reglas_df <- do.call(rbind, reglas_list)
rownames(reglas_df) <- NULL

# Ordenar alfabéticamente por Gene
reglas_df <- reglas_df %>%
  arrange(Gene)

output_file <- file.path(ruta_base_output, paste0("reglas_BoolNet_", "bestfit", "_", nombre_archivo))
write.table(reglas_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  ✓ Reglas guardadas: %s\n", basename(output_file)))

# CREAR RED ÚNICA
singlenet <- chooseNetwork(red,
                           functionIndices = rep(1, 13),
                           readableFunctions = TRUE)

# PLOT DE LA RED (guardar como PNG)
output_network_plot <- file.path(ruta_base_output, paste0("network_plot_", 
                                                           sub("\\.tsv$", ".png", nombre_archivo)))
png(output_network_plot, width = 800, height = 800, res = 100)
plot_obj <- plotNetworkWiring(singlenet, plotIt = FALSE)
plot.igraph(plot_obj,
            vertex.label.cex = 2,
            vertex.label.dist = 2,
            vertex.label.font = 2.5,
            edge.arrow.size = 1,
            edge.color = "black",
            loop.size = 0.5,
            layout = layout.kamada.kawai,
            vertex.label.color = "black")
dev.off()
cat(sprintf("  ✓ Network plot guardado: %s\n", basename(output_network_plot)))

# CALCULAR ATRACTORES
cat("  Calculando atractores...\n")
attractors <- getAttractors(network = singlenet,
                            type = "synchronous",
                            method = "exhaustive",
                            returnTable = TRUE)

# GUARDAR INFORMACIÓN DE ATRACTORES COMO TEXTO
output_attractors_txt <- file.path(ruta_base_output, 
                                   paste0("attractors_info_", 
                                          sub("\\.tsv$", ".txt", nombre_archivo)))

sink(output_attractors_txt)
cat("=", rep("=", 70), "\n", sep = "")
cat("INFORMACIÓN DE ATRACTORES\n")
cat("=", rep("=", 70), "\n\n", sep = "")
cat(sprintf("Archivo de entrada: %s\n", nombre_archivo))
cat(sprintf("Fecha de análisis: %s\n\n", Sys.time()))
print(attractors)
cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
sink()

cat(sprintf("  ✓ Información de atractores guardada: %s\n", basename(output_attractors_txt)))

# PLOT DEL STATE GRAPH (guardar como PNG)
output_state_graph <- file.path(ruta_base_output, 
                                paste0("state_graph_", 
                                       sub("\\.tsv$", ".png", nombre_archivo)))
png(output_state_graph, width = 1200, height = 1000, res = 100)
p <- plotStateGraph(attractors, plotIt = TRUE)
dev.off()
cat(sprintf("  ✓ State graph guardado: %s\n", basename(output_state_graph)))

# RESUMEN FINAL
cat("\n", rep("=", 70), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS\n")
cat(rep("=", 70), "\n\n", sep = "")
cat(sprintf("✓ Reglas: %s\n", basename(output_file)))
cat(sprintf("✓ Network plot: %s\n", basename(output_network_plot)))
cat(sprintf("✓ Attractors info: %s\n", basename(output_attractors_txt)))
cat(sprintf("✓ State graph: %s\n", basename(output_state_graph)))
cat("\nAnálisis completado exitosamente.\n")