library("BoolNet")

# ==== PARÁMETROS GLOBALES ====
# Run from this folder (BoolNet/)
ruta_base_input  <- "../Selected_genes"
ruta_base_output <- "results"

datasets <- c(
  Control = file.path(ruta_base_input, "ConteoNormalizado_mean_Control_aba_genes_selected_nodes_symbol.tsv"),
  Drought = file.path(ruta_base_input, "ConteoNormalizado_mean_Drought_aba_genes_selected_nodes_symbol.tsv")
)

# k-means tiene inicializaciones aleatorias; la semilla solo fija la reproducibilidad
set.seed(42)

if (!dir.exists(ruta_base_output)) dir.create(ruta_base_output, recursive = TRUE)

cat("========== INFERENCIA DE REDES CON BOOLNET - SLY DROUGHT LEAVES ==========\n\n")
cat("Binarización: binarizeTimeSeries() por defecto (kmeans)\n")
cat("Inferencia:   reconstructNetwork() por defecto (bestfit, maxK = 5)\n")
cat(sprintf("Salida: %s\n\n", ruta_base_output))

# ==== FUNCIÓN PARA PROCESAR UN DATASET ====
procesar_dataset <- function(condicion, ruta_conteos, output_dir) {

  cat(sprintf("\n=== Procesando: %s (%s) ===\n", condicion, basename(ruta_conteos)))

  tryCatch({
    # 1. LEER DATOS
    timeseries <- read.table(ruta_conteos, sep = "\t", header = TRUE,
                             row.names = "ID", check.names = FALSE)
    timeseries <- as.matrix(timeseries)
    cat(sprintf("  Dimensiones: %d genes × %d timepoints\n",
                nrow(timeseries), ncol(timeseries)))

    # 2. BINARIZACIÓN (por defecto: method = "kmeans")
    bin_result <- binarizeTimeSeries(timeseries)

    bin_df <- as.data.frame(bin_result$binarizedMeasurements)
    bin_df <- cbind(Gen = rownames(bin_df), bin_df)
    rownames(bin_df) <- NULL
    write.table(bin_df, file.path(output_dir, paste0("bin_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    write.table(data.frame(Gen = rownames(timeseries), Umbral = bin_result$thresholds),
                file.path(output_dir, paste0("umbrales_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # 3. RECONSTRUIR RED (por defecto: method = "bestfit", maxK = 5)
    # readableFunctions = TRUE solo cambia el formato de la expresión
    red <- reconstructNetwork(bin_result$binarizedMeasurements,
                              readableFunctions = TRUE)

    # 4. EXTRAER REGLAS
    reglas_list <- list()
    todas_list  <- list()
    for (gene in names(red$interactions)) {
      gene_rules <- red$interactions[[gene]]
      if (length(gene_rules) > 0) {
        # Primera regla (la predicha por BoolNet)
        reglas_list[[gene]] <- data.frame(
          Gen = gene,
          Regla = gene_rules[[1]]$expression,
          Error = gene_rules[[1]]$error,
          N_reguladores = sum(gene_rules[[1]]$input != 0),
          N_reglas_alternativas = length(gene_rules),
          stringsAsFactors = FALSE
        )
        # Todas las reglas con el mismo error mínimo
        todas_list[[gene]] <- data.frame(
          Gen = gene,
          Indice = seq_along(gene_rules),
          Regla = sapply(gene_rules, function(r) r$expression),
          Error = sapply(gene_rules, function(r) r$error),
          stringsAsFactors = FALSE
        )
      }
    }
    reglas_final <- do.call(rbind, reglas_list); rownames(reglas_final) <- NULL
    todas_final  <- do.call(rbind, todas_list);  rownames(todas_final)  <- NULL

    # 5. GUARDAR RESULTADOS
    write.table(reglas_final,
                file.path(output_dir, paste0("reglas_BoolNet_bestfit_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(todas_final,
                file.path(output_dir, paste0("reglas_BoolNet_bestfit_", condicion, "_todas.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    saveRDS(list(binarizacion = bin_result, red = red),
            file.path(output_dir, paste0("BoolNet_", condicion, ".rds")))

    # 6. RESUMEN
    cat(sprintf("  ✓ Genes con reglas: %d/%d\n", nrow(reglas_final), nrow(timeseries)))
    cat(sprintf("  ✓ Error promedio: %.4f (min %.4f | max %.4f)\n",
                mean(reglas_final$Error), min(reglas_final$Error), max(reglas_final$Error)))
    cat(sprintf("  ✓ Reglas con error=0: %d (%.1f%%)\n",
                sum(reglas_final$Error == 0),
                100 * mean(reglas_final$Error == 0)))
    TRUE
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    FALSE
  })
}

# ==== PROCESAR CONTROL Y SEQUÍA ====
resultados <- sapply(names(datasets), function(cond)
  procesar_dataset(cond, datasets[[cond]], ruta_base_output))

linea <- strrep("=", 70)
cat("\n", linea, "\nRESUMEN FINAL\n", linea, "\n", sep = "")
cat(sprintf("Datasets procesados con éxito: %d | Errores: %d\n",
            sum(resultados), sum(!resultados)))
cat(sprintf("Resultados guardados en: %s\n", ruta_base_output))
cat(linea, "\nPROCESO COMPLETADO\n", linea, "\n", sep = "")
