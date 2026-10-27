library("BoolNet")

# ==== GLOBAL PARAMETERS ====
# Run from this folder (BoolNet/)
ruta_base_input  <- "../Selected_genes"
ruta_base_output <- "results"

datasets <- c(
  Control = file.path(ruta_base_input, "ConteoNormalizado_mean_Control_aba_genes_selected_nodes_symbol.tsv"),
  Drought = file.path(ruta_base_input, "ConteoNormalizado_mean_Drought_aba_genes_selected_nodes_symbol.tsv")
)

# k-means has random initializations; the seed only makes the run reproducible
set.seed(42)

if (!dir.exists(ruta_base_output)) dir.create(ruta_base_output, recursive = TRUE)

cat("========== NETWORK INFERENCE WITH BOOLNET - SLY DROUGHT LEAVES ==========\n\n")
cat("Binarization: binarizeTimeSeries() default (kmeans)\n")
cat("Inference:     reconstructNetwork() default (bestfit, maxK = 5)\n")
cat(sprintf("Output: %s\n\n", ruta_base_output))

# ==== FUNCTION TO PROCESS ONE DATASET ====
procesar_dataset <- function(condicion, ruta_conteos, output_dir) {

  cat(sprintf("\n=== Processing: %s (%s) ===\n", condicion, basename(ruta_conteos)))

  tryCatch({
    # 1. READ DATA
    timeseries <- read.table(ruta_conteos, sep = "\t", header = TRUE,
                             row.names = "ID", check.names = FALSE)
    timeseries <- as.matrix(timeseries)
    cat(sprintf("  Dimensions: %d genes × %d timepoints\n",
                nrow(timeseries), ncol(timeseries)))

    # 2. BINARIZATION (default: method = "kmeans")
    bin_result <- binarizeTimeSeries(timeseries)

    bin_df <- as.data.frame(bin_result$binarizedMeasurements)
    bin_df <- cbind(Gen = rownames(bin_df), bin_df)
    rownames(bin_df) <- NULL
    write.table(bin_df, file.path(output_dir, paste0("bin_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    write.table(data.frame(Gen = rownames(timeseries), Umbral = bin_result$thresholds),
                file.path(output_dir, paste0("umbrales_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # 3. RECONSTRUCT NETWORK (default: method = "bestfit", maxK = 5)
    # readableFunctions = TRUE only changes the format of the expression
    red <- reconstructNetwork(bin_result$binarizedMeasurements,
                              readableFunctions = TRUE)

    # 4. EXTRACT RULES
    reglas_list <- list()
    todas_list  <- list()
    for (gene in names(red$interactions)) {
      gene_rules <- red$interactions[[gene]]
      if (length(gene_rules) > 0) {
        # First rule (the one predicted by BoolNet)
        reglas_list[[gene]] <- data.frame(
          Gen = gene,
          Regla = gene_rules[[1]]$expression,
          Error = gene_rules[[1]]$error,
          N_reguladores = sum(gene_rules[[1]]$input != 0),
          N_reglas_alternativas = length(gene_rules),
          stringsAsFactors = FALSE
        )
        # All rules with the same minimum error
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

    # 5. SAVE RESULTS
    write.table(reglas_final,
                file.path(output_dir, paste0("reglas_BoolNet_bestfit_", condicion, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(todas_final,
                file.path(output_dir, paste0("reglas_BoolNet_bestfit_", condicion, "_todas.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    saveRDS(list(binarizacion = bin_result, red = red),
            file.path(output_dir, paste0("BoolNet_", condicion, ".rds")))

    # 6. SUMMARY
    cat(sprintf("  ✓ Genes with rules: %d/%d\n", nrow(reglas_final), nrow(timeseries)))
    cat(sprintf("  ✓ Mean error: %.4f (min %.4f | max %.4f)\n",
                mean(reglas_final$Error), min(reglas_final$Error), max(reglas_final$Error)))
    cat(sprintf("  ✓ Rules with error=0: %d (%.1f%%)\n",
                sum(reglas_final$Error == 0),
                100 * mean(reglas_final$Error == 0)))
    TRUE
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    FALSE
  })
}

# ==== PROCESS CONTROL AND DROUGHT ====
resultados <- sapply(names(datasets), function(cond)
  procesar_dataset(cond, datasets[[cond]], ruta_base_output))

linea <- strrep("=", 70)
cat("\n", linea, "\nFINAL SUMMARY\n", linea, "\n", sep = "")
cat(sprintf("Datasets processed successfully: %d | Errors: %d\n",
            sum(resultados), sum(!resultados)))
cat(sprintf("Results saved in: %s\n", ruta_base_output))
cat(linea, "\nDONE\n", linea, "\n", sep = "")
