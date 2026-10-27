library("BoolNet")
library("dplyr")

# BoolNet (reconstructNetwork, bestfit) on the DREAM4 size10 series.
# Usage: Rscript run_boolnet_dream4.r [output_dir]
#   output_dir defaults to ../rerun/BoolNet so the paper's results (../results/BoolNet)
#   are not overwritten. Evaluate afterwards with eval_gold_standard_boolnet.r.

# ==== GLOBAL PARAMETERS ====
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if(length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
dream4_dir <- dirname(script_dir)
args <- commandArgs(trailingOnly = TRUE)
ruta_base_raw <- file.path(dream4_dir, "data", "time_series")
ruta_base_gold <- file.path(dream4_dir, "data", "gold_standard")
ruta_base_output <- if(length(args) >= 1) args[1] else file.path(dream4_dir, "rerun", "BoolNet")
metodo <- "bestfit"

# Create the base output folder if needed
if(!dir.exists(ruta_base_output)) {
  dir.create(ruta_base_output, recursive = TRUE)
}

cat("========== NETWORK RECONSTRUCTION WITH BOOLNET - SIZE10 ==========\n\n")
cat(sprintf("Method: %s\n", metodo))
cat(sprintf("Output: %s\n\n", ruta_base_output))

# ==== FUNCTION TO PROCESS ONE DATASET ====
procesar_dataset <- function(ruta_timeseries, size_num, output_dir) {
  
  nombre_archivo <- basename(ruta_timeseries)
  
  cat(sprintf("\n=== Processing: %s ===\n", nombre_archivo))
  
  tryCatch({
    # 1. READ DATA
    cat("  Reading time series...\n")
    timeseries <- read.table(ruta_timeseries,
                             sep = "\t",
                             header = TRUE,
                             row.names = "ID")
    
    cat(sprintf("  Dimensions: %d genes × %d timepoints\n", 
                nrow(timeseries), ncol(timeseries)))
    
    # 2. BINARIZATION
    cat("  Binarizing with k-means...\n")
    bin_result <- binarizeTimeSeries(timeseries, method = "kmeans")
    
    output_binary <- file.path(output_dir, paste0("bin_", nombre_archivo, ".tsv"))
    bin_df <- as.data.frame(bin_result$binarizedMeasurements)
    bin_df <- cbind(Gen = rownames(bin_df), bin_df)
    rownames(bin_df) <- NULL
    write.table(bin_df, output_binary, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Binarized matrix saved: %s\n", basename(output_binary)))
    
    # 3. RECONSTRUCT NETWORK
    cat("  Reconstructing the network with BoolNet...\n")
    
    red <- reconstructNetwork(bin_result$binarizedMeasurements,
                              method = metodo,
                              readableFunctions = TRUE)
    
    # 4. EXTRACT RULES (ONE PER GENE)
    cat("  Extracting rules...\n")
    reglas_list <- list()
    
    for (gene in names(red$interactions)) {
      gene_rules <- red$interactions[[gene]]
      
      if (length(gene_rules) > 0) {
        # Take only the first rule (the one predicted by BoolNet)
        regla_expression <- gene_rules[[1]]$expression
        regla_error <- gene_rules[[1]]$error
        n_inputs <- length(gene_rules[[1]]$input)
        
        reglas_list[[gene]] <- data.frame(
          Gen = gene,
          Regla = regla_expression,
          Error = regla_error,
          N_reguladores = n_inputs,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Combine all rules
    reglas_df <- do.call(rbind, reglas_list)
    rownames(reglas_df) <- NULL
    
    # 5. SORT BY GENE NUMBER (G1, G2, ..., G10)
    extract_gene_number <- function(gene_names) {
      as.numeric(sub("^G", "", gene_names))
    }
    
    reglas_final <- reglas_df %>%
      mutate(Gen_numero = extract_gene_number(Gen)) %>%
      arrange(Gen_numero) %>%
      select(-Gen_numero)
    
    # 6. SAVE RESULTS
    output_file <- file.path(output_dir, 
                            paste0("reglas_BoolNet_", metodo, "_", nombre_archivo, ".tsv"))
    write.table(reglas_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("  ✓ Rules saved: %s\n", basename(output_file)))
    
    # 7. SUMMARY
    cat(sprintf("  ✓ Genes with rules: %d/%d\n", 
                nrow(reglas_final), nrow(timeseries)))
    cat(sprintf("  ✓ Mean error: %.4f\n", mean(reglas_final$Error)))
    cat(sprintf("  ✓ Min error: %.4f | Max error: %.4f\n",
                min(reglas_final$Error), max(reglas_final$Error)))
    cat(sprintf("  ✓ Rules with error=0: %d (%.1f%%)\n",
                sum(reglas_final$Error == 0),
                100 * sum(reglas_final$Error == 0) / nrow(reglas_final)))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    return(FALSE)
  })
}

# ==== PROCESS ALL DATASETS ====
cat("Scanning datasets...\n\n")

total_procesados <- 0
total_errores <- 0

# Loop over each network (1 to 5)
for(size_num in 1:5) {
  size_dir_raw <- file.path(ruta_base_raw, paste0("size10_", size_num))
  
  if(!dir.exists(size_dir_raw)) {
    cat(sprintf("WARNING: directory does not exist: %s\n", size_dir_raw))
    next
  }
  
  # Output folder for this network
  output_dir <- file.path(ruta_base_output, paste0("size10_", size_num))
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("\n========== PROCESSING SIZE10_%d ==========\n", size_num))
  
  # All time series files in this folder
  timeseries_files <- list.files(size_dir_raw, full.names = TRUE)
  
  for(ts_file in timeseries_files) {
    if(file.info(ts_file)$isdir) {
      next  # Skip directories
    }
    
    exito <- procesar_dataset(ts_file, size_num, output_dir)
    
    if(exito) {
      total_procesados <- total_procesados + 1
    } else {
      total_errores <- total_errores + 1
    }
  }
}

# ==== FINAL SUMMARY ====
cat("\n")
cat(strrep("=", 70), "\n", sep="")
cat("FINAL SUMMARY\n")
cat(strrep("=", 70), "\n", sep="")
cat(sprintf("Datasets processed successfully: %d\n", total_procesados))
cat(sprintf("Errors: %d\n", total_errores))
cat(sprintf("Success rate: %.1f%%\n", 
            100 * total_procesados / (total_procesados + total_errores)))
cat("\n")
cat(sprintf("Results saved in: %s\n", ruta_base_output))
cat("\n")
cat("Output layout:\n")
cat("  size10_1/\n")
cat("    ├── bin_size10_1_1.tsv              (binarized matrix)\n")
cat("    ├── reglas_BoolNet_bestfit_size10_1_1.tsv\n")
cat("    ├── bin_size10_1_2.tsv\n")
cat("    ├── reglas_BoolNet_bestfit_size10_1_2.tsv\n")
cat("    └── ...\n")
cat("  size10_2/\n")
cat("    └── ...\n")
cat("  size10_3/\n")
cat("    └── ...\n")
cat("  size10_4/\n")
cat("    └── ...\n")
cat("  size10_5/\n")
cat("    └── ...\n")
cat("\n")
cat("Rules file columns:\n")
cat("  Gen | Regla | Error | N_reguladores\n")
cat("  G1  | G2&G5 | 0.125 | 2\n")
cat("  G2  | G1|G4 | 0.000 | 2\n")
cat("  ...\n")
cat("\n")
cat(strrep("=", 70), "\n", sep="")
cat("DONE\n")
cat(strrep("=", 70), "\n", sep="")