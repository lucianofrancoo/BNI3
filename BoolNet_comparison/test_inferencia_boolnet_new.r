library("BoolNet")
library("dplyr")

# ==== PARÁMETROS GLOBALES ====
ruta_base_raw <- "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/raw_Data"
ruta_base_gold <- "/home/lahumada/disco1/BNI3/DREAM4/challenge_data/DREAM4_InSilicoNetworks_GoldStandard/DREAM4_Challenge2_GoldStandards/Size 10"
ruta_base_output <- "/home/lahumada/disco1/BNI3/BoolNet_comparison/size10"
metodo <- "bestfit"

# Crear directorio base de salida si no existe
if(!dir.exists(ruta_base_output)) {
  dir.create(ruta_base_output, recursive = TRUE)
}

cat("========== RECONSTRUCCIÓN DE REDES CON BOOLNET - SIZE10 ==========\n\n")
cat(sprintf("Método: %s\n", metodo))
cat(sprintf("Salida: %s\n\n", ruta_base_output))

# ==== FUNCIÓN PARA PROCESAR UN DATASET ====
procesar_dataset <- function(ruta_timeseries, size_num, output_dir) {
  
  nombre_archivo <- basename(ruta_timeseries)
  
  cat(sprintf("\n=== Procesando: %s ===\n", nombre_archivo))
  
  tryCatch({
    # 1. LEER DATOS
    cat("  Leyendo datos de time-series...\n")
    timeseries <- read.table(ruta_timeseries,
                             sep = "\t",
                             header = TRUE,
                             row.names = "ID")
    
    cat(sprintf("  Dimensiones: %d genes × %d timepoints\n", 
                nrow(timeseries), ncol(timeseries)))
    
    # 2. BINARIZACIÓN
    cat("  Binarizando con k-means...\n")
    bin_result <- binarizeTimeSeries(timeseries, method = "kmeans")
    
    output_binary <- file.path(output_dir, paste0("bin_", nombre_archivo, ".tsv"))
    bin_df <- as.data.frame(bin_result$binarizedMeasurements)
    bin_df <- cbind(Gen = rownames(bin_df), bin_df)
    rownames(bin_df) <- NULL
    write.table(bin_df, output_binary, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Matriz binarizada guardada: %s\n", basename(output_binary)))
    
    # 3. RECONSTRUIR RED
    cat("  Reconstruyendo red con BoolNet...\n")
    
    red <- reconstructNetwork(bin_result$binarizedMeasurements,
                              method = metodo,
                              readableFunctions = TRUE)
    
    # 4. EXTRAER REGLAS (UNA POR GEN)
    cat("  Extrayendo reglas...\n")
    reglas_list <- list()
    
    for (gene in names(red$interactions)) {
      gene_rules <- red$interactions[[gene]]
      
      if (length(gene_rules) > 0) {
        # Tomar solo la primera regla (la predicha por BoolNet)
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
    
    # Combinar todas las reglas
    reglas_df <- do.call(rbind, reglas_list)
    rownames(reglas_df) <- NULL
    
    # 5. ORDENAR POR NÚMERO DE GEN (G1, G2, ..., G10)
    extract_gene_number <- function(gene_names) {
      as.numeric(sub("^G", "", gene_names))
    }
    
    reglas_final <- reglas_df %>%
      mutate(Gen_numero = extract_gene_number(Gen)) %>%
      arrange(Gen_numero) %>%
      select(-Gen_numero)
    
    # 6. GUARDAR RESULTADOS
    output_file <- file.path(output_dir, 
                            paste0("reglas_BoolNet_", metodo, "_", nombre_archivo, ".tsv"))
    write.table(reglas_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("  ✓ Reglas guardadas: %s\n", basename(output_file)))
    
    # 7. RESUMEN
    cat(sprintf("  ✓ Genes con reglas: %d/%d\n", 
                nrow(reglas_final), nrow(timeseries)))
    cat(sprintf("  ✓ Error promedio: %.4f\n", mean(reglas_final$Error)))
    cat(sprintf("  ✓ Error mínimo: %.4f | Error máximo: %.4f\n",
                min(reglas_final$Error), max(reglas_final$Error)))
    cat(sprintf("  ✓ Reglas con error=0: %d (%.1f%%)\n",
                sum(reglas_final$Error == 0),
                100 * sum(reglas_final$Error == 0) / nrow(reglas_final)))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    return(FALSE)
  })
}

# ==== PROCESAR TODOS LOS DATASETS ====
cat("Escaneando datasets...\n\n")

total_procesados <- 0
total_errores <- 0

# Iterar por cada size (1 al 5)
for(size_num in 1:5) {
  size_dir_raw <- file.path(ruta_base_raw, paste0("size10_", size_num))
  
  if(!dir.exists(size_dir_raw)) {
    cat(sprintf("ADVERTENCIA: No existe directorio %s\n", size_dir_raw))
    next
  }
  
  # Crear directorio de salida para este size
  output_dir <- file.path(ruta_base_output, paste0("size10_", size_num))
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("\n========== PROCESANDO SIZE10_%d ==========\n", size_num))
  
  # Obtener todos los archivos de timeseries en este directorio
  timeseries_files <- list.files(size_dir_raw, full.names = TRUE)
  
  for(ts_file in timeseries_files) {
    if(file.info(ts_file)$isdir) {
      next  # Saltar si es un directorio
    }
    
    exito <- procesar_dataset(ts_file, size_num, output_dir)
    
    if(exito) {
      total_procesados <- total_procesados + 1
    } else {
      total_errores <- total_errores + 1
    }
  }
}

# ==== RESUMEN FINAL ====
cat("\n")
cat("="*70, "\n", sep="")
cat("RESUMEN FINAL\n")
cat("="*70, "\n", sep="")
cat(sprintf("Total datasets procesados con éxito: %d\n", total_procesados))
cat(sprintf("Total errores: %d\n", total_errores))
cat(sprintf("Tasa de éxito: %.1f%%\n", 
            100 * total_procesados / (total_procesados + total_errores)))
cat("\n")
cat(sprintf("Resultados guardados en: %s\n", ruta_base_output))
cat("\n")
cat("Estructura de salida:\n")
cat("  size10_1/\n")
cat("    ├── bin_size10_1_1.tsv              (matriz binarizada)\n")
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
cat("Formato de archivo de reglas:\n")
cat("  Gen | Regla | Error | N_reguladores\n")
cat("  G1  | G2&G5 | 0.125 | 2\n")
cat("  G2  | G1|G4 | 0.000 | 2\n")
cat("  ...\n")
cat("\n")
cat("="*70, "\n", sep="")
cat("PROCESO COMPLETADO\n")
cat("="*70, "\n", sep="")