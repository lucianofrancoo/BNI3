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
    cat("  Reconstruyendo red...\n")
    
    red <- reconstructNetwork(bin_result$binarizedMeasurements,
                              method = metodo,
                              readableFunctions = TRUE)
    
    # 4. EXTRAER REGLAS
    cat("  Extrayendo reglas...\n")
    reglas_list <- list()
    
    for (gene in names(red$interactions)) {
      gene_rules <- red$interactions[[gene]]
      
      if (length(gene_rules) > 0) {
        regla_expression <- gene_rules[[1]]$expression
        regla_error <- gene_rules[[1]]$error
        n_inputs <- length(gene_rules[[1]]$input)
        
        reglas_list[[gene]] <- data.frame(
          Gen = gene,
          Regla_ID = paste0(gene, "(1)"),
          Regla = regla_expression,
          Error = regla_error,
          N_reguladores = n_inputs,
          stringsAsFactors = FALSE
        )
      }
    }
    
    reglas_df <- do.call(rbind, reglas_list)
    rownames(reglas_df) <- NULL
    
    # 5. FILTRAR REGLAS INVÁLIDAS
    cat("  Filtrando reglas no interpretables...\n")
    
    is_valid_rule <- function(rule_text) {
      !grepl("<f\\(", rule_text) &
      !grepl("^\\(\\)$", rule_text) &
      rule_text != ""
    }
    
    reglas_validas <- reglas_df %>% filter(is_valid_rule(Regla))
    
    cat(sprintf("  Reglas totales: %d\n", nrow(reglas_df)))
    cat(sprintf("  Reglas válidas: %d\n", nrow(reglas_validas)))
    
    # 6. ORDENAR POR NÚMERO DE GEN
    extract_gene_number <- function(gene_names) {
      as.numeric(sub("^G", "", gene_names))
    }
    
    reglas_final <- reglas_validas %>%
      mutate(Gen_numero = extract_gene_number(Gen)) %>%
      arrange(Gen_numero) %>%
      select(-Gen_numero)
    
    # 7. GUARDAR RESULTADOS
    output_file <- file.path(output_dir, 
                            paste0("reglas_Boolnet_", metodo, "_", nombre_archivo, ".tsv"))
    write.table(reglas_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("  Reglas guardadas: %s\n", basename(output_file)))
    
    # 8. RESUMEN
    cat(sprintf("  Genes con reglas: %d | Error promedio: %.2f | Error=0: %d\n",
                nrow(reglas_final),
                mean(reglas_final$Error),
                sum(reglas_final$Error == 0)))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
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
cat("\n========== RESUMEN FINAL ==========\n")
cat(sprintf("Total datasets procesados: %d\n", total_procesados))
cat(sprintf("Total errores: %d\n", total_errores))
cat(sprintf("\nResultados guardados en: %s\n", ruta_base_output))
cat("\nEstructura de salida:\n")
cat("  size10_1/\n")
cat("    ├── bin_size10_1_1.tsv\n")
cat("    ├── reglas_Boolnet_bestfit_size10_1_1.tsv\n")
cat("    └── ...\n")
cat("  size10_2/\n")
cat("    └── ...\n")
cat("  size10_3/\n")
cat("    └── ...\n")
cat("  size10_4/\n")
cat("    └── ...\n")
cat("  size10_5/\n")
cat("    └── ...\n")

cat("\n========== PROCESO COMPLETADO ==========\n")