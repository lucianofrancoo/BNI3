library("dplyr")

# ==== PARÁMETROS GLOBALES ====
SIZE <- 10  # <<< CAMBIAR AQUÍ PARA EVALUAR DIFERENTES TAMAÑOS

ruta_base_reglas <- sprintf("/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size%d/rules", SIZE)
ruta_base_gold <- sprintf("/home/lahumada/disco1/BNI3/DREAM4/challenge_data/DREAM4_InSilicoNetworks_GoldStandard/DREAM4_Challenge2_GoldStandards/Size %d", SIZE)
ruta_base_bin <- sprintf("/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size%d/bin_Data", SIZE)

# Establecer directorio de trabajo
setwd(ruta_base_reglas)

# Usar la carpeta rules directamente como output
output_dir <- ruta_base_reglas

cat(sprintf("========== EVALUACIÓN CONSOLIDADA BNI3 - DREAM4 SIZE%d ==========\n", SIZE))
cat(sprintf("Ruta base reglas: %s\n", ruta_base_reglas))
cat(sprintf("Ruta base gold:   %s\n", ruta_base_gold))
cat(sprintf("Ruta base bin:    %s\n", ruta_base_bin))
cat(sprintf("Output dir:       %s\n\n", output_dir))

# ==== FUNCIONES AUXILIARES ====

extract_regulators <- function(rule_text, target_gene) {
  cleaned <- gsub("~", " ", rule_text)
  cleaned <- gsub("!", " ", cleaned)
  cleaned <- gsub("\\(|\\)|&|\\|", " ", cleaned)
  genes <- unlist(strsplit(cleaned, "\\s+"))
  genes <- genes[grepl("^G[0-9]+$", genes)]
  genes <- unique(genes)
  genes <- genes[genes != target_gene]
  return(genes)
}

evaluate_boolean_rule <- function(rule_text, state_vector) {
  expression <- rule_text
  expression <- gsub("~", "!", expression)
  
  for(gene_name in names(state_vector)) {
    value <- state_vector[gene_name]
    expression <- gsub(paste0("\\b", gene_name, "\\b"), 
                      as.character(value), 
                      expression)
  }
  
  expression <- gsub("&", "&&", expression, fixed = TRUE)
  expression <- gsub("|", "||", expression, fixed = TRUE)
  
  tryCatch({
    result <- eval(parse(text = expression))
    return(as.integer(result))
  }, error = function(e) {
    return(NA)
  })
}

# Función para extraer información del nombre del método
parse_method_name <- function(method_name, size) {
  # Extraer método de binarización (SSD o WCSS)
  bin_method <- ifelse(grepl("SSD", method_name), "SSD", "WCSS")
  
  # Patrón genérico que se adapta al tamaño
  # Formato: rules_[METHOD]_size[SIZE]_[NETWORK]_[REP]
  pattern <- sprintf("rules_(SSD|WCSS)_size%d_(\\d+)_(\\d+)(_5times)?$", size)
  matches <- regexec(pattern, method_name)
  match_data <- regmatches(method_name, matches)[[1]]
  
  if(length(match_data) >= 4) {
    network <- as.integer(match_data[3])
    rep_base <- as.integer(match_data[4])
    has_5times <- length(match_data) >= 5 && match_data[5] == "_5times"
    
    # Construir el string de Rep
    if(has_5times) {
      rep_str <- paste0(rep_base, "_5times")
    } else {
      rep_str <- as.character(rep_base)
    }
    
    return(list(
      bin_method = bin_method,
      network = network,
      rep = rep_str,
      rep_numeric = rep_base  # Para ordenamiento
    ))
  }
  
  return(NULL)
}

# ==== FUNCIÓN PRINCIPAL DE EVALUACIÓN ====
evaluar_red <- function(ruta_reglas, size, size_num, method_name, ruta_gold_standard, ruta_matriz_bin) {
  
  cat(sprintf("\n--- Procesando: size%d_%d - %s ---\n", size, size_num, method_name))
  
  # Verificar que existe el archivo
  if(!file.exists(ruta_reglas)) {
    cat(sprintf("  ADVERTENCIA: No se encontró %s\n", ruta_reglas))
    return(NULL)
  }
  
  tryCatch({
    # 1. LEER REGLAS
    reglas_completas <- read.table(ruta_reglas,
                                   sep = "\t",
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
    
    reglas <- reglas_completas %>%
      filter(Position == 1) %>%
      select(Gene, Rule) %>%
      rename(Gen = Gene, Regla = Rule)
    
    # Extraer estadísticas de las reglas si existen
    score_prom <- NA
    mse_prom <- NA
    correct_prom <- NA
    
    if("Score" %in% colnames(reglas_completas)) {
      reglas_stats <- reglas_completas %>% filter(Position == 1)
      score_prom <- mean(reglas_stats$Score, na.rm = TRUE)
      
      # Manejo especial para MSE (puede ser no numérico en archivos _5times)
      if("MSE" %in% colnames(reglas_stats)) {
        if(is.numeric(reglas_stats$MSE)) {
          mse_prom <- mean(reglas_stats$MSE, na.rm = TRUE)
        } else {
          mse_prom <- NA
        }
      }
      
      correct_prom <- mean(reglas_stats$Correct, na.rm = TRUE)
    }
    
    # 2. EXTRAER REGULADORES
    reglas$Reguladores <- lapply(seq_len(nrow(reglas)), function(i) {
      extract_regulators(reglas$Regla[i], reglas$Gen[i])
    })
    
    n_reguladores <- sapply(reglas$Reguladores, length)
    
    # 3. LEER GOLD STANDARD
    gold_standard <- read.table(ruta_gold_standard,
                                sep = "\t",
                                header = FALSE,
                                col.names = c("Regulador", "Target", "Interaccion"),
                                stringsAsFactors = FALSE)
    
    gold_standard_positive <- gold_standard %>%
      filter(Interaccion == 1) %>%
      select(Regulador, Target)
    
    # 4. CONSTRUIR PREDICCIONES
    predicciones <- data.frame()
    
    for(i in seq_len(nrow(reglas))) {
      target <- reglas$Gen[i]
      reguladores <- reglas$Reguladores[[i]]
      
      if(length(reguladores) > 0) {
        for(reg in reguladores) {
          predicciones <- rbind(predicciones, 
                               data.frame(Regulador = reg, 
                                         Target = target,
                                         stringsAsFactors = FALSE))
        }
      }
    }
    
    # 5. CALCULAR MÉTRICAS ESTRUCTURALES
    all_genes <- unique(c(gold_standard$Regulador, gold_standard$Target, reglas$Gen))
    
    TP <- 0
    FP <- 0
    FN <- 0
    TN <- 0
    
    for(regulador in all_genes) {
      for(target in all_genes) {
        if(regulador == target) next
        
        en_gold <- any(gold_standard_positive$Regulador == regulador & 
                       gold_standard_positive$Target == target)
        
        en_pred <- any(predicciones$Regulador == regulador & 
                       predicciones$Target == target)
        
        if(en_gold && en_pred) {
          TP <- TP + 1
        } else if(!en_gold && en_pred) {
          FP <- FP + 1
        } else if(en_gold && !en_pred) {
          FN <- FN + 1
        } else {
          TN <- TN + 1
        }
      }
    }
    
    precision <- if((TP + FP) > 0) TP / (TP + FP) else 0
    recall <- if((TP + FN) > 0) TP / (TP + FN) else 0
    f1_score <- if((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    specificity <- if((TN + FP) > 0) TN / (TN + FP) else 0
    
    # 6. DYNAMIC ACCURACY (CORREGIDO)
    matriz_bin_raw <- read.table(ruta_matriz_bin,
                                  sep = "\t",
                                  header = TRUE,
                                  stringsAsFactors = FALSE)
    
    # La matriz viene con filas=timepoints, columnas=genes
    # Transponer para tener filas=genes, columnas=timepoints
    matriz_bin <- t(matriz_bin_raw)
    
    n_genes <- nrow(matriz_bin)
    n_timepoints <- ncol(matriz_bin)
    
    # Asignar nombres de genes como rownames
    rownames(matriz_bin) <- colnames(matriz_bin_raw)
    
    total_errors <- 0
    total_predictions <- 0
    
    for(t in seq_len(n_timepoints - 1)) {
      estado_actual <- as.numeric(matriz_bin[, t])
      names(estado_actual) <- rownames(matriz_bin)
      estado_real <- as.numeric(matriz_bin[, t + 1])
      
      for(i in seq_len(nrow(reglas))) {
        gen <- reglas$Gen[i]
        regla <- reglas$Regla[i]
        
        # Verificar que el gen existe en la matriz
        if(!(gen %in% rownames(matriz_bin))) {
          next
        }
        
        prediccion <- evaluate_boolean_rule(regla, estado_actual)
        
        if(!is.na(prediccion)) {
          gen_idx <- which(rownames(matriz_bin) == gen)
          error <- abs(estado_real[gen_idx] - prediccion)
          total_errors <- total_errors + error
          total_predictions <- total_predictions + 1
        }
      }
    }
    
    # Calcular Dynamic Accuracy usando predicciones realmente hechas
    dynacc_global <- if(total_predictions > 0) {
      1 - (total_errors / total_predictions)
    } else {
      0
    }
    
    # 7. PARSEAR INFORMACIÓN DEL MÉTODO
    method_info <- parse_method_name(method_name, size)
    bin_method <- if(!is.null(method_info)) method_info$bin_method else "UNKNOWN"
    rep_str <- if(!is.null(method_info)) method_info$rep else NA
    rep_numeric <- if(!is.null(method_info)) method_info$rep_numeric else NA
    
    # 8. RETORNAR RESULTADOS
    cat(sprintf("  Precision: %.4f | Recall: %.4f | F1: %.4f | DynAcc: %.4f\n", 
                precision, recall, f1_score, dynacc_global))
    
    return(data.frame(
      Size = paste0("size", size, "_", size_num),
      Size_Numeric = size_num,  # Para ordenamiento
      Method = method_name,
      Binarization = bin_method,
      Rep = rep_str,
      Rep_Numeric = rep_numeric,  # Para ordenamiento
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      Precision = precision,
      Recall = recall,
      F1_score = f1_score,
      Accuracy = accuracy,
      Specificity = specificity,
      DynamicAccuracy = dynacc_global,
      N_Predicciones = nrow(predicciones),
      N_GoldStandard = nrow(gold_standard_positive),
      Reguladores_Promedio = mean(n_reguladores),
      Reguladores_Mediana = median(n_reguladores),
      Reguladores_Max = max(n_reguladores),
      Score_Promedio = score_prom,
      MSE_Promedio = mse_prom,
      Correct_Promedio = correct_prom,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
}

# ==== ESCANEAR TODAS LAS CARPETAS Y ARCHIVOS ====
cat("Escaneando directorios...\n\n")

resultados_totales <- list()

# Iterar por cada carpeta size[SIZE]_X (1 al 5)
for(size_num in 1:5) {
  size_dir <- file.path(ruta_base_reglas, sprintf("size%d_%d", SIZE, size_num))
  
  if(!dir.exists(size_dir)) {
    cat(sprintf("ADVERTENCIA: No existe %s\n", size_dir))
    next
  }
  
  # Gold standard correspondiente a este size
  gold_file <- file.path(ruta_base_gold, 
                        sprintf("DREAM4_GoldStandard_InSilico_Size%d_%d.tsv", SIZE, size_num))
  
  if(!file.exists(gold_file)) {
    cat(sprintf("ADVERTENCIA: No existe gold standard %s\n", gold_file))
    next
  }
  
  cat(sprintf("=== Procesando size%d_%d ===\n", SIZE, size_num))
  
  # Obtener todas las subcarpetas de métodos
  method_dirs <- list.dirs(size_dir, full.names = TRUE, recursive = FALSE)
  
  for(method_dir in method_dirs) {
    method_name <- basename(method_dir)
    rules_file <- file.path(method_dir, "rules_by_gene.tsv")
    
    if(!file.exists(rules_file)) {
      cat(sprintf("  Saltando %s (no tiene rules_by_gene.tsv)\n", method_name))
      next
    }
    
    # Parsear información del método
    method_info <- parse_method_name(method_name, SIZE)
    
    if(is.null(method_info)) {
      cat(sprintf("  ADVERTENCIA: No se pudo parsear el nombre del método: %s\n", method_name))
      next
    }
    
    # Construir nombre del archivo binarizado
    bin_prefix <- if(method_info$bin_method == "SSD") "bin_SSD" else "bin_WCSS"
    
    # Detectar si el nombre de la carpeta tiene "_5times"
    has_5times_suffix <- grepl("_5times$", method_name)
    
    if(has_5times_suffix) {
      # Si tiene _5times, usar ese sufijo directamente en el nombre del archivo
      pattern_5times <- "_(\\d+)_5times$"
      matches <- regexec(pattern_5times, method_name)
      match_data <- regmatches(method_name, matches)[[1]]
      
      if(length(match_data) >= 2) {
        threshold_base <- match_data[2]
        threshold_str <- paste0(threshold_base, "_5times")
      } else {
        # Fallback
        threshold_str <- paste0(method_info$rep_numeric, "_5times")
      }
    } else {
      # Formatear threshold normalmente
      threshold_str <- method_info$rep
    }
    
    bin_file <- file.path(ruta_base_bin,
                         sprintf("size%d_%d", SIZE, size_num),
                         sprintf("%s_size%d_%d_%s.tsv", 
                                bin_prefix, 
                                SIZE,
                                method_info$network, 
                                threshold_str))
    
    if(!file.exists(bin_file)) {
      cat(sprintf("  ADVERTENCIA: No existe matriz binarizada %s\n", bin_file))
      next
    }
    
    # Evaluar esta red
    resultado <- evaluar_red(rules_file, SIZE, size_num, method_name, gold_file, bin_file)
    
    if(!is.null(resultado)) {
      resultados_totales[[length(resultados_totales) + 1]] <- resultado
    }
  }
  
  cat("\n")  # Separador entre sizes
}

# ==== CONSOLIDAR RESULTADOS ====
cat("\n========== CONSOLIDANDO RESULTADOS ==========\n")

if(length(resultados_totales) == 0) {
  cat("ERROR: No se pudo procesar ningún archivo\n")
  quit(status = 1)
}

resultados_df <- do.call(rbind, resultados_totales)

# Ordenar por Size y luego por Rep (numérico)
resultados_df <- resultados_df %>%
  arrange(Size_Numeric, Rep_Numeric) %>%
  select(-Size_Numeric, -Rep_Numeric)  # Eliminar columnas auxiliares de ordenamiento

cat(sprintf("Total de evaluaciones: %d\n", nrow(resultados_df)))

# Guardar tabla completa
output_file_complete <- file.path(output_dir, sprintf("BNI3_size%d_complete_evaluation_results.tsv", SIZE))
write.table(resultados_df,
            output_file_complete,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat(sprintf("\nResultados completos guardados en:\n  %s\n", output_file_complete))

# ==== MOSTRAR VISTA PREVIA ====
cat("\n========== VISTA PREVIA DE RESULTADOS ==========\n")
cat("\nPrimeras 15 filas:\n")
print(head(resultados_df, 15))

cat("\n========== PROCESO COMPLETADO ==========\n")
cat(sprintf("Resultados guardados en: %s\n", output_dir))
cat(sprintf("Archivo generado: BNI3_size%d_complete_evaluation_results.tsv\n", SIZE))
cat(sprintf("Total de evaluaciones: %d\n", nrow(resultados_df)))