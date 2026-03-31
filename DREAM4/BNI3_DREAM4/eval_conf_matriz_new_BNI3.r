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

# Función para calcular métricas estructurales y dinámicas de UNA regla específica
calcular_metricas_regla <- function(regla_text, gen_name, matriz_bin, gold_standard_positive, all_genes) {
  
  # 1. EXTRAER REGULADORES DE ESTA REGLA
  reguladores <- extract_regulators(regla_text, gen_name)
  
  # 2. CONSTRUIR PREDICCIONES PARA ESTA REGLA
  predicciones <- data.frame()
  if(length(reguladores) > 0) {
    for(reg in reguladores) {
      predicciones <- rbind(predicciones, 
                           data.frame(Regulador = reg, 
                                     Target = gen_name,
                                     stringsAsFactors = FALSE))
    }
  }
  
  # 3. CALCULAR MÉTRICAS ESTRUCTURALES (solo para este gen como target)
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  
  for(regulador in all_genes) {
    if(regulador == gen_name) next
    
    en_gold <- any(gold_standard_positive$Regulador == regulador & 
                   gold_standard_positive$Target == gen_name)
    
    en_pred <- any(predicciones$Regulador == regulador & 
                   predicciones$Target == gen_name)
    
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
  
  precision <- if((TP + FP) > 0) TP / (TP + FP) else 0
  recall <- if((TP + FN) > 0) TP / (TP + FN) else 0
  f1_score <- if((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  # 4. CALCULAR DYNAMIC ACCURACY
  n_timepoints <- ncol(matriz_bin)
  total_errors <- 0
  total_predictions <- 0
  
  for(t in seq_len(n_timepoints - 1)) {
    estado_actual <- as.numeric(matriz_bin[, t])
    names(estado_actual) <- rownames(matriz_bin)
    estado_real <- as.numeric(matriz_bin[, t + 1])
    
    if(!(gen_name %in% rownames(matriz_bin))) {
      next
    }
    
    prediccion <- evaluate_boolean_rule(regla_text, estado_actual)
    
    if(!is.na(prediccion)) {
      gen_idx <- which(rownames(matriz_bin) == gen_name)
      error <- abs(estado_real[gen_idx] - prediccion)
      total_errors <- total_errors + error
      total_predictions <- total_predictions + 1
    }
  }
  
  dynacc <- if(total_predictions > 0) {
    1 - (total_errors / total_predictions)
  } else {
    NA
  }
  
  return(list(
    # Métricas estructurales
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    accuracy = accuracy,
    # Métricas dinámicas
    dynacc = dynacc,
    n_predictions = total_predictions,
    n_errors = total_errors,
    # Info adicional
    n_regulators = length(reguladores)
  ))
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

# ==== FUNCIÓN PRINCIPAL DE EVALUACIÓN (MODIFICADA) ====
evaluar_red <- function(ruta_reglas, size, size_num, method_name, ruta_gold_standard, ruta_matriz_bin) {
  
  cat(sprintf("\n--- Procesando: size%d_%d - %s ---\n", size, size_num, method_name))
  
  # Verificar que existe el archivo
  if(!file.exists(ruta_reglas)) {
    cat(sprintf("  ADVERTENCIA: No se encontró %s\n", ruta_reglas))
    return(NULL)
  }
  
  tryCatch({
    # 1. LEER TODAS LAS REGLAS (SIN FILTRAR POR POSITION)
    reglas_completas <- read.table(ruta_reglas,
                                   sep = "\t",
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
    
    # 2. LEER MATRIZ BINARIZADA (necesaria para evaluaciones)
    matriz_bin_raw <- read.table(ruta_matriz_bin,
                                  sep = "\t",
                                  header = TRUE,
                                  stringsAsFactors = FALSE)
    
    # Transponer para tener filas=genes, columnas=timepoints
    matriz_bin <- t(matriz_bin_raw)
    rownames(matriz_bin) <- colnames(matriz_bin_raw)
    
    # 3. LEER GOLD STANDARD (necesario para evaluar cada regla)
    gold_standard <- read.table(ruta_gold_standard,
                                sep = "\t",
                                header = FALSE,
                                col.names = c("Regulador", "Target", "Interaccion"),
                                stringsAsFactors = FALSE)
    
    gold_standard_positive <- gold_standard %>%
      filter(Interaccion == 1) %>%
      select(Regulador, Target)
    
    all_genes <- unique(c(gold_standard$Regulador, gold_standard$Target))
    
    # 4. PARA CADA GEN, ENCONTRAR REGLAS CON SCORE MÁXIMO Y EVALUAR TODAS
    genes_unicos <- unique(reglas_completas$Gene)
    
    reglas_evaluadas_detalle <- list()
    reglas_seleccionadas <- data.frame()
    
    for(gen in genes_unicos) {
      reglas_gen <- reglas_completas %>% filter(Gene == gen)
      
      # Encontrar el Score máximo para este gen
      max_score <- max(reglas_gen$Score, na.rm = TRUE)
      
      # Filtrar solo las reglas con Score máximo
      reglas_max_score <- reglas_gen %>% filter(Score == max_score)
      
      # Evaluar todas las métricas para cada regla con Score máximo
      resultados_gen <- data.frame()
      
      for(i in seq_len(nrow(reglas_max_score))) {
        regla_info <- reglas_max_score[i, ]
        
        metricas <- calcular_metricas_regla(
          regla_info$Rule,
          gen,
          matriz_bin,
          gold_standard_positive,
          all_genes
        )
        
        # Parsear información del método para el detalle
        method_info_detail <- parse_method_name(method_name, size)
        
        # Almacenar resultado con información del archivo fuente
        resultado_temp <- data.frame(
          Size = size,
          Network = sprintf("size%d_%d", size, size_num),
          Method = method_name,
          Binarization = if(!is.null(method_info_detail)) method_info_detail$bin_method else "UNKNOWN",
          Rep = if(!is.null(method_info_detail)) method_info_detail$rep else NA,
          Source_File = basename(ruta_reglas),
          Gene = gen,
          Position = regla_info$Position,
          Rule = regla_info$Rule,
          Score = regla_info$Score,
          MSE = ifelse("MSE" %in% colnames(regla_info), regla_info$MSE, NA),
          Correct = regla_info$Correct,
          TP = metricas$TP,
          FP = metricas$FP,
          FN = metricas$FN,
          TN = metricas$TN,
          Precision = metricas$precision,
          Recall = metricas$recall,
          F1_score = metricas$f1_score,
          Accuracy = metricas$accuracy,
          DynamicAccuracy = metricas$dynacc,
          N_Predictions = metricas$n_predictions,
          N_Errors = metricas$n_errors,
          N_Regulators = metricas$n_regulators,
          stringsAsFactors = FALSE
        )
        
        resultados_gen <- rbind(resultados_gen, resultado_temp)
      }
      
      # Guardar detalle de todas las reglas evaluadas
      reglas_evaluadas_detalle[[gen]] <- resultados_gen
      
      # Seleccionar la regla con mejor Accuracy estructural
      if(nrow(resultados_gen) > 0) {
        # Ordenar por Accuracy (de mayor a menor) y tomar la primera
        mejor_regla <- resultados_gen %>%
          arrange(desc(Accuracy)) %>%
          slice(1)
        
        reglas_seleccionadas <- rbind(reglas_seleccionadas, mejor_regla)
      }
    }
    
    # 4. GUARDAR RESUMEN DETALLADO POR GEN
    resumen_file <- file.path(
      dirname(ruta_reglas),
      sprintf("detailed_rules_evaluation_size%d_%d_%s.tsv", size, size_num, method_name)
    )
    
    # Consolidar todas las evaluaciones
    todas_evaluaciones <- do.call(rbind, reglas_evaluadas_detalle)
    write.table(todas_evaluaciones,
                resumen_file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    
    cat(sprintf("  [DETALLE] Guardado resumen de %d reglas evaluadas en:\n", nrow(todas_evaluaciones)))
    cat(sprintf("            %s\n", basename(resumen_file)))
    
    # 5. CONSTRUIR TABLA DE REGLAS SELECCIONADAS (formato antiguo)
    reglas <- reglas_seleccionadas %>%
      select(Gene, Rule) %>%
      rename(Gen = Gene, Regla = Rule)
    
    # Extraer estadísticas de las reglas seleccionadas
    score_prom <- mean(reglas_seleccionadas$Score, na.rm = TRUE)
    
    mse_prom <- NA
    if("MSE" %in% colnames(reglas_seleccionadas)) {
      if(is.numeric(reglas_seleccionadas$MSE)) {
        mse_prom <- mean(reglas_seleccionadas$MSE, na.rm = TRUE)
      }
    }
    
    correct_prom <- mean(reglas_seleccionadas$Correct, na.rm = TRUE)
    accuracy_prom <- mean(reglas_seleccionadas$Accuracy, na.rm = TRUE)
    dynacc_prom <- mean(reglas_seleccionadas$DynamicAccuracy, na.rm = TRUE)
    
    # 6. EXTRAER REGULADORES DE LAS REGLAS SELECCIONADAS
    reglas$Reguladores <- lapply(seq_len(nrow(reglas)), function(i) {
      extract_regulators(reglas$Regla[i], reglas$Gen[i])
    })
    
    n_reguladores <- sapply(reglas$Reguladores, length)
    
    # 7. CONSTRUIR PREDICCIONES (ya no necesitamos leer gold standard de nuevo)
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
    
    # 8. CALCULAR MÉTRICAS ESTRUCTURALES GLOBALES (toda la red)
    # Aunque ya calculamos por gen, necesitamos la evaluación global de la red completa
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
    
    # 9. PARSEAR INFORMACIÓN DEL MÉTODO
    method_info <- parse_method_name(method_name, size)
    bin_method <- if(!is.null(method_info)) method_info$bin_method else "UNKNOWN"
    rep_str <- if(!is.null(method_info)) method_info$rep else NA
    rep_numeric <- if(!is.null(method_info)) method_info$rep_numeric else NA
    
    # 10. GENERAR RESUMEN DESPLEGABLE
    cat("\n  ╔═══════════════════════════════════════════════════════════════╗\n")
    cat("  ║           RESUMEN DE SELECCIÓN DE REGLAS                     ║\n")
    cat("  ╚═══════════════════════════════════════════════════════════════╝\n")
    cat(sprintf("  Total de genes:                    %d\n", length(genes_unicos)))
    cat(sprintf("  Reglas con Score máximo evaluadas: %d\n", nrow(todas_evaluaciones)))
    cat(sprintf("  Reglas seleccionadas (mejores):    %d\n", nrow(reglas_seleccionadas)))
    cat(sprintf("  Accuracy promedio (por gen):       %.4f\n", accuracy_prom))
    cat(sprintf("  DynamicAccuracy promedio:          %.4f\n", dynacc_prom))
    
    # Mostrar algunos ejemplos de genes con múltiples reglas
    genes_con_multiples <- sapply(reglas_evaluadas_detalle, nrow)
    genes_con_multiples <- genes_con_multiples[genes_con_multiples > 1]
    
    if(length(genes_con_multiples) > 0) {
      cat(sprintf("\n  Genes con múltiples reglas (Score máximo): %d\n", length(genes_con_multiples)))
      cat("  Ejemplos (ordenados por Accuracy):\n")
      
      for(i in seq_len(min(5, length(genes_con_multiples)))) {
        gen_ejemplo <- names(genes_con_multiples)[i]
        n_reglas <- genes_con_multiples[i]
        # ORDENAR por Accuracy descendente antes de mostrar
        reglas_ejemplo <- reglas_evaluadas_detalle[[gen_ejemplo]] %>%
          arrange(desc(Accuracy))
        
        cat(sprintf("\n    - %s: %d reglas evaluadas\n", gen_ejemplo, n_reglas))
        for(j in seq_len(nrow(reglas_ejemplo))) {
          # La primera fila (j==1) después de ordenar es la seleccionada
          seleccionada <- if(j == 1) " ← SELECCIONADA (mejor Acc)" else ""
          cat(sprintf("      Pos %d: Acc=%.4f, DynAcc=%.4f, Score=%.4f%s\n",
                     reglas_ejemplo$Position[j],
                     reglas_ejemplo$Accuracy[j],
                     reglas_ejemplo$DynamicAccuracy[j],
                     reglas_ejemplo$Score[j],
                     seleccionada))
        }
      }
    }
    
    cat("\n  ───────────────────────────────────────────────────────────────\n\n")
    
    # 11. RETORNAR RESULTADOS
    cat(sprintf("  Precision: %.4f | Recall: %.4f | F1: %.4f | Acc: %.4f | DynAcc: %.4f\n", 
                precision, recall, f1_score, accuracy, dynacc_prom))
    
    return(list(
      resultado_final = data.frame(
        Size = paste0("size", size, "_", size_num),
        Size_Numeric = size_num,
        Method = method_name,
        Binarization = bin_method,
        Rep = rep_str,
        Rep_Numeric = rep_numeric,
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        Precision = precision,
        Recall = recall,
        F1_score = f1_score,
        Accuracy = accuracy,
        Specificity = specificity,
        DynamicAccuracy = dynacc_prom,
        N_Predicciones = nrow(predicciones),
        N_GoldStandard = nrow(gold_standard_positive),
        Reguladores_Promedio = mean(n_reguladores),
        Reguladores_Mediana = median(n_reguladores),
        Reguladores_Max = max(n_reguladores),
        Score_Promedio = score_prom,
        MSE_Promedio = mse_prom,
        Correct_Promedio = correct_prom,
        N_Genes = length(genes_unicos),
        N_Reglas_Evaluadas = nrow(todas_evaluaciones),
        N_Genes_Multiple_Rules = length(genes_con_multiples),
        stringsAsFactors = FALSE
      ),
      resumen_detallado = todas_evaluaciones
    ))
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
}

# ==== ESCANEAR TODAS LAS CARPETAS Y ARCHIVOS ====
cat("Escaneando directorios...\n\n")

resultados_totales <- list()
resumenes_detallados <- list()

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
      pattern_5times <- "_(\\d+)_5times$"
      matches <- regexec(pattern_5times, method_name)
      match_data <- regmatches(method_name, matches)[[1]]
      
      if(length(match_data) >= 2) {
        threshold_base <- match_data[2]
        threshold_str <- paste0(threshold_base, "_5times")
      } else {
        threshold_str <- paste0(method_info$rep_numeric, "_5times")
      }
    } else {
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
      resultados_totales[[length(resultados_totales) + 1]] <- resultado$resultado_final
      
      # Guardar resumen detallado con identificador
      id_resumen <- sprintf("size%d_%d_%s", SIZE, size_num, method_name)
      resumenes_detallados[[id_resumen]] <- resultado$resumen_detallado
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
  select(-Size_Numeric, -Rep_Numeric)

cat(sprintf("Total de evaluaciones: %d\n", nrow(resultados_df)))

# Guardar tabla completa
output_file_complete <- file.path(output_dir, sprintf("BNI3_size%d_complete_evaluation_results_BEST_ACCURACY.tsv", SIZE))
write.table(resultados_df,
            output_file_complete,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat(sprintf("\nResultados completos guardados en:\n  %s\n", output_file_complete))

# ==== GUARDAR RESUMEN CONSOLIDADO DE TODAS LAS REGLAS EVALUADAS ====
if(length(resumenes_detallados) > 0) {
  output_file_detailed <- file.path(output_dir, sprintf("BNI3_size%d_ALL_RULES_DETAILED.tsv", SIZE))
  
  # Consolidar todos los resúmenes (ya tienen las columnas identificadoras)
  todos_resumenes <- do.call(rbind, resumenes_detallados)
  
  write.table(todos_resumenes,
              output_file_detailed,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  cat(sprintf("\nResumen detallado de TODAS las reglas evaluadas guardado en:\n  %s\n", output_file_detailed))
  cat(sprintf("Total de reglas evaluadas: %d\n", nrow(todos_resumenes)))
}

# ==== MOSTRAR VISTA PREVIA ====
cat("\n========== VISTA PREVIA DE RESULTADOS ==========\n")
cat("\nPrimeras 15 filas:\n")
print(head(resultados_df, 15))

# ==== ESTADÍSTICAS FINALES ====
cat("\n========== ESTADÍSTICAS FINALES ==========\n")
cat(sprintf("Total de redes evaluadas:          %d\n", nrow(resultados_df)))
cat(sprintf("Total de reglas evaluadas:         %d\n", sum(resultados_df$N_Reglas_Evaluadas)))
cat(sprintf("Promedio de reglas por red:        %.2f\n", mean(resultados_df$N_Reglas_Evaluadas)))
cat(sprintf("Redes con genes con múltiples:     %d\n", sum(resultados_df$N_Genes_Multiple_Rules > 0)))

cat("\n========== PROCESO COMPLETADO ==========\n")
cat(sprintf("Resultados guardados en: %s\n", output_dir))
cat(sprintf("Archivo principal: BNI3_size%d_complete_evaluation_results_BEST_ACCURACY.tsv\n", SIZE))
cat(sprintf("Archivo detallado: BNI3_size%d_ALL_RULES_DETAILED.tsv\n", SIZE))
cat(sprintf("Total de evaluaciones: %d\n", nrow(resultados_df)))