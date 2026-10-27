library("dplyr")

# Structural evaluation of BNI3 networks against the DREAM4 size10 gold standard.
#
# Usage: Rscript eval_gold_standard.r [rules_dir] [binarized_dir] [gold_dir] [output_file]
#   Defaults (relative to this script) evaluate the results used in the paper:
#     ../results/BNI3/rules   ../results/BNI3/binarized   ../data/gold_standard
#     ../results/BNI3_size10_complete_evaluation_results.tsv
#
# Expected layout:
#   <binarized_dir>/size10_X/<series>/<series>_binarized_<SSD|WCSS>.tsv   (BNI3_Binarizator.py)
#   <rules_dir>/size10_X/rules_<SSD|WCSS>_<series>/rules_by_gene_evaluated.tsv
#     (falls back to rules_by_gene.tsv when the evaluation was skipped)
#
# Self-regulation is excluded on both sides: the target gene is dropped from its own
# regulators and the diagonal is skipped, so 90 ordered pairs are scored per network.

SIZE <- 10

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if(length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
dream4_dir <- dirname(script_dir)

ruta_base_reglas <- normalizePath(if(length(args) >= 1) args[1] else file.path(dream4_dir, "results", "BNI3", "rules"))
ruta_base_bin <- normalizePath(if(length(args) >= 2) args[2] else file.path(dream4_dir, "results", "BNI3", "binarized"))
ruta_base_gold <- normalizePath(if(length(args) >= 3) args[3] else file.path(dream4_dir, "data", "gold_standard"))
output_file_complete <- if(length(args) >= 4) args[4] else
  file.path(dream4_dir, "results", "BNI3_size10_complete_evaluation_results.tsv")
output_dir <- dirname(output_file_complete)

cat(sprintf("========== CONSOLIDATED EVALUATION BNI3 - DREAM4 SIZE%d ==========\n", SIZE))
cat(sprintf("Rules base dir:   %s\n", ruta_base_reglas))
cat(sprintf("Gold base dir:    %s\n", ruta_base_gold))
cat(sprintf("Binarized dir:    %s\n", ruta_base_bin))
cat(sprintf("Output dir:       %s\n\n", output_dir))

# ==== HELPER FUNCTIONS ====

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

# Matthews correlation coefficient
calculate_mcc <- function(TP, FP, FN, TN) {
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  
  if(denominator == 0) {
    return(0)
  }
  
  return(numerator / denominator)
}

# Parse the method folder name
parse_method_name <- function(method_name, size) {
  # Binarization method (SSD or WCSS)
  bin_method <- ifelse(grepl("SSD", method_name), "SSD", "WCSS")
  
  # Generic pattern that adapts to the network size
  # Format: rules_[METHOD]_size[SIZE]_[NETWORK]_[REP]
  pattern <- sprintf("rules_(SSD|WCSS)_size%d_(\\d+)_(\\d+)(_5times)?$", size)
  matches <- regexec(pattern, method_name)
  match_data <- regmatches(method_name, matches)[[1]]
  
  if(length(match_data) >= 4) {
    network <- as.integer(match_data[3])
    rep_base <- as.integer(match_data[4])
    has_5times <- length(match_data) >= 5 && match_data[5] == "_5times"
    
    # Build the Rep string
    if(has_5times) {
      rep_str <- paste0(rep_base, "_5times")
    } else {
      rep_str <- as.character(rep_base)
    }
    
    return(list(
      bin_method = bin_method,
      network = network,
      rep = rep_str,
      rep_numeric = rep_base  # For sorting
    ))
  }
  
  return(NULL)
}

# ==== MAIN EVALUATION FUNCTION ====
evaluar_red <- function(ruta_reglas, size, size_num, method_name, ruta_gold_standard, ruta_matriz_bin) {
  
  cat(sprintf("\n--- Processing: size%d_%d - %s ---\n", size, size_num, method_name))
  
  # Check that the file exists
  if(!file.exists(ruta_reglas)) {
    cat(sprintf("  WARNING: not found: %s\n", ruta_reglas))
    return(NULL)
  }
  
  tryCatch({
    # 1. READ RULES
    reglas_completas <- read.table(ruta_reglas,
                                   sep = "\t",
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   na.strings = c("NA", "N/A", ""))
    
    # Take the first rule of each gene (already sorted by descending Score)
    reglas <- reglas_completas %>%
      group_by(Gene) %>%
      slice(1) %>%
      ungroup() %>%
      select(Gene, Rule) %>%
      rename(Gen = Gene, Regla = Rule)
    
    # Statistics of the selected rules
    score_prom <- NA
    mse_prom <- NA
    correct_prom <- NA
    
    if("Score" %in% colnames(reglas_completas)) {
      reglas_stats <- reglas_completas %>% 
        group_by(Gene) %>%
        slice(1) %>%
        ungroup()
      
      score_prom <- mean(reglas_stats$Score, na.rm = TRUE)
      
      # Special handling for MSE
      if("MSE" %in% colnames(reglas_stats)) {
        mse_prom <- mean(reglas_stats$MSE, na.rm = TRUE)
      }
      
      correct_prom <- mean(reglas_stats$Correct, na.rm = TRUE)
    }
    
    # 2. EXTRACT REGULATORS
    reglas$Reguladores <- lapply(seq_len(nrow(reglas)), function(i) {
      extract_regulators(reglas$Regla[i], reglas$Gen[i])
    })
    
    n_reguladores <- sapply(reglas$Reguladores, length)
    
    # 3. READ GOLD STANDARD
    gold_standard <- read.table(ruta_gold_standard,
                                sep = "\t",
                                header = FALSE,
                                col.names = c("Regulador", "Target", "Interaccion"),
                                stringsAsFactors = FALSE)
    
    gold_standard_positive <- gold_standard %>%
      filter(Interaccion == 1) %>%
      select(Regulador, Target)
    
    # 4. BUILD PREDICTIONS
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
    
    # 5. STRUCTURAL METRICS
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
    
    # ADDITIONAL METRICS
    # Matthews Correlation Coefficient
    mcc <- calculate_mcc(TP, FP, FN, TN)
    
    # Balanced Accuracy
    balanced_accuracy <- (recall + specificity) / 2
    
    # Coverage Ratio (predicted edges relative to the gold standard)
    coverage_ratio <- if(nrow(gold_standard_positive) > 0) {
      nrow(predicciones) / nrow(gold_standard_positive)
    } else {
      0
    }
    
    # 6. DYNAMIC ACCURACY
    matriz_bin_raw <- read.table(ruta_matriz_bin,
                                  sep = "\t",
                                  header = TRUE,
                                  stringsAsFactors = FALSE)
    
    # The matrix comes with rows=timepoints, columns=genes
    # Transpose to rows=genes, columns=timepoints
    matriz_bin <- t(matriz_bin_raw)
    
    n_genes <- nrow(matriz_bin)
    n_timepoints <- ncol(matriz_bin)
    
    # Gene names as rownames
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
        
        # Check that the gene is in the matrix
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
    
    # Dynamic Accuracy over the predictions actually made
    dynacc_global <- if(total_predictions > 0) {
      1 - (total_errors / total_predictions)
    } else {
      0
    }
    
    # 7. PARSE THE METHOD NAME
    method_info <- parse_method_name(method_name, size)
    bin_method <- if(!is.null(method_info)) method_info$bin_method else "UNKNOWN"
    rep_str <- if(!is.null(method_info)) method_info$rep else NA
    rep_numeric <- if(!is.null(method_info)) method_info$rep_numeric else NA
    
    # 8. RETURN RESULTS
    cat(sprintf("  Precision: %.4f | Recall: %.4f | F1: %.4f | MCC: %.4f | DynAcc: %.4f\n", 
                precision, recall, f1_score, mcc, dynacc_global))
    
    return(data.frame(
      Size = paste0("size", size, "_", size_num),
      Size_Numeric = size_num,  # For sorting
      Method = method_name,
      Binarization = bin_method,
      Rep = rep_str,
      Rep_Numeric = rep_numeric,  # For sorting
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      Precision = precision,
      Recall = recall,
      F1_score = f1_score,
      Accuracy = accuracy,
      Specificity = specificity,
      MCC = mcc,
      Balanced_Accuracy = balanced_accuracy,
      Coverage_Ratio = coverage_ratio,
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

# ==== SCAN ALL FOLDERS AND FILES ====
cat("Scanning directories...\n\n")

resultados_totales <- list()

# Loop over each size[SIZE]_X folder (1 to 5)
for(size_num in 1:5) {
  size_dir <- file.path(ruta_base_reglas, sprintf("size%d_%d", SIZE, size_num))
  
  if(!dir.exists(size_dir)) {
    cat(sprintf("WARNING: does not exist: %s\n", size_dir))
    next
  }
  
  # Gold standard for this network
  gold_file <- file.path(ruta_base_gold, 
                        sprintf("DREAM4_GoldStandard_InSilico_Size%d_%d.tsv", SIZE, size_num))
  
  if(!file.exists(gold_file)) {
    cat(sprintf("WARNING: gold standard does not exist: %s\n", gold_file))
    next
  }
  
  cat(sprintf("=== Processing size%d_%d ===\n", SIZE, size_num))
  
  # All method subfolders
  method_dirs <- list.dirs(size_dir, full.names = TRUE, recursive = FALSE)
  
  for(method_dir in method_dirs) {
    method_name <- basename(method_dir)
    rules_file <- file.path(method_dir, "rules_by_gene_evaluated.tsv")
    if(!file.exists(rules_file)) {
      rules_file <- file.path(method_dir, "rules_by_gene.tsv")
    }
    
    if(!file.exists(rules_file)) {
      cat(sprintf("  Skipping %s (no rules_by_gene[_evaluated].tsv)\n", method_name))
      next
    }
    
    # Parse the method name
    method_info <- parse_method_name(method_name, SIZE)
    
    if(is.null(method_info)) {
      cat(sprintf("  WARNING: could not parse the method name: %s\n", method_name))
      next
    }
    
    # Binarized matrix written by BNI3_Binarizator.py:
    #   bin_Data/sizeN_X/<stem>/<stem>_binarized_<METHOD>.tsv
    stem <- sub("^rules_(SSD|WCSS)_", "", method_name)
    bin_file <- file.path(ruta_base_bin,
                         sprintf("size%d_%d", SIZE, size_num),
                         stem,
                         sprintf("%s_binarized_%s.tsv", stem, method_info$bin_method))
    
    if(!file.exists(bin_file)) {
      cat(sprintf("  WARNING: binarized matrix does not exist: %s\n", bin_file))
      next
    }
    
    # Evaluate this network
    resultado <- evaluar_red(rules_file, SIZE, size_num, method_name, gold_file, bin_file)
    
    if(!is.null(resultado)) {
      resultados_totales[[length(resultados_totales) + 1]] <- resultado
    }
  }
  
  cat("\n")  # Separator between networks
}

# ==== CONSOLIDATE RESULTS ====
cat("\n========== CONSOLIDATING RESULTS ==========\n")

if(length(resultados_totales) == 0) {
  cat("ERROR: no file could be processed\n")
  quit(status = 1)
}

resultados_df <- do.call(rbind, resultados_totales)

# Sort by Size, then by Rep (numeric)
resultados_df <- resultados_df %>%
  arrange(Size_Numeric, Rep_Numeric) %>%
  select(-Size_Numeric, -Rep_Numeric)  # Drop the helper sorting columns

cat(sprintf("Total evaluations: %d\n", nrow(resultados_df)))

# Save the full table
write.table(resultados_df,
            output_file_complete,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat(sprintf("\nFull results saved in:\n  %s\n", output_file_complete))

# ==== PREVIEW ====
cat("\n========== RESULTS PREVIEW ==========\n")
cat("\nFirst 15 rows:\n")
print(head(resultados_df, 15))

# ==== KEY METRICS SUMMARY ====
cat("\n========== SUMMARY BY BINARIZATION METHOD AND SAMPLING ==========\n")

# Split full datasets vs 5times
resultados_df$Tipo_Muestreo <- ifelse(grepl("_5times$", resultados_df$Rep), 
                                      "5_timepoints", 
                                      "21_timepoints")

resumen <- resultados_df %>%
  group_by(Binarization, Tipo_Muestreo) %>%
  summarise(
    N_casos = n(),
    F1_medio = mean(F1_score, na.rm = TRUE),
    MCC_medio = mean(MCC, na.rm = TRUE),
    Coverage_medio = mean(Coverage_Ratio, na.rm = TRUE),
    Accuracy_medio = mean(Accuracy, na.rm = TRUE),
    Balanced_Acc_medio = mean(Balanced_Accuracy, na.rm = TRUE),
    DynAcc_medio = mean(DynamicAccuracy, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nComparison by binarization method and timepoints:\n")
print(resumen)

# Additional summary by sampling only (both methods pooled)
cat("\n========== OVERALL SUMMARY BY SAMPLING ==========\n")

resumen_global <- resultados_df %>%
  group_by(Tipo_Muestreo) %>%
  summarise(
    N_casos = n(),
    F1_medio = mean(F1_score, na.rm = TRUE),
    MCC_medio = mean(MCC, na.rm = TRUE),
    Coverage_medio = mean(Coverage_Ratio, na.rm = TRUE),
    Accuracy_medio = mean(Accuracy, na.rm = TRUE),
    Balanced_Acc_medio = mean(Balanced_Accuracy, na.rm = TRUE),
    DynAcc_medio = mean(DynamicAccuracy, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nComparison 21 vs 5 timepoints (all methods):\n")
print(resumen_global)

cat("\n========== DONE ==========\n")
cat(sprintf("Results saved in: %s\n", output_dir))
cat(sprintf("Output file: %s\n", basename(output_file_complete)))
cat(sprintf("Total evaluations: %d\n", nrow(resultados_df)))