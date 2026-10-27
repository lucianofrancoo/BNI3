library("dplyr")

# Structural evaluation of BoolNet networks against the DREAM4 size10 gold standard,
# with the same scoring as eval_gold_standard.r (self-regulation excluded).
# Usage: Rscript eval_gold_standard_boolnet.r [boolnet_dir] [gold_dir] [output_file]
#   Defaults evaluate the paper's results: ../results/BoolNet  ../data/gold_standard
#   ../results/BoolNet_size10_complete_evaluation_results.tsv

SIZE <- 10

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if(length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
dream4_dir <- dirname(script_dir)
args <- commandArgs(trailingOnly = TRUE)
ruta_base_reglas <- normalizePath(if(length(args) >= 1) args[1] else file.path(dream4_dir, "results", "BoolNet"))
ruta_base_bin <- ruta_base_reglas
ruta_base_gold <- normalizePath(if(length(args) >= 2) args[2] else file.path(dream4_dir, "data", "gold_standard"))
output_file_complete <- if(length(args) >= 3) args[3] else
  file.path(dream4_dir, "results", "BoolNet_size10_complete_evaluation_results.tsv")
output_dir <- dirname(output_file_complete)

cat(sprintf("========== CONSOLIDATED EVALUATION BOOLNET - DREAM4 SIZE%d ==========\n", SIZE))
cat(sprintf("Rules base dir:   %s\n", ruta_base_reglas))
cat(sprintf("Gold base dir:    %s\n", ruta_base_gold))
cat(sprintf("Binarized dir:    %s\n", ruta_base_bin))
cat(sprintf("Output dir:       %s\n\n", output_dir))

# ==== HELPER FUNCTIONS ====

extract_regulators <- function(rule_text, target_gene) {
  # Strip negations and operators from the rule
  cleaned <- gsub("!", " ", rule_text)
  cleaned <- gsub("\\(|\\)|&|\\|", " ", cleaned)
  genes <- unlist(strsplit(cleaned, "\\s+"))
  genes <- genes[grepl("^G[0-9]+$", genes)]
  genes <- unique(genes)
  genes <- genes[genes != target_gene]
  return(genes)
}

evaluate_boolean_rule <- function(rule_text, state_vector) {
  expression <- rule_text
  
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

# Parse the file name
parse_file_name <- function(file_name, size) {
  # Format: reglas_Boolnet_bestfit_sizeX_Y_Z.tsv or sizeX_Y_Z_5times.tsv
  pattern <- paste0("size", size, "_(\\d+)_(\\d+)(_5times)?")
  matches <- regexec(pattern, file_name)
  match_data <- regmatches(file_name, matches)[[1]]
  
  if(length(match_data) >= 3) {
    size_num <- as.integer(match_data[2])
    rep_base <- as.integer(match_data[3])
    has_5times <- length(match_data) >= 4 && match_data[4] == "_5times"
    
    # Build the Rep string
    if(has_5times) {
      rep_str <- paste0(rep_base, "_5times")
    } else {
      rep_str <- as.character(rep_base)
    }
    
    return(list(
      size_num = size_num,
      rep = rep_str,
      rep_numeric = rep_base,
      has_5times = has_5times
    ))
  }
  
  return(NULL)
}

# ==== MAIN EVALUATION FUNCTION ====
evaluar_red <- function(ruta_reglas, size_num, file_name, ruta_gold_standard, ruta_matriz_bin, size) {
  
  cat(sprintf("\n--- Processing: size%d_%d - %s ---\n", size, size_num, file_name))
  
  if(!file.exists(ruta_reglas)) {
    cat(sprintf("  WARNING: not found: %s\n", ruta_reglas))
    return(NULL)
  }
  
  tryCatch({
    # 1. READ RULES
    reglas_completas <- read.table(ruta_reglas,
                                   sep = "\t",
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
    
    reglas <- reglas_completas %>%
      select(Gen, Regla, Error, N_reguladores)
    
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
    matriz_bin <- read.table(ruta_matriz_bin,
                             sep = "\t",
                             header = TRUE,
                             stringsAsFactors = FALSE)
    
    rownames(matriz_bin) <- matriz_bin$Gen
    matriz_bin <- matriz_bin[, -1]
    
    n_genes <- nrow(matriz_bin)
    n_timepoints <- ncol(matriz_bin)
    
    total_errors <- 0
    
    for(t in seq_len(n_timepoints - 1)) {
      estado_actual <- as.numeric(matriz_bin[, t])
      names(estado_actual) <- rownames(matriz_bin)
      estado_real <- as.numeric(matriz_bin[, t + 1])
      
      for(i in seq_len(nrow(reglas))) {
        gen <- reglas$Gen[i]
        regla <- reglas$Regla[i]
        prediccion <- evaluate_boolean_rule(regla, estado_actual)
        
        if(!is.na(prediccion)) {
          gen_idx <- which(rownames(matriz_bin) == gen)
          error <- abs(estado_real[gen_idx] - prediccion)
          total_errors <- total_errors + error
        }
      }
    }
    
    dynacc_global <- 1 - (total_errors / (n_genes * (n_timepoints - 1)))
    
    # 7. PARSE THE FILE NAME
    file_info <- parse_file_name(file_name, size)
    rep_str <- if(!is.null(file_info)) file_info$rep else NA
    rep_numeric <- if(!is.null(file_info)) file_info$rep_numeric else NA
    
    # 8. RETURN RESULTS
    cat(sprintf("  Precision: %.4f | Recall: %.4f | F1: %.4f | MCC: %.4f | DynAcc: %.4f\n", 
                precision, recall, f1_score, mcc, dynacc_global))
    
    return(data.frame(
      Size = paste0("size", size, "_", size_num),
      Size_Numeric = size_num,
      Method = "BoolNet_bestfit",
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
      MCC = mcc,
      Balanced_Accuracy = balanced_accuracy,
      Coverage_Ratio = coverage_ratio,
      DynamicAccuracy = dynacc_global,
      N_Predicciones = nrow(predicciones),
      N_GoldStandard = nrow(gold_standard_positive),
      Reguladores_Promedio = mean(n_reguladores),
      Reguladores_Mediana = median(n_reguladores),
      Reguladores_Max = max(n_reguladores),
      Error_Promedio = mean(reglas$Error, na.rm = TRUE),
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

# Loop over each sizeX_Y folder (1 to 5)
for(size_num in 1:5) {
  size_dir <- file.path(ruta_base_reglas, paste0("size", SIZE, "_", size_num))
  
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
  
  # All rule files in this folder
  rules_files <- list.files(size_dir, 
                            pattern = "^reglas_BoolNet_bestfit_.*\\.tsv$",
                            full.names = TRUE)
  
  for(rules_file in rules_files) {
    file_name <- basename(rules_file)
    
    # Parse the file name
    file_info <- parse_file_name(file_name, SIZE)
    
    if(is.null(file_info)) {
      cat(sprintf("  WARNING: could not parse the file name: %s\n", file_name))
      next
    }
    
    # Name of the matching binarized file
    bin_name <- gsub("^reglas_BoolNet_bestfit_", "bin_", file_name)
    bin_file <- file.path(size_dir, bin_name)
    
    if(!file.exists(bin_file)) {
      cat(sprintf("  WARNING: binarized matrix does not exist: %s\n", bin_name))
      next
    }
    
    # Evaluate this network
    resultado <- evaluar_red(rules_file, size_num, file_name, gold_file, bin_file, SIZE)
    
    if(!is.null(resultado)) {
      resultados_totales[[length(resultados_totales) + 1]] <- resultado
    }
  }
  
  cat("\n")
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
  select(-Size_Numeric, -Rep_Numeric)

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
cat("\n========== SUMMARY BY SAMPLING ==========\n")

# Split full datasets vs 5times
resultados_df$Tipo_Muestreo <- ifelse(grepl("_5times$", resultados_df$Rep), 
                                      "5_timepoints", 
                                      "21_timepoints")

resumen <- resultados_df %>%
  group_by(Tipo_Muestreo) %>%
  summarise(
    N_casos = n(),
    F1_medio = mean(F1_score, na.rm = TRUE),
    MCC_medio = mean(MCC, na.rm = TRUE),
    Coverage_medio = mean(Coverage_Ratio, na.rm = TRUE),
    Accuracy_medio = mean(Accuracy, na.rm = TRUE),
    Balanced_Acc_medio = mean(Balanced_Accuracy, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nComparison 21 vs 5 timepoints:\n")
print(resumen)

cat("\n========== DONE ==========\n")
cat(sprintf("Results saved in: %s\n", output_dir))
cat(sprintf("Output file: BoolNet_complete_evaluation_results.tsv\n"))
cat(sprintf("Total evaluations: %d\n", nrow(resultados_df)))