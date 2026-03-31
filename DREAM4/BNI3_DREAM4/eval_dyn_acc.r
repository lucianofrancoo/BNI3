# ==== SCRIPT DE DIAGNÓSTICO PARA DYNAMIC ACCURACY ====

library(dplyr)

# Usar un caso real
ruta_reglas <- "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/rules/size100_1/rules_WCSS_size100_1_1/rules_by_gene.tsv"
ruta_matriz_bin <- "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/bin_Data/size100_1/bin_WCSS_size100_1_1.tsv"

# Función de evaluación
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

# Leer reglas
cat("Leyendo reglas...\n")
reglas_completas <- read.table(ruta_reglas,
                               sep = "\t",
                               header = TRUE,
                               stringsAsFactors = FALSE)

reglas <- reglas_completas %>%
  filter(Position == 1) %>%
  select(Gene, Rule, Correct) %>%
  rename(Gen = Gene, Regla = Rule)

cat(sprintf("Total genes con reglas: %d\n", nrow(reglas)))
cat(sprintf("Promedio Correct: %.2f\n", mean(reglas$Correct)))

# Leer matriz binarizada
cat("\nLeyendo matriz binarizada...\n")
matriz_bin <- read.table(ruta_matriz_bin,
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE)

# La matriz está en formato: filas = timepoints, columnas = genes
# Necesitamos transponer para tener: filas = genes, columnas = timepoints
cat("Dimensiones originales (antes de transponer):\n")
cat(sprintf("  Filas (timepoints): %d\n", nrow(matriz_bin)))
cat(sprintf("  Columnas (genes): %d\n", ncol(matriz_bin)))

# Transponer la matriz
matriz_bin <- t(matriz_bin)

n_genes <- nrow(matriz_bin)
n_timepoints <- ncol(matriz_bin)

cat(sprintf("Genes en matriz: %d\n", n_genes))
cat(sprintf("Timepoints: %d\n", n_timepoints))
cat(sprintf("Total transiciones a evaluar: %d\n", n_timepoints - 1))

# Verificar que los genes coinciden
genes_reglas <- reglas$Gen
genes_matriz <- rownames(matriz_bin)

cat(sprintf("\nGenes en reglas: %d\n", length(genes_reglas)))
cat(sprintf("Genes en matriz: %d\n", length(genes_matriz)))
cat(sprintf("Genes en común: %d\n", length(intersect(genes_reglas, genes_matriz))))

if(length(setdiff(genes_reglas, genes_matriz)) > 0) {
  cat("\nADVERTENCIA: Genes en reglas pero no en matriz:\n")
  print(setdiff(genes_reglas, genes_matriz))
}

if(length(setdiff(genes_matriz, genes_reglas)) > 0) {
  cat("\nADVERTENCIA: Genes en matriz pero no en reglas:\n")
  print(head(setdiff(genes_matriz, genes_reglas), 10))
}

# Calcular Dynamic Accuracy con diagnóstico detallado
cat("\n========== CALCULANDO DYNAMIC ACCURACY ==========\n")

total_errors <- 0
total_predictions <- 0
total_na <- 0
errors_por_gen <- rep(0, n_genes)
predictions_por_gen <- rep(0, n_genes)
names(errors_por_gen) <- rownames(matriz_bin)
names(predictions_por_gen) <- rownames(matriz_bin)

# Analizar solo las primeras 3 transiciones en detalle
cat("\n=== ANÁLISIS DETALLADO DE LAS PRIMERAS 3 TRANSICIONES ===\n")

for(t in 1:min(3, n_timepoints - 1)) {
  cat(sprintf("\n--- Transición %d -> %d ---\n", t, t + 1))
  
  estado_actual <- as.numeric(matriz_bin[, t])
  names(estado_actual) <- rownames(matriz_bin)
  
  estado_real <- as.numeric(matriz_bin[, t + 1])
  
  errors_this_t <- 0
  
  # Analizar solo los primeros 5 genes en detalle
  for(i in 1:min(5, nrow(reglas))) {
    gen <- reglas$Gen[i]
    regla <- reglas$Regla[i]
    
    # Verificar si el gen existe en la matriz
    if(!(gen %in% names(estado_actual))) {
      cat(sprintf("  Gen %s no está en la matriz\n", gen))
      next
    }
    
    prediccion <- evaluate_boolean_rule(regla, estado_actual)
    
    gen_idx <- which(rownames(matriz_bin) == gen)
    real <- estado_real[gen_idx]
    
    if(!is.na(prediccion)) {
      error <- abs(real - prediccion)
      errors_this_t <- errors_this_t + error
      
      cat(sprintf("  %s: Real=%d, Pred=%d, Error=%d | Regla: %s\n", 
                  gen, real, prediccion, error, 
                  substr(regla, 1, 50)))
    } else {
      cat(sprintf("  %s: Predicción = NA | Regla: %s\n", 
                  gen, substr(regla, 1, 50)))
    }
  }
  
  cat(sprintf("  Errores en esta transición (primeros 5 genes): %d\n", errors_this_t))
}

# Ahora calcular todo
cat("\n=== CALCULANDO PARA TODAS LAS TRANSICIONES ===\n")

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
      errors_por_gen[gen] <- errors_por_gen[gen] + error
      predictions_por_gen[gen] <- predictions_por_gen[gen] + 1
    } else {
      total_na <- total_na + 1
    }
  }
}

# Resultados
cat("\n========== RESULTADOS ==========\n")
cat(sprintf("Total predicciones: %d\n", total_predictions))
cat(sprintf("Total errores: %d\n", total_errors))
cat(sprintf("Total NA: %d\n", total_na))
cat(sprintf("Tasa de error: %.4f\n", total_errors / total_predictions))

# AQUÍ ESTÁ EL PROBLEMA - veamos ambas formas de calcular
dynacc_method1 <- 1 - (total_errors / total_predictions)
dynacc_method2 <- 1 - (total_errors / (n_genes * (n_timepoints - 1)))

cat(sprintf("\nDynamic Accuracy (Método 1 - predicciones válidas): %.6f\n", dynacc_method1))
cat(sprintf("Dynamic Accuracy (Método 2 - total teórico): %.6f\n", dynacc_method2))

cat(sprintf("\nDenominador método 1: %d (predicciones realmente hechas)\n", total_predictions))
cat(sprintf("Denominador método 2: %d (genes × transiciones)\n", n_genes * (n_timepoints - 1)))

# Ver genes con más errores
cat("\n=== TOP 10 GENES CON MÁS ERRORES ===\n")
dynacc_por_gen <- 1 - (errors_por_gen / predictions_por_gen)
worst_genes <- sort(dynacc_por_gen)[1:min(10, length(dynacc_por_gen))]
for(i in seq_along(worst_genes)) {
  gen <- names(worst_genes)[i]
  cat(sprintf("%s: DynAcc=%.4f (Errores=%d, Pred=%d)\n", 
              gen, worst_genes[i], 
              errors_por_gen[gen], 
              predictions_por_gen[gen]))
}