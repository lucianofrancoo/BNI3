# ==== DIAGNÓSTICO ESPECÍFICO SIZE100 ====

library(dplyr)

ruta_reglas <- "/home/lahumada/disco1/BNI3/BoolNet_comparison/size10/size10_1/reglas_Boolnet_bestfit_size10_1_1.tsv"
ruta_matriz_bin <- "/home/lahumada/disco1/BNI3/BoolNet_comparison/size10/size10_1/bin_size10_1_1.tsv"

cat("========== ANÁLISIS DE REGLAS ==========\n")

# Leer reglas
reglas <- read.table(ruta_reglas,
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = FALSE)

cat(sprintf("Total genes con reglas: %d\n", nrow(reglas)))
cat(sprintf("Columnas en reglas: %s\n", paste(colnames(reglas), collapse=", ")))

# Ver distribución de errores
cat("\nDistribución de errores:\n")
print(table(reglas$Error))

cat("\nPrimeras 10 reglas:\n")
print(reglas[1:10, c("Gen", "Regla", "Error", "N_reguladores")])

# Verificar si hay genes sin reglas válidas o con errores altos
genes_con_error_alto <- reglas %>% filter(Error > 5)
cat(sprintf("\nGenes con Error > 5: %d\n", nrow(genes_con_error_alto)))

if(nrow(genes_con_error_alto) > 0) {
  cat("Ejemplos:\n")
  print(head(genes_con_error_alto[, c("Gen", "Error", "Regla")]))
}

cat("\n========== ANÁLISIS DE MATRIZ BINARIZADA ==========\n")

# Leer matriz
matriz_bin_raw <- read.table(ruta_matriz_bin,
                              sep = "\t",
                              header = TRUE,
                              stringsAsFactors = FALSE)

rownames(matriz_bin_raw) <- matriz_bin_raw$Gen
matriz_bin <- matriz_bin_raw[, -1]

cat(sprintf("Genes en matriz: %d\n", nrow(matriz_bin)))
cat(sprintf("Timepoints: %d\n", ncol(matriz_bin)))

# Verificar coincidencia
genes_reglas <- reglas$Gen
genes_matriz <- rownames(matriz_bin)

cat(sprintf("\nGenes en reglas: %d\n", length(genes_reglas)))
cat(sprintf("Genes en matriz: %d\n", length(genes_matriz)))
cat(sprintf("Genes en común: %d\n", length(intersect(genes_reglas, genes_matriz))))

genes_solo_reglas <- setdiff(genes_reglas, genes_matriz)
genes_solo_matriz <- setdiff(genes_matriz, genes_reglas)

if(length(genes_solo_reglas) > 0) {
  cat("\nGenes en reglas pero NO en matriz:\n")
  print(head(genes_solo_reglas, 10))
}

if(length(genes_solo_matriz) > 0) {
  cat("\nGenes en matriz pero NO en reglas:\n")
  print(head(genes_solo_matriz, 10))
}

cat("\n========== SIMULACIÓN DE PREDICCIONES ==========\n")

# Función de evaluación
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

# Simular primera transición en detalle
cat("\nPrimera transición (timepoint 1 -> 2):\n")

estado_t1 <- as.numeric(matriz_bin[, 1])
names(estado_t1) <- rownames(matriz_bin)
estado_t2_real <- as.numeric(matriz_bin[, 2])

errores_t1 <- 0
predicciones_t1 <- 0

for(i in 1:min(20, nrow(reglas))) {
  gen <- reglas$Gen[i]
  regla <- reglas$Regla[i]
  error_reportado <- reglas$Error[i]
  
  prediccion <- evaluate_boolean_rule(regla, estado_t1)
  
  gen_idx <- which(rownames(matriz_bin) == gen)
  real <- estado_t2_real[gen_idx]
  
  if(!is.na(prediccion)) {
    error_calc <- abs(real - prediccion)
    errores_t1 <- errores_t1 + error_calc
    predicciones_t1 <- predicciones_t1 + 1
    
    if(error_calc > 0 || error_reportado > 0) {
      cat(sprintf("%s: Real=%d, Pred=%d, Error=%d (Reportado=%d)\n", 
                  gen, real, prediccion, error_calc, error_reportado))
    }
  }
}

cat(sprintf("\nErrores en primera transición: %d/%d predicciones\n", 
            errores_t1, predicciones_t1))

# Ahora todas las transiciones
cat("\n========== CÁLCULO COMPLETO ==========\n")

total_errors <- 0
total_predictions <- 0

for(t in 1:(ncol(matriz_bin) - 1)) {
  estado_actual <- as.numeric(matriz_bin[, t])
  names(estado_actual) <- rownames(matriz_bin)
  estado_real <- as.numeric(matriz_bin[, t + 1])
  
  for(i in 1:nrow(reglas)) {
    gen <- reglas$Gen[i]
    regla <- reglas$Regla[i]
    
    prediccion <- evaluate_boolean_rule(regla, estado_actual)
    
    if(!is.na(prediccion)) {
      gen_idx <- which(rownames(matriz_bin) == gen)
      error <- abs(estado_real[gen_idx] - prediccion)
      total_errors <- total_errors + error
      total_predictions <- total_predictions + 1
    }
  }
}

dynacc_calculado <- 1 - (total_errors / total_predictions)

cat(sprintf("Total predicciones: %d\n", total_predictions))
cat(sprintf("Total errores: %d\n", total_errors))
cat(sprintf("Dynamic Accuracy calculado: %.6f\n", dynacc_calculado))

# Calcular Dynamic Accuracy esperado basado en Error promedio
error_promedio <- mean(reglas$Error)
n_timepoints <- ncol(matriz_bin)
errores_esperados <- error_promedio * nrow(reglas)
dynacc_esperado <- 1 - (errores_esperados / (nrow(reglas) * (n_timepoints - 1)))

cat(sprintf("\n=== COMPARACIÓN ===\n"))
cat(sprintf("Error promedio reportado: %.2f\n", error_promedio))
cat(sprintf("Errores esperados por transición: %.2f\n", errores_esperados))
cat(sprintf("Dynamic Accuracy esperado: %.6f\n", dynacc_esperado))
cat(sprintf("Dynamic Accuracy calculado: %.6f\n", dynacc_calculado))
cat(sprintf("Diferencia: %.6f\n", abs(dynacc_esperado - dynacc_calculado)))