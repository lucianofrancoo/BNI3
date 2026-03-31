library(BoolNet)

setwd("/home/lahumada/Desktop/BNI3/synthetic_cellcycle/binary")

# Cargar el modelo
data(cellcycle)

# Ver qué contiene
print(cellcycle)
cellcycle$genes  
cellcycle$interactions

# ============================================
# PASO 1: CONFIGURAR PARÁMETROS
# ============================================
set.seed(123)  # Para que siempre generes los mismos datos

n_replicates <- 3      # Número de réplicas biológicas
n_timepoints <- 10      # Número de timepoints
n_genes <- length(cellcycle$genes)
gene_names <- cellcycle$genes

# Ver cuántos genes tenemos
print(paste("Número de genes:", n_genes))

# ============================================
# PASO 2: GENERAR TRAYECTORIAS BOOLEANAS
# ============================================

all_boolean_trajectories <- list()

for (rep in 1:n_replicates) {
  
  cat("\n=== Réplica", rep, "===\n")
  
  # Estado inicial aleatorio
  initial_state <- sample(c(0, 1), n_genes, replace = TRUE)
  names(initial_state) <- gene_names
  
  cat("Estado inicial:\n")
  print(initial_state)
  
  # Simular - devuelve un data.frame directamente
  trajectory_df <- getPathToAttractor(cellcycle, initial_state, includeAttractorStates = "all")
  
  cat("Longitud de trayectoria:", nrow(trajectory_df), "pasos\n")
  
  # Convertir a matriz para trabajar más fácil
  boolean_states <- as.matrix(trajectory_df)
  
  # Ajustar longitud a n_timepoints
  current_length <- nrow(boolean_states)
  
  if (current_length < n_timepoints) {
    # Repetir último estado (el atractor)
    last_state <- boolean_states[current_length, ]
    for (i in (current_length + 1):n_timepoints) {
      boolean_states <- rbind(boolean_states, last_state)
    }
  } else if (current_length > n_timepoints) {
    # Tomar primeros n_timepoints
    boolean_states <- boolean_states[1:n_timepoints, ]
  }
  
  all_boolean_trajectories[[rep]] <- boolean_states
  
  cat("Estados Boolean guardados:\n")
  print(boolean_states)
}

# Parámetros
low_mean <- 3
high_mean <- 10
noise_sd <- 2.5

cat("Expresión baja (OFF):", low_mean, "±", noise_sd, "\n")
cat("Expresión alta (ON):", high_mean, "±", noise_sd, "\n")

all_continuous_expression <- list()

for (rep in 1:n_replicates) {
  
  boolean_matrix <- all_boolean_trajectories[[rep]]
  continuous_matrix <- boolean_matrix
  
  for (t in seq_len(nrow(boolean_matrix))) {
    for (g in seq_len(ncol(boolean_matrix))) {
      
      if (boolean_matrix[t, g] == 0) {
        value <- rnorm(1, mean = low_mean, sd = noise_sd)
      } else {
        value <- rnorm(1, mean = high_mean, sd = noise_sd)
      }
      
      continuous_matrix[t, g] <- max(0, value)
    }
  }
  
  all_continuous_expression[[rep]] <- continuous_matrix
}

cat("\nEjemplo réplica 1 (continuo):\n")
print(all_continuous_expression[[1]])

# Promediar réplicas
avg_expression <- matrix(0, nrow = n_timepoints, ncol = n_genes)
colnames(avg_expression) <- gene_names
rownames(avg_expression) <- paste0("T", 0:(n_timepoints - 1))

for (t in 1:n_timepoints) {
  for (g in 1:n_genes) {
    values <- numeric(n_replicates)
    for (rep in 1:n_replicates) {
      values[rep] <- all_continuous_expression[[rep]][t, g]
    }
    avg_expression[t, g] <- mean(values)
  }
}

# ============================================
# TRANSPONER: Genes en filas, Tiempos en columnas
# ============================================

avg_expression_transposed <- t(avg_expression)

cat("\n=== Matriz transpuesta ===\n")
print(avg_expression_transposed)
print(paste("Dimensiones:", nrow(avg_expression_transposed), "genes x", 
            ncol(avg_expression_transposed), "timepoints"))

# Convertir a data.frame con columna ID
df_expression <- as.data.frame(avg_expression_transposed)

# Añadir columna ID con los nombres de genes
df_expression <- cbind(ID = rownames(df_expression), df_expression)

# Eliminar rownames
rownames(df_expression) <- NULL

cat("\nData.frame final:\n")
print(df_expression)

# Guardar
write.table(df_expression, 
            "synthetic_cellcycle_expression.tsv",
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE,
            col.names=TRUE)

cat("\n✓ Archivo guardado: synthetic_cellcycle_expression.tsv\n")
cat("Formato: genes en filas, tiempos en columnas, ID en primera columna\n")
