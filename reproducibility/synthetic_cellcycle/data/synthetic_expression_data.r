library(BoolNet)

# Run from this folder (data/)

# Load the model
data(cellcycle)

# Inspect it
print(cellcycle)
cellcycle$genes  
cellcycle$interactions

# ============================================
# STEP 1: PARAMETERS
# ============================================
set.seed(123)  # So the same data are always generated

n_replicates <- 3      # Number of biological replicates
n_timepoints <- 10      # Number of time points
n_genes <- length(cellcycle$genes)
gene_names <- cellcycle$genes

# Number of genes
print(paste("Number of genes:", n_genes))

# ============================================
# STEP 2: GENERATE BOOLEAN TRAJECTORIES
# ============================================

all_boolean_trajectories <- list()

for (rep in 1:n_replicates) {
  
  cat("\n=== Replicate", rep, "===\n")
  
  # Random initial state
  initial_state <- sample(c(0, 1), n_genes, replace = TRUE)
  names(initial_state) <- gene_names
  
  cat("Initial state:\n")
  print(initial_state)
  
  # Simulate - returns a data.frame directly
  trajectory_df <- getPathToAttractor(cellcycle, initial_state, includeAttractorStates = "all")
  
  cat("Trajectory length:", nrow(trajectory_df), "steps\n")
  
  # Convert to a matrix to make it easier to work with
  boolean_states <- as.matrix(trajectory_df)
  
  # Adjust the length to n_timepoints
  current_length <- nrow(boolean_states)
  
  if (current_length < n_timepoints) {
    # Repeat the last state (the attractor)
    last_state <- boolean_states[current_length, ]
    for (i in (current_length + 1):n_timepoints) {
      boolean_states <- rbind(boolean_states, last_state)
    }
  } else if (current_length > n_timepoints) {
    # Take the first n_timepoints
    boolean_states <- boolean_states[1:n_timepoints, ]
  }
  
  all_boolean_trajectories[[rep]] <- boolean_states
  
  cat("Boolean states saved:\n")
  print(boolean_states)
}

# Parameters
low_mean <- 3
high_mean <- 10
noise_sd <- 2.5

cat("Low expression (OFF):", low_mean, "±", noise_sd, "\n")
cat("High expression (ON):", high_mean, "±", noise_sd, "\n")

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

cat("\nExample, replicate 1 (continuous):\n")
print(all_continuous_expression[[1]])

# Average the replicates
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
# TRANSPOSE: genes in rows, time points in columns
# ============================================

avg_expression_transposed <- t(avg_expression)

cat("\n=== Transposed matrix ===\n")
print(avg_expression_transposed)
print(paste("Dimensions:", nrow(avg_expression_transposed), "genes x", 
            ncol(avg_expression_transposed), "timepoints"))

# Convert to a data.frame with an ID column
df_expression <- as.data.frame(avg_expression_transposed)

# Add an ID column with the gene names
df_expression <- cbind(ID = rownames(df_expression), df_expression)

# Remove rownames
rownames(df_expression) <- NULL

cat("\nFinal data.frame:\n")
print(df_expression)

# Save
write.table(df_expression, 
            "synthetic_cellcycle_expression.tsv",
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE,
            col.names=TRUE)

cat("\n✓ File saved: synthetic_cellcycle_expression.tsv\n")
cat("Format: genes in rows, time points in columns, ID in the first column\n")
