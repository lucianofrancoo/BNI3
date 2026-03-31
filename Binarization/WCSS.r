# Librerías necesarias
library(dplyr)

# Establecer el directorio de trabajo
setwd("/home/lahumada/Desktop/BNI3/DREAM4/BNI3_DREAM4/Size10/bin_Data")

ruta_df <- "/home/lahumada/Desktop/BNI3/DREAM4/BNI3_DREAM4/Size10/raw_Data/size10_1_5"
# Leer el archivo de datos en un dataframe
df <- read.table(ruta_df,
                 header = TRUE,
                 sep = "\t",
                 row.names = "ID")

# Transponer el dataframe
t(df)

# Función para discretizar y reemplazar en el dataframe original
discretizar_y_reemplazar <- function(exp_matrix) {
  
  # Crear una copia del dataframe original
  df_discretizado <- exp_matrix
  
  # Iterar sobre cada fila del dataframe
  for (i in 1:nrow(exp_matrix)) { # nolint
    # Convertir la fila a un vector numérico
    fila <- as.numeric(exp_matrix[i, ])
    # Ordenar los valores de la fila
    fila_ordenada <- sort(fila)
    # Obtener el número de columnas
    num_cols <- length(fila)
    # Inicializar el mejor WCSS (Within-Cluster Sum of Squares) como infinito
    # Inicializamos con infinito para asegurarnos de encontrar un mejor valor
    mejor_wcss <- Inf
    # Inicializar el mejor clúster como NULL
    # Aún no tenemos una partición definida
    mejor_cluster <- NULL
    
    # Encontrar la mejor división en dos clústeres
    for (k in 1:(num_cols - 1)) {
      # Dividir la fila ordenada en dos clústeres
      cluster_1 <- fila_ordenada[1:k]
      cluster_2 <- fila_ordenada[(k + 1):num_cols]
      
      # Calcular la media de cada clúster
      mean_1 <- mean(cluster_1)
      mean_2 <- mean(cluster_2)
      
      # Calcular el WCSS para cada clúster
      wcss_1 <- sum((cluster_1 - mean_1)^2)
      wcss_2 <- sum((cluster_2 - mean_2)^2)
      
      # Calcular el WCSS total
      total_wcss <- wcss_1 + wcss_2
      
      # Actualizar el mejor WCSS y el mejor clúster si se encuentra una mejor división
      if (total_wcss < mejor_wcss) {
        mejor_wcss <- total_wcss
        mejor_cluster <- list(cluster_1 = cluster_1,
                              cluster_2 = cluster_2)
      }
    }
    
    # Reemplazar los valores en la matriz original con 0 y 1 según el mejor clúster encontrado
    df_discretizado[i, ] <- ifelse(exp_matrix[i, ] %in% mejor_cluster$cluster_1, 0, 1)
  }
  
  # Devolver el dataframe discretizado
  return(df_discretizado) # nolint
}

# Ejecutar la función y obtener la versión discretizada del dataframe
df_discretizado <- discretizar_y_reemplazar(df)

# Mostrar el resultado transpuesto
t(df_discretizado)

nombre_archivo <- basename(ruta_df)
write.table(t(df_discretizado),
            file = paste0("Bin_WCSS_", nombre_archivo),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
