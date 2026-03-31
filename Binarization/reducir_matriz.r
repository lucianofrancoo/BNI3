library("dplyr")

setwd("/home/lahumada/Desktop/redes_neuronales/BNI/Sequia_Avance/Ath/binarizacion")
# Cargar matriz de conteo
counts <- read.table("/home/lahumada/Desktop/redes_neuronales/BNI/Sequia_Avance/Ath/binarizacion/matrices_conteo/Conteo_genes_Aerts.tsv",
                     sep = "\t",
                     header = TRUE,
                     row.names = "ID")
head(counts)
dim(counts)

coldata <- read.table("/home/lahumada/Desktop/redes_neuronales/BNI/Sequia_Avance/Ath/binarizacion/matrices_conteo/coldata_Aerts",
                      sep = "\t",
                      header = TRUE,
                      row.names = "Name")
coldata

# Filtrar coldata para obtener solo las filas con Organ igual a "SeedlingShoot"
coldata <- coldata %>% filter(Tissue == "Rosette_leaf")
coldata <- coldata %>% filter(Treatment == "50_uM_ABA" | Treatment == "None")
coldata

# Reducir la matriz de conteo para tener solo el tratamiento
columns <- rownames(coldata)

# Filtrar countdata para obtener solo las columnas correspondientes a SeedlingShoot
counts <- counts[, columns]
head(counts)
dim(counts)

# Leer el archivo que contiene los genes de interés
genes_interes1 <- read.table("/home/lahumada/Desktop/redes_neuronales/BNI/Sequia_Avance/Ath/binarizacion/DGEs/lista_final",
                             header = TRUE)
head(genes_interes1)
dim(genes_interes1)
# genes_interes2 <- read.table("/home/lahumada/Desktop/redes_neuronales/BNI/ensayo_N_ABA/binarizacion/lista_genes_ABA",
#                              header = TRUE)
# head(genes_interes2)
# dim(genes_interes2)

# # Concatenar ambas tablas
# genes_total <- rbind(genes_interes1, genes_interes2)
# head(genes_total)
# dim(genes_total)

genes_total <- unique(genes_interes1)
dim(genes_total)

# Filtrar la matriz de conteo para obtener solo las filas correspondientes a los genes de interés
counts <- counts %>% filter(rownames(counts) %in% genes_total$ID)
head(counts)
dim(counts)

# Renombrar las filas de la matriz de conteo con los nombres de los genes de interés
#rownames(counts) <- genes_total$Name
rownames(counts) <- genes_total$Name[match(rownames(counts), genes_total$ID)]
head(t(counts))
dim(counts)

# Promediar los datos de conteo por condicion
# Transponer countdata para que las filas sean las muestras
counts_t <- t(counts)

# Convertir countdata_t en un data frame
counts_t <- as.data.frame(counts_t)

# Agregar la columna de Condition a countdata_df
counts_t$Condition <- coldata$Condition
head(counts_t)

# Calcular el promedio por Condition
counts_t_average <- counts_t %>%
  group_by(Condition) %>%
  summarise(across(everything(), mean))
head(counts_t_average)

# Transponer de nuevo para que las filas sean los genes
counts_average <- t(counts_t_average)

# Convertir de nuevo en data frame y ajustar nombres de columnas
counts_average <- as.data.frame(counts_average)
# Asignar la primera fila como nombres de las columnas
colnames(counts_average) <- counts_average[1, ]

# Eliminar la primera fila
counts_average <- counts_average[-1, ]
head(t(counts_average))

# Ordenar las columnas
column_numbers <- as.numeric(sub(".*_([0-9]+)$", "\\1", names(counts_average)))

# Ordenar las columnas basándose en esos números
counts_average <- counts_average[, order(column_numbers)]

# Guardar los nombres originales de las filas
rownames_originales <- rownames(counts_average)

# Convertir los valores de counts_average a numéricos
counts_average <- as.data.frame(lapply(counts_average, function(x) as.numeric(as.character(x))))

# Restaurar los nombres originales de las filas
rownames(counts_average) <- rownames_originales

# Mostrar las primeras filas de la matriz transpuesta
head(t(counts_average))
dim(counts_average)

counts_average$ID <- rownames(counts_average)
counts_average <- counts_average %>% select(ID, everything())

write.table(counts_average,
            paste0("matrices_conteo/counts_Aerts_genes95_ABA_Rosette_leaf.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

