#!/bin/bash

numeros=(1 2 3 4 5)

ruta_algoritmo="/home/lahumada/disco1/BNI3/Rules_Inference"
ruta_raw="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/raw_Data/randomSubset"
ruta_bin="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/bin_Data_randomSubset"
ruta_out="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/rules_randomSubset"
summary_file="${ruta_bin}/summary_patterns.tsv"
size=10

# Verificar que existe el archivo de resumen
if [ ! -f "$summary_file" ]; then
    echo "Error: No se encuentra el archivo $summary_file"
    echo "Ejecuta primero el script de análisis de patrones"
    exit 1
fi

echo "Iniciando inferencia de reglas booleanas para random subset..."
echo "========================================="
echo ""

total_procesados=0
total_errores=0

for numero in "${numeros[@]}"
do
    directorio_raw="${ruta_raw}/size${size}_${numero}"
    directorio_bin="${ruta_bin}/size${size}_${numero}"
    directorio_out="${ruta_out}/size${size}_${numero}"
    
    # Verificar que existen los directorios
    if [ ! -d "$directorio_raw" ]; then
        echo "ADVERTENCIA: No existe $directorio_raw, saltando..."
        continue
    fi
    
    if [ ! -d "$directorio_bin" ]; then
        echo "ADVERTENCIA: No existe $directorio_bin, saltando..."
        continue
    fi
    
    # Crear directorio de salida
    mkdir -p "$directorio_out"
    
    echo "Procesando size${size}_${numero}..."
    echo ""
    
    # Procesar cada archivo .tsv en raw
    for raw_file in ${directorio_raw}/*.tsv; do
        # Verificar que el archivo existe
        if [ ! -f "$raw_file" ]; then
            continue
        fi
        
        # Extraer el nombre base del archivo
        basename_file=$(basename "$raw_file" .tsv)
        
        echo "  Procesando: $basename_file"
        
        # Buscar unique_ratio para SSD y WCSS en el summary
        ssd_ratio=$(awk -F'\t' -v file="bin_SSD_${basename_file}.tsv" \
            'NR>1 && $4==file {print $9; exit}' "$summary_file")
        
        wcss_ratio=$(awk -F'\t' -v file="bin_WCSS_${basename_file}.tsv" \
            'NR>1 && $4==file {print $9; exit}' "$summary_file")
        
        # Verificar que se encontraron los valores
        if [ -z "$ssd_ratio" ] || [ -z "$wcss_ratio" ]; then
            echo "    ✗ No se encontraron datos en summary para $basename_file"
            echo "    (Buscando: bin_SSD_${basename_file}.tsv y bin_WCSS_${basename_file}.tsv)"
            total_errores=$((total_errores + 1))
            echo ""
            continue
        fi
        
        # Comparar y seleccionar método (en caso de empate, WCSS gana)
        if (( $(echo "$wcss_ratio >= $ssd_ratio" | bc -l) )); then
            method="WCSS"
            selected_ratio=$wcss_ratio
        else
            method="SSD"
            selected_ratio=$ssd_ratio
        fi
        
        echo "    SSD ratio: ${ssd_ratio}, WCSS ratio: ${wcss_ratio}"
        echo "    ✓ Seleccionado: ${method} (ratio: ${selected_ratio})"
        
        # Definir archivo binarizado
        bin_file="${directorio_bin}/bin_${method}_${basename_file}.tsv"
        
        # Verificar que existe el archivo binarizado
        if [ ! -f "$bin_file" ]; then
            echo "    ✗ Archivo binarizado no encontrado: $bin_file"
            total_errores=$((total_errores + 1))
            echo ""
            continue
        fi
        
        # Ejecutar inferencia de reglas
        echo "    → Ejecutando inferencia..."
        python3 ${ruta_algoritmo}/BNI3_Boolean_Rules_Inference.py \
            -i "$raw_file" \
            -i_binary "$bin_file" \
            -o "${directorio_out}/rules_${method}_${basename_file}" \
            -p 50 -n_rep 10 -pop 50 -gen 50
        
        if [ $? -eq 0 ]; then
            echo "    ✓ Completado exitosamente"
            total_procesados=$((total_procesados + 1))
        else
            echo "    ✗ Error en la ejecución"
            total_errores=$((total_errores + 1))
        fi
        
        echo ""
    done
done

echo "========================================="
echo "Proceso de inferencia finalizado"
echo ""
echo "Resumen:"
echo "  Archivos procesados exitosamente: $total_procesados"
echo "  Errores: $total_errores"
echo ""
echo "Archivos de reglas guardados en: ${ruta_out}"
echo "========================================="