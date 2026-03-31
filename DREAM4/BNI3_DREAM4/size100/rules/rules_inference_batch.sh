#!/bin/bash

numeros=(2 3 4 5)
repeticiones=(1 2 3 4 5)

ruta_algoritmo="/home/lahumada/disco1/BNI3/Rules_Inference"
ruta_raw="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/raw_Data"
ruta_bin="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/bin_Data"
ruta_out="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/rules"
summary_file="${ruta_bin}/summary_patterns.tsv"
size=100

# Verificar que existe el archivo de resumen
if [ ! -f "$summary_file" ]; then
    echo "Error: No se encuentra el archivo $summary_file"
    echo "Ejecuta primero el script de análisis de patrones"
    exit 1
fi

echo "Iniciando inferencia de reglas booleanas..."
echo ""

for numero in "${numeros[@]}"
do
    for rep in "${repeticiones[@]}"
    do
        echo "Procesando: Dataset ${numero}, Repetición ${rep}"
        
        # Extraer unique_ratio para SSD y WCSS usando awk
        ssd_ratio=$(awk -F'\t' -v ds="$numero" -v r="$rep" \
            'NR>1 && $1==ds && $2==r && $3=="SSD" {print $9}' "$summary_file")
        
        wcss_ratio=$(awk -F'\t' -v ds="$numero" -v r="$rep" \
            'NR>1 && $1==ds && $2==r && $3=="WCSS" {print $9}' "$summary_file")
        
        # Verificar que se encontraron los valores
        if [ -z "$ssd_ratio" ] || [ -z "$wcss_ratio" ]; then
            echo "  ✗ No se encontraron datos para Dataset ${numero}, Rep ${rep}"
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
        
        echo "  SSD ratio: ${ssd_ratio}, WCSS ratio: ${wcss_ratio}"
        echo "  ✓ Seleccionado: ${method} (ratio: ${selected_ratio})"
        
        # Definir archivos de entrada
        raw_file="${ruta_raw}/size${size}_${numero}/size${size}_${numero}_${rep}"
        bin_file="${ruta_bin}/size${size}_${numero}/bin_${method}_size${size}_${numero}_${rep}.tsv"
        
        # Verificar que existen los archivos
        if [ ! -f "$raw_file" ]; then
            echo "  ✗ Archivo raw no encontrado: $raw_file"
            echo ""
            continue
        fi
        
        if [ ! -f "$bin_file" ]; then
            echo "  ✗ Archivo binarizado no encontrado: $bin_file"
            echo ""
            continue
        fi
        
        # Crear directorio de salida si no existe
        mkdir -p "${ruta_out}/size${size}_${numero}"
        
        # Ejecutar inferencia de reglas
        echo "  Ejecutando inferencia..."
        python3 ${ruta_algoritmo}/BNI3_Boolean_Rules_Inference.py \
            -i "$raw_file" \
            -i_binary "$bin_file" \
            -o "${ruta_out}/size${size}_${numero}/rules_${method}_size${size}_${numero}_${rep}" \
            -p 20 -n_rep 1 -pop 50 -gen 50
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Completado exitosamente"
        else
            echo "  ✗ Error en la ejecución"
        fi
        
        echo ""
    done
done

echo "Proceso de inferencia finalizado"