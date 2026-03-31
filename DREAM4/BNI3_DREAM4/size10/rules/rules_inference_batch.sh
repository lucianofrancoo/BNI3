#!/bin/bash

numeros=(1 2 3 4 5)
repeticiones=(1 2 3 4 5)

ruta_algoritmo="/home/lahumada/disco1/BNI3/Rules_Inference"
ruta_raw="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/raw_Data"
ruta_bin="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/bin_Data"
ruta_out="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/rules"
summary_file="${ruta_bin}/summary_patterns.tsv"
size=10

# Verificar que existe el archivo de resumen
if [ ! -f "$summary_file" ]; then
    echo "Error: No se encuentra el archivo $summary_file"
    echo "Ejecuta primero el script de análisis de patrones"
    exit 1
fi

echo "Iniciando inferencia de reglas booleanas (excluyendo archivos _5times)..."
echo ""

for numero in "${numeros[@]}"
do
    for rep in "${repeticiones[@]}"
    do
        # Buscar archivos raw disponibles para este dataset y repetición (excluyendo _5times)
        raw_pattern="${ruta_raw}/size${size}_${numero}/size${size}_${numero}_${rep}"
        
        # Encontrar todos los archivos que NO terminen en _5times
        for raw_file in ${raw_pattern}*; do
            # Verificar que el archivo existe y NO contiene _5times
            if [ ! -f "$raw_file" ] || [[ "$raw_file" == *"_5times"* ]]; then
                continue
            fi
            
            # Extraer el sufijo del archivo (todo lo que viene después de size10_X_Y)
            basename_file=$(basename "$raw_file")
            suffix="${basename_file#size${size}_${numero}_${rep}}"
            
            # Si no hay sufijo (archivo base sin extensión), usar string vacío
            if [ "$suffix" == "$basename_file" ]; then
                suffix=""
            fi
            
            echo "Procesando: Dataset ${numero}, Repetición ${rep}${suffix}"
            
            # Extraer unique_ratio para SSD y WCSS (excluyendo _5times)
            # Construir el nombre del archivo binarizado esperado
            if [ -z "$suffix" ]; then
                file_pattern="bin_(SSD|WCSS)_size${size}_${numero}_${rep}.tsv"
            else
                file_pattern="bin_(SSD|WCSS)_size${size}_${numero}_${rep}${suffix}.tsv"
            fi
            
            ssd_ratio=$(awk -F'\t' -v ds="$numero" -v r="$rep" -v suf="$suffix" \
                'NR>1 && $1==ds && $2==r && $3=="SSD" && $4 !~ /_5times\.tsv$/ && $4 ~ suf {print $9; exit}' "$summary_file")
            
            wcss_ratio=$(awk -F'\t' -v ds="$numero" -v r="$rep" -v suf="$suffix" \
                'NR>1 && $1==ds && $2==r && $3=="WCSS" && $4 !~ /_5times\.tsv$/ && $4 ~ suf {print $9; exit}' "$summary_file")
            
            # Verificar que se encontraron los valores
            if [ -z "$ssd_ratio" ] || [ -z "$wcss_ratio" ]; then
                echo "  ✗ No se encontraron datos para Dataset ${numero}, Rep ${rep}${suffix}"
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
            
            # Definir archivo binarizado
            bin_file="${ruta_bin}/size${size}_${numero}/bin_${method}_size${size}_${numero}_${rep}${suffix}.tsv"
            
            # Verificar que existe el archivo binarizado
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
                -o "${ruta_out}/size${size}_${numero}/rules_${method}_size${size}_${numero}_${rep}${suffix}" \
                -p 50 -n_rep 10 -pop 50 -gen 50
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Completado exitosamente"
            else
                echo "  ✗ Error en la ejecución"
            fi
            
            echo ""
        done
    done
done

echo "========================================="
echo "Proceso de inferencia finalizado"
echo "Archivos de reglas guardados en: ${ruta_out}"
echo "========================================="