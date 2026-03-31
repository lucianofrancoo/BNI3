#!/bin/bash

numeros=(1 2 3 4 5)

ruta_algoritmos="/home/lahumada/disco1/BNI3/Binarization"
ruta_raw="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/raw_Data/randomSubset"
ruta_out="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/bin_Data_randomSubset"
size=10

echo "Ejecutando algoritmos de binarización para archivos random subset..."
echo "========================================="

# Procesar cada size10_X
for numero in "${numeros[@]}"
do
    directorio_entrada="${ruta_raw}/size${size}_${numero}"
    directorio_salida="${ruta_out}/size${size}_${numero}"
    
    # Verificar que el directorio de entrada existe
    if [ ! -d "$directorio_entrada" ]; then
        echo "  ADVERTENCIA: No existe $directorio_entrada, saltando..."
        continue
    fi
    
    # Crear directorio de salida si no existe
    mkdir -p ${directorio_salida}
    
    echo ""
    echo "Procesando size${size}_${numero}..."
    
    # Contar archivos a procesar
    total_archivos=$(ls ${directorio_entrada}/*.tsv 2>/dev/null | wc -l)
    
    if [ $total_archivos -eq 0 ]; then
        echo "  ADVERTENCIA: No hay archivos .tsv en $directorio_entrada"
        continue
    fi
    
    echo "  Total de archivos a procesar: $total_archivos"
    
    contador=0
    
    # Procesar cada archivo .tsv en el directorio
    for archivo_entrada in ${directorio_entrada}/*.tsv
    do
        # Extraer el nombre base del archivo sin extensión
        nombre_base=$(basename "$archivo_entrada" .tsv)
        
        contador=$((contador + 1))
        echo ""
        echo "  [$contador/$total_archivos] Procesando: $nombre_base"
        
        # Ejecutar SSD
        echo "    → Ejecutando SSD..."
        python ${ruta_algoritmos}/BNI3_SSD.py \
            -i ${archivo_entrada} \
            -o ${directorio_salida}/bin_SSD_${nombre_base}.tsv
        
        if [ $? -ne 0 ]; then
            echo "    ✗ Error en SSD para $nombre_base"
        else
            echo "    ✓ SSD completado"
        fi
        
        # Ejecutar WCSS
        echo "    → Ejecutando WCSS..."
        python ${ruta_algoritmos}/BNI3_WCSS.py \
            -i ${archivo_entrada} \
            -o ${directorio_salida}/bin_WCSS_${nombre_base}.tsv
        
        if [ $? -ne 0 ]; then
            echo "    ✗ Error en WCSS para $nombre_base"
        else
            echo "    ✓ WCSS completado"
        fi
    done
    
    echo ""
    echo "  ✓ Completado size${size}_${numero}: $contador archivos procesados"
done

echo ""
echo "========================================="
echo "Proceso finalizado"
echo "Archivos procesados guardados en: ${ruta_out}"
echo ""
echo "Resumen:"
for numero in "${numeros[@]}"
do
    if [ -d "${ruta_out}/size${size}_${numero}" ]; then
        total=$(ls ${ruta_out}/size${size}_${numero}/*.tsv 2>/dev/null | wc -l)
        echo "  size${size}_${numero}: $total archivos binarizados"
    fi
done
echo "========================================="