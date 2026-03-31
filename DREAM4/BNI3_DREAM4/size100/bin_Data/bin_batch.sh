#!/bin/bash

numeros=(5)
repeticiones=(5)

ruta_algoritmos="/home/lahumada/disco1/BNI3/Binarization"
ruta_raw="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/raw_Data"
ruta_out="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/bin_Data"
size=100

# Primero cambiar "Time" por "ID" en todos los archivos raw
echo "Cambiando 'Time' por 'ID' en archivos raw..."
for numero in "${numeros[@]}"
do
    for rep in "${repeticiones[@]}"
    do
        archivo="${ruta_raw}/size${size}_${numero}/size${size}_${numero}_${rep}"
        
        if [ -f "$archivo" ]; then
            sed -i '1s/^Time/ID/' "$archivo"
            echo "  Modificado: size${size}_${numero}_${rep}"
        else
            echo "  Archivo no encontrado: $archivo"
        fi
    done
done

echo ""
echo "Ejecutando algoritmos de binarización..."

# Luego ejecutar los algoritmos de binarización
for numero in "${numeros[@]}"
do
    for rep in "${repeticiones[@]}"
    do
        # Crear directorio de salida si no existe
        mkdir -p ${ruta_out}/size${size}_${numero}
        
        # Ejecutar SSD
        python ${ruta_algoritmos}/BNI3_SSD.py \
            -i ${ruta_raw}/size${size}_${numero}/size${size}_${numero}_${rep} \
            -o ${ruta_out}/size${size}_${numero}/bin_SSD_size${size}_${numero}_${rep}.tsv
        
        # Ejecutar WCSS
        python ${ruta_algoritmos}/BNI3_WCSS.py \
            -i ${ruta_raw}/size${size}_${numero}/size${size}_${numero}_${rep} \
            -o ${ruta_out}/size${size}_${numero}/bin_WCSS_size${size}_${numero}_${rep}.tsv
        
        echo "  Completado: dataset ${numero}, repetición ${rep}"
    done
done

echo ""
echo "Proceso finalizado"