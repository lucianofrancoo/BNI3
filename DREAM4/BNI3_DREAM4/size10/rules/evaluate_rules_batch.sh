#!/bin/bash

# Script de evaluación para datasets DREAM4 size10
# Procesa todas las carpetas EXCEPTO las que terminan en _5times

# ============================================
# CONFIGURACIÓN
# ============================================

# Datasets a procesar (size10_1 hasta size10_5)
datasets=(5)

# Rutas
ruta_algoritmo="/home/lahumada/disco1/BNI3/Rules_Inference"
ruta_bin="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/bin_Data"
ruta_rules="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/rules"
ruta_eval="/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/evaluations"

# Parámetros de evaluación
N_PROCESSES=10        # Ajusta según tus cores disponibles

# ============================================
# INICIO
# ============================================

# Crear directorio de evaluaciones
mkdir -p "${ruta_eval}"

# Log file
timestamp=$(date '+%Y%m%d_%H%M%S')
log_file="${ruta_eval}/evaluation_log_no5times_${timestamp}.txt"

echo "==========================================="
echo "EVALUACIÓN DREAM4 SIZE10 (NO _5times)"
echo "==========================================="
echo ""
echo "Configuración:"
echo "  - Datasets: size10_${datasets[@]}"
echo "  - Excluye carpetas: *_5times"
echo "  - Total a procesar: ~75-100 evaluaciones"
echo "  - Procesos paralelos: ${N_PROCESSES}"
echo ""
echo "Salida: ${ruta_eval}"
echo "Log: ${log_file}"
echo ""

# Inicializar log
echo "Evaluación DREAM4 size10 (NO _5times) - $(date)" > "$log_file"
echo "==========================================" >> "$log_file"

# Contadores
total=0
exitosos=0
fallidos=0
excluidos=0

# ============================================
# PROCESAMIENTO
# ============================================

for dataset_num in "${datasets[@]}"
do
    dataset_dir="${ruta_rules}/size10_${dataset_num}"
    
    # Verificar que existe el directorio del dataset
    if [ ! -d "$dataset_dir" ]; then
        echo "⚠️  Directorio no encontrado: size10_${dataset_num}"
        echo "   Saltando..."
        echo ""
        continue
    fi
    
    echo "-------------------------------------------"
    echo "DATASET: size10_${dataset_num}"
    echo "-------------------------------------------"
    
    # Buscar todas las carpetas que empiezan con rules_
    for rules_dir in ${dataset_dir}/rules_*; do
        
        # Verificar que es un directorio
        if [ ! -d "$rules_dir" ]; then
            continue
        fi
        
        dirname=$(basename "$rules_dir")
        
        # EXCLUIR carpetas que terminan en _5times
        if [[ "$dirname" == *"_5times" ]]; then
            ((excluidos++))
            echo "  ⏭️  Excluido: $(basename $rules_dir) (_5times)"
            continue
        fi
        
        # Verificar que existe rules_by_gene.tsv
        rules_file="${rules_dir}/rules_by_gene.tsv"
        if [ ! -f "$rules_file" ]; then
            echo "  ⚠️  No se encuentra rules_by_gene.tsv en $(basename $rules_dir)"
            continue
        fi
        
        # Extraer información del nombre de la carpeta
        # Formato: rules_{METHOD}_size10_{NUM}_{REP}
        # Usar regex para extraer METHOD, NUM, REP (sin _5times)
        if [[ "$dirname" =~ rules_(SSD|WCSS)_size10_([0-9]+)_([0-9]+)$ ]]; then
            method="${BASH_REMATCH[1]}"
            num="${BASH_REMATCH[2]}"
            rep="${BASH_REMATCH[3]}"
        else
            echo "  ⚠️  Formato inesperado: $dirname"
            continue
        fi
        
        ((total++))
        
        echo ""
        echo "  📁 $(basename $rules_dir)"
        echo "     Método: ${method}, Dataset: ${num}, Repetición: ${rep}"
        
        # Construir path del archivo binarizado correspondiente (sin _5times)
        bin_file="${ruta_bin}/size10_${num}/bin_${method}_size10_${num}_${rep}.tsv"
        
        # Verificar que existe el archivo binarizado
        if [ ! -f "$bin_file" ]; then
            echo "     ❌ Archivo binarizado no encontrado:"
            echo "        $(basename $bin_file)"
            echo "     Saltando..."
            ((fallidos++))
            echo "size10_${num}, rep ${rep}, ${method}: FALLIDO (sin binarizado)" >> "$log_file"
            continue
        fi
        
        echo "     ✓ Archivo binarizado: $(basename $bin_file)"
        
        # Guardar directamente en la carpeta de reglas (NO crear carpeta nueva)
        eval_output_dir="${rules_dir}"
        
        echo "     📂 Salida: Misma carpeta de reglas"
        echo ""
        echo "     🚀 Ejecutando evaluación con defaults inteligentes..."
        
        # Ejecutar evaluación (los archivos se crearán en rules_dir)
        python3 ${ruta_algoritmo}/BNI3_Evaluate_rules.py \
            -i "$rules_file" \
            -m "$bin_file" \
            -o "${eval_output_dir}/evaluation_results.tsv" \
            -n ${N_PROCESSES} \
            -v
        
        exit_code=$?
        
        # Procesar resultado
        if [ $exit_code -eq 0 ]; then
            echo ""
            echo "     ✅ Evaluación completada"
            
            # Verificar archivos generados (en la misma carpeta)
            if [ -f "${eval_output_dir}/evaluation_results.tsv" ] && \
               [ -f "${eval_output_dir}/rules_by_gene_evaluated.tsv" ]; then
                
                # Extraer mejor score
                best_score=$(tail -n +2 "${eval_output_dir}/evaluation_results.tsv" | head -n 1 | cut -f2)
                n_combos=$(tail -n +2 "${eval_output_dir}/evaluation_results.tsv" | wc -l)
                
                echo "     📊 Mejor final_score: ${best_score}"
                echo "     📊 Combinaciones evaluadas: ${n_combos}"
                
                ((exitosos++))
                echo "size10_${num}, rep ${rep}, ${method}: EXITOSO (score=${best_score}, combos=${n_combos})" >> "$log_file"
            else
                echo "     ⚠️  Archivos incompletos"
                ((fallidos++))
                echo "size10_${num}, rep ${rep}, ${method}: PARCIAL" >> "$log_file"
            fi
        else
            echo ""
            echo "     ❌ Error (código: ${exit_code})"
            ((fallidos++))
            echo "size10_${num}, rep ${rep}, ${method}: FALLIDO (error=${exit_code})" >> "$log_file"
        fi
        
        echo ""
    done
done

# ============================================
# RESUMEN FINAL
# ============================================

echo "==========================================="
echo "RESUMEN DE EVALUACIÓN"
echo "==========================================="
echo ""
echo "📊 Estadísticas:"
echo "  - Total encontrados: ${total}"
echo "  - Excluidos (_5times): ${excluidos}"
echo "  - Procesados: ${total}"
echo "  - Exitosos: ${exitosos}"
echo "  - Fallidos: ${fallidos}"

if [ $total -gt 0 ]; then
    success_rate=$(echo "scale=1; ($exitosos * 100) / $total" | bc)
    echo "  - Tasa de éxito: ${success_rate}%"
fi

echo ""
echo "📁 Resultados en carpetas de reglas"
echo "📄 Log: ${log_file}"
echo ""

# Generar resumen CSV
summary_csv="${ruta_eval}/evaluation_summary_no5times_${timestamp}.csv"
echo "Dataset,Repeticion,Metodo,FinalScore,Combinaciones,CarpetaReglas" > "$summary_csv"

# Buscar en carpetas de reglas (no en evaluations) - EXCLUYENDO _5times
for dataset_num in "${datasets[@]}"
do
    dataset_dir="${ruta_rules}/size10_${dataset_num}"
    
    if [ ! -d "$dataset_dir" ]; then
        continue
    fi
    
    for rules_dir in ${dataset_dir}/rules_*; do
        # Saltar carpetas _5times
        if [[ "$(basename $rules_dir)" == *"_5times" ]]; then
            continue
        fi
        
        if [ -f "${rules_dir}/evaluation_results.tsv" ]; then
            dirname=$(basename "$rules_dir")
            
            # Extraer info del nombre: rules_{METHOD}_size10_{NUM}_{REP}
            if [[ "$dirname" =~ rules_(SSD|WCSS)_size10_([0-9]+)_([0-9]+)$ ]]; then
                method="${BASH_REMATCH[1]}"
                ds="${BASH_REMATCH[2]}"
                rep="${BASH_REMATCH[3]}"
                
                score=$(tail -n +2 "${rules_dir}/evaluation_results.tsv" | head -n 1 | cut -f2)
                combos=$(tail -n +2 "${rules_dir}/evaluation_results.tsv" | wc -l)
                
                echo "${ds},${rep},${method},${score},${combos},${dirname}" >> "$summary_csv"
            fi
        fi
    done
done

if [ -f "$summary_csv" ]; then
    echo "📊 Resumen CSV: ${summary_csv}"
    echo ""
    
    # Mostrar preview del CSV
    echo "Preview del resumen:"
    head -n 11 "$summary_csv" | column -t -s','
    echo "..."
    
    # Mostrar estadísticas del CSV
    n_rows=$(tail -n +2 "$summary_csv" | wc -l)
    echo ""
    echo "Total de evaluaciones en CSV: ${n_rows}"
fi

echo ""
echo "✅ Proceso finalizado: $(date)"
echo "==========================================="

# Log final
echo "" >> "$log_file"
echo "==========================================" >> "$log_file"
echo "RESUMEN FINAL" >> "$log_file"
echo "Total encontrados: ${total}" >> "$log_file"
echo "Excluidos (_5times): ${excluidos}" >> "$log_file"
echo "Exitosos: ${exitosos}" >> "$log_file"
echo "Fallidos: ${fallidos}" >> "$log_file"
echo "Finalizado: $(date)" >> "$log_file"