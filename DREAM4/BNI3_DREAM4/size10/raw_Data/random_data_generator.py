#!/usr/bin/env python3
"""
Script para generar subconjuntos aleatorios de cursos temporales
de tablas de expresión génica manteniendo el orden temporal.
"""

import os
import random
import pandas as pd
from pathlib import Path
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Genera subconjuntos aleatorios de cursos temporales'
    )
    parser.add_argument(
        '--input-dir',
        default='/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/raw_Data',
        help='Directorio base con los datos'
    )
    parser.add_argument(
        '--n-samples',
        type=int,
        default=10,
        help='Número de muestras aleatorias a generar por archivo (default: 10)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Semilla para reproducibilidad (opcional)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Mostrar información de debug sobre archivos encontrados'
    )
    return parser.parse_args()

def get_timepoint_columns(df):
    """Obtiene las columnas que son timepoints (excluyendo 'ID')"""
    return [col for col in df.columns if col != 'ID']

def sample_timepoints(timepoints, n_points):
    """
    Selecciona n_points timepoints aleatorios manteniendo el orden temporal
    
    Args:
        timepoints: Lista de timepoints disponibles
        n_points: Número de timepoints a seleccionar
    
    Returns:
        Lista ordenada de timepoints seleccionados
    """
    # Convertir a índices numéricos para ordenar
    indices = list(range(len(timepoints)))
    selected_indices = sorted(random.sample(indices, n_points))
    return [timepoints[i] for i in selected_indices]

def is_data_file(filepath):
    """
    Verifica si un archivo es un archivo de datos válido
    (no es directorio, no termina en _5times, y es un archivo regular)
    """
    if not filepath.is_file():
        return False
    
    if '_5times' in filepath.name:
        return False
    
    # Verificar si parece ser un archivo de datos
    # (podemos intentar leerlo como TSV)
    try:
        # Intentar leer las primeras líneas
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            # Verificar que tenga formato de tabla con tabs
            if '\t' in first_line:
                return True
    except:
        pass
    
    return False

def process_file(input_file, output_dir, n_samples, base_name):
    """
    Procesa un archivo y genera múltiples versiones con diferentes
    números de timepoints
    
    Args:
        input_file: Ruta al archivo de entrada
        output_dir: Directorio de salida
        n_samples: Número de muestras aleatorias por combinación
        base_name: Nombre base para los archivos de salida
    """
    # Leer el archivo
    df = pd.read_csv(input_file, sep='\t')
    timepoints = get_timepoint_columns(df)
    
    print(f"  Procesando: {input_file.name}")
    print(f"    Timepoints disponibles: {len(timepoints)}")
    
    # Generar para 3, 4 y 5 timepoints
    for n_times in [3, 4, 5]:
        for sample_idx in range(n_samples):
            # Seleccionar timepoints aleatorios
            selected_times = sample_timepoints(timepoints, n_times)
            
            # Crear nuevo dataframe con timepoints seleccionados
            new_df = df[['ID'] + selected_times].copy()
            
            # Generar nombre de archivo de salida
            output_name = f"{base_name}_{n_times}times_random_{sample_idx+1}.tsv"
            output_path = output_dir / output_name
            
            # Guardar
            new_df.to_csv(output_path, sep='\t', index=False)
            
            if sample_idx == 0:  # Solo mostrar el primero para no saturar
                print(f"    Generado: {output_name} (timepoints: {', '.join(selected_times)})")

def find_data_files(directory, debug=False):
    """
    Busca recursivamente archivos de datos en un directorio
    excluyendo aquellos con '_5times' en el nombre
    """
    data_files = []
    
    if debug:
        print(f"\n  Buscando en: {directory}")
    
    for root, dirs, files in os.walk(directory):
        # Filtrar directorios con _5times
        dirs[:] = [d for d in dirs if '_5times' not in d]
        
        for file in files:
            full_path = Path(root) / file
            
            if is_data_file(full_path):
                data_files.append(full_path)
                if debug:
                    print(f"    Encontrado: {full_path}")
    
    return data_files

def main():
    args = parse_arguments()
    
    # Establecer semilla si se proporciona
    if args.seed is not None:
        random.seed(args.seed)
        print(f"Usando semilla: {args.seed}")
    
    input_dir = Path(args.input_dir)
    output_base = input_dir / 'randomSubset'
    
    # Crear directorio de salida
    output_base.mkdir(exist_ok=True)
    print(f"Directorio de salida: {output_base}\n")
    
    # Procesar cada size10_X
    total_processed = 0
    for size_dir_name in ['size10_1', 'size10_2', 'size10_3', 'size10_4', 'size10_5']:
        size_dir = input_dir / size_dir_name
        
        if not size_dir.exists():
            print(f"Advertencia: {size_dir} no existe, saltando...")
            continue
        
        print(f"\nProcesando {size_dir_name}:")
        
        # Crear subdirectorio en la salida
        output_subdir = output_base / size_dir_name
        output_subdir.mkdir(exist_ok=True)
        
        # Buscar todos los archivos de datos recursivamente
        data_files = find_data_files(size_dir, debug=args.debug)
        
        if not data_files:
            print(f"  ⚠ No se encontraron archivos de datos en {size_dir_name}")
            continue
        
        print(f"  Encontrados {len(data_files)} archivos de datos")
        
        for data_file in data_files:
            # Extraer nombre base del subdirectorio
            # Por ejemplo: size10_1/size10_1_1/archivo -> size10_1_1
            relative_path = data_file.relative_to(size_dir)
            if len(relative_path.parts) > 1:
                base_name = relative_path.parts[0]
            else:
                base_name = size_dir_name
            
            try:
                process_file(
                    data_file,
                    output_subdir,
                    args.n_samples,
                    base_name
                )
                total_processed += 1
            except Exception as e:
                print(f"  ✗ Error procesando {data_file}: {e}")
    
    print(f"\n✓ Proceso completado. Archivos guardados en: {output_base}")
    print(f"✓ Archivos originales procesados: {total_processed}")
    
    # Resumen
    total_files = sum(1 for _ in output_base.rglob('*.tsv'))
    print(f"✓ Total de archivos generados: {total_files}")
    print(f"✓ Archivos por dataset original: {args.n_samples * 3} (3, 4 y 5 timepoints)")

if __name__ == '__main__':
    main()