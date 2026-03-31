import os
import pandas as pd
from pathlib import Path

# === Configuración ===
base_dir = "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size100/bin_Data"
output_summary = os.path.join(base_dir, "summary_patterns.tsv")

# === Función para procesar un archivo ===
def process_binarized_file(input_file):
    """
    Procesa un archivo binarizado y retorna los patrones y sus frecuencias
    """
    # Cargar matriz
    df = pd.read_csv(input_file, sep="\t")
    
    # Convertir columnas en secuencias temporales
    patterns = {}
    for gene in df.columns:
        if gene != "ID":  # Saltar columna ID si existe
            sequence = ''.join(df[gene].astype(str).tolist())
            patterns[gene] = sequence
    
    # Contar frecuencia de cada comportamiento
    pattern_counts = pd.Series(patterns).value_counts()
    
    # Construir tabla de salida
    output = pd.DataFrame([
        {"Gene": gene, 
         "Pattern": pat, 
         "Frequency": pattern_counts[pat]} 
        for gene, pat in patterns.items()
    ])
    
    # Ordenar: únicos arriba, compartidos abajo
    output = output.sort_values(by=["Frequency", "Pattern"], ascending=[True, True])
    
    # Guardar archivo individual
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_pattern.tsv"
    output.to_csv(output_file, sep="\t", index=False)
    
    return output, os.path.basename(input_file)

# === Buscar todos los archivos binarizados ===
all_files = []
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.startswith("bin_") and file.endswith(".tsv") and "_pattern" not in file:
            all_files.append(os.path.join(root, file))

print(f"Encontrados {len(all_files)} archivos binarizados")
print()

# === Procesar todos los archivos y crear resumen ===
summary_data = []

for file_path in sorted(all_files):
    print(f"Procesando: {os.path.basename(file_path)}")
    
    try:
        patterns_df, filename = process_binarized_file(file_path)
        
        # Extraer información del nombre del archivo
        parts = filename.replace("bin_", "").replace(".tsv", "").split("_")
        method = parts[0]  # SSD o WCSS
        size = parts[1]    # size10
        dataset = parts[2]  # número del dataset
        rep = parts[3]      # número de repetición
        
        # Calcular estadísticas
        n_genes = len(patterns_df)
        n_unique_patterns = len(patterns_df[patterns_df["Frequency"] == 1])
        n_shared_patterns = len(patterns_df["Pattern"].unique())
        max_frequency = patterns_df["Frequency"].max()
        
        summary_data.append({
            "Dataset": dataset,
            "Replicate": rep,
            "Method": method,
            "File": filename,
            "Total_Genes": n_genes,
            "Unique_Patterns": n_unique_patterns,
            "Distinct_Patterns": n_shared_patterns,
            "Max_Frequency": max_frequency,
            "Unique_Ratio": round(n_unique_patterns / n_genes, 3)
        })
        
        print(f"  ✓ Genes: {n_genes}, Únicos: {n_unique_patterns}, Patrones distintos: {n_shared_patterns}")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
    
    print()

# === Crear tabla resumen ===
summary_df = pd.DataFrame(summary_data)
# Ordenar por Dataset, Replicate y luego Method (alfabéticamente SSD antes que WCSS)
summary_df = summary_df.sort_values(by=["Dataset", "Replicate", "Method"])

# Guardar resumen
summary_df.to_csv(output_summary, sep="\t", index=False)
print(f"Resumen guardado en: {output_summary}")
print()
print(summary_df.to_string(index=False))

# === Estadísticas por método ===
print("\n" + "="*60)
print("ESTADÍSTICAS POR MÉTODO")
print("="*60)

for method in sorted(summary_df["Method"].unique()):
    method_data = summary_df[summary_df["Method"] == method]
    print(f"\n{method}:")
    print(f"  Archivos procesados: {len(method_data)}")
    print(f"  Promedio de patrones únicos: {method_data['Unique_Patterns'].mean():.1f}")
    print(f"  Promedio de patrones distintos: {method_data['Distinct_Patterns'].mean():.1f}")
    print(f"  Promedio ratio único: {method_data['Unique_Ratio'].mean():.3f}")
    print(f"  Frecuencia máxima promedio: {method_data['Max_Frequency'].mean():.1f}")

# === Comparación lado a lado ===
print("\n" + "="*60)
print("COMPARACIÓN DIRECTA (Promedio por Dataset-Replicate)")
print("="*60)

# Pivotar para comparación directa
comparison = summary_df.pivot_table(
    index=["Dataset", "Replicate"],
    columns="Method",
    values=["Unique_Patterns", "Distinct_Patterns", "Unique_Ratio"]
)
print("\n", comparison)