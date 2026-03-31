import os
import pandas as pd

# === Configuración ===
base_dir = "/home/lahumada/disco1/BNI3/Binarization/binary_expression_matrix/marin"
output_summary = os.path.join(base_dir, "summary_patterns.tsv")

# === Función para procesar un archivo ===
def process_binarized_file(input_file):
    """
    Procesa un archivo binarizado y retorna los patrones y sus frecuencias
    """
    df = pd.read_csv(input_file, sep="\t")

    # Detectar columna ID si existe
    gene_cols = [c for c in df.columns if c.lower() != "id"]

    patterns = {}
    for gene in gene_cols:
        sequence = ''.join(df[gene].astype(str).tolist())
        patterns[gene] = sequence

    pattern_counts = pd.Series(patterns).value_counts()

    output = pd.DataFrame([
        {
            "Gene": gene,
            "Pattern": pat,
            "Frequency": pattern_counts[pat]
        }
        for gene, pat in patterns.items()
    ])

    output = output.sort_values(
        by=["Frequency", "Pattern"],
        ascending=[True, True]
    )

    base, _ = os.path.splitext(input_file)
    output_file = f"{base}_pattern.tsv"
    output.to_csv(output_file, sep="\t", index=False)

    return output, os.path.basename(input_file)

# === Buscar archivos binarizados ===
all_files = []
for root, _, files in os.walk(base_dir):
    for file in files:
        if file.lower().startswith("bin_") and file.endswith(".tsv") and "_pattern" not in file:
            all_files.append(os.path.join(root, file))

print(f"Encontrados {len(all_files)} archivos binarizados\n")

# === Procesar archivos ===
summary_data = []

for file_path in sorted(all_files):
    filename = os.path.basename(file_path)
    print(f"Procesando: {filename}")

    try:
        patterns_df, fname = process_binarized_file(file_path)

        # === Parsing ROBUSTO del nombre ===
        # Acepta: bin_SSD.tsv, bin_WCSS.tsv
        name = fname.replace(".tsv", "")
        parts = name.split("_")

        method = parts[1] if len(parts) > 1 else "NA"
        size = "NA"
        dataset = "NA"
        rep = "NA"

        # === Estadísticas ===
        n_genes = len(patterns_df)
        n_unique_patterns = (patterns_df["Frequency"] == 1).sum()
        n_distinct_patterns = patterns_df["Pattern"].nunique()
        max_frequency = patterns_df["Frequency"].max()

        summary_data.append({
            "Dataset": dataset,
            "Replicate": rep,
            "Method": method,
            "File": fname,
            "Total_Genes": n_genes,
            "Unique_Patterns": n_unique_patterns,
            "Distinct_Patterns": n_distinct_patterns,
            "Max_Frequency": max_frequency,
            "Unique_Ratio": round(n_unique_patterns / n_genes, 3)
        })

        print(
            f"  ✓ Genes: {n_genes}, "
            f"Únicos: {n_unique_patterns}, "
            f"Patrones distintos: {n_distinct_patterns}"
        )

    except Exception as e:
        print(f"  ✗ Error procesando {filename}: {e}")

    print()

# === Crear resumen ===
summary_df = pd.DataFrame(summary_data)

if summary_df.empty:
    raise RuntimeError(
        "No se generó el resumen. "
        "Revisa el formato y contenido de los archivos binarizados."
    )

summary_df = summary_df.sort_values(
    by=["Dataset", "Replicate", "Method"]
)

summary_df.to_csv(output_summary, sep="\t", index=False)

print(f"Resumen guardado en: {output_summary}\n")
print(summary_df.to_string(index=False))

# === Estadísticas por método ===
print("\n" + "=" * 60)
print("ESTADÍSTICAS POR MÉTODO")
print("=" * 60)

for method in sorted(summary_df["Method"].unique()):
    md = summary_df[summary_df["Method"] == method]
    print(f"\n{method}:")
    print(f"  Archivos procesados: {len(md)}")
    print(f"  Promedio patrones únicos: {md['Unique_Patterns'].mean():.1f}")
    print(f"  Promedio patrones distintos: {md['Distinct_Patterns'].mean():.1f}")
    print(f"  Promedio ratio único: {md['Unique_Ratio'].mean():.3f}")
    print(f"  Frecuencia máxima promedio: {md['Max_Frequency'].mean():.1f}")

# === Comparación lado a lado ===
print("\n" + "=" * 60)
print("COMPARACIÓN DIRECTA (SSD vs WCSS)")
print("=" * 60)

comparison = summary_df.pivot_table(
    index="File",
    columns="Method",
    values=["Unique_Patterns", "Distinct_Patterns", "Unique_Ratio"]
)

print("\n", comparison)
