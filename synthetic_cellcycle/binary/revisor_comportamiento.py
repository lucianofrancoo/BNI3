import os
import pandas as pd

os.chdir("/home/lahumada/Desktop/BNI3/synthetic_cellcycle")

# === 1. Cargar matriz ===
input = "/home/lahumada/Desktop/BNI3/synthetic_cellcycle/binary/binary_data.tsv"
df = pd.read_csv(input, sep="\t")
df.head()

# === 2. Convertir columnas en secuencias temporales ===
patterns = {}
for gene in df.columns:
    sequence = ''.join(df[gene].astype(str).tolist())  # ej: "1110"
    patterns[gene] = sequence
print(patterns)

# === 3. Contar frecuencia de cada comportamiento ===
pattern_counts = pd.Series(patterns).value_counts() 
pattern_counts.head()

# === 4. Construir tabla de salida ===
output = pd.DataFrame([
    {"Gene": gene, 
     "Pattern": pat, 
     "Frequency": pattern_counts[pat]} 
    for gene, pat in patterns.items()
])
output.head()

# === 5. Ordenar: únicos arriba, compartidos abajo ===
output = output.sort_values(by=["Frequency", "Pattern"], ascending=[True, True])
print(output)

# === 6. Guardar a archivo ===
base, ext = os.path.splitext(input)
output_file = f"{base}_pattern.tsv"
output.to_csv(output_file, sep="\t", index=False)