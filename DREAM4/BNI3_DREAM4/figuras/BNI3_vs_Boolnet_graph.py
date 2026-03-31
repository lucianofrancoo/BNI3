#!/usr/bin/env python3
"""
Comparación estadística y visual entre BNI3 y BoolNet para datasets size10
Autor: Luciano Humada
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
from pathlib import Path

# ============================================================================
# CONFIGURACIÓN
# ============================================================================
BOOLNET_PATH = "/home/lahumada/disco1/BNI3/BoolNet_comparison/size10/BoolNet_complete_evaluation_results.tsv"
BNI3_PATH = "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/size10/rules/BNI3_size10_complete_evaluation_results.tsv"
OUTPUT_DIR = "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/figuras"

# Crear directorio de salida si no existe
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# ============================================================================
# CARGA Y PREPARACIÓN DE DATOS
# ============================================================================
print("Cargando datos...")
boolnet_df = pd.read_csv(BOOLNET_PATH, sep='\t')
bni3_df = pd.read_csv(BNI3_PATH, sep='\t')

# Preparar datos de BoolNet
boolnet_df['Method_Type'] = 'BoolNet'
boolnet_df['Condition'] = boolnet_df['Rep'].apply(
    lambda x: '5times' if '5times' in str(x) else 'Normal'
)
boolnet_df['Rep_Number'] = boolnet_df['Rep'].apply(
    lambda x: str(x).replace('_5times', '') if '5times' in str(x) else str(x)
)

# Preparar datos de BNI3
bni3_df['Method_Type'] = 'BNI3'
bni3_df['Condition'] = bni3_df['Rep'].apply(
    lambda x: '5times' if '5times' in str(x) else 'Normal'
)
bni3_df['Rep_Number'] = bni3_df['Rep'].apply(
    lambda x: str(x).replace('_5times', '') if '5times' in str(x) else str(x)
)

# Combinar datasets
combined_df = pd.concat([boolnet_df, bni3_df], ignore_index=True)

# Separar por condición
normal_df = combined_df[combined_df['Condition'] == 'Normal']
fivex_df = combined_df[combined_df['Condition'] == '5times']

print(f"\nData loaded:")
print(f"  BoolNet Normal (21 timepoints): {len(boolnet_df[boolnet_df['Condition'] == 'Normal'])} measurements")
print(f"  BoolNet 5 timepoints: {len(boolnet_df[boolnet_df['Condition'] == '5times'])} measurements")
print(f"  BNI3 Normal (21 timepoints): {len(bni3_df[bni3_df['Condition'] == 'Normal'])} measurements")
print(f"  BNI3 5 timepoints: {len(bni3_df[bni3_df['Condition'] == '5times'])} measurements")

# ============================================================================
# ESTADÍSTICAS DESCRIPTIVAS
# ============================================================================
print("\n" + "="*80)
print("DESCRIPTIVE STATISTICS")
print("="*80)

for condition in ['Normal', '5times']:
    condition_label = "21 Time Points" if condition == 'Normal' else "5 Time Points"
    print(f"\n{condition_label.upper()}:")
    condition_df = combined_df[combined_df['Condition'] == condition]
    
    for method in ['BNI3', 'BoolNet']:
        method_data = condition_df[condition_df['Method_Type'] == method]['Balanced_Accuracy']
        print(f"  {method}:")
        print(f"    Mean: {method_data.mean():.4f} ± {method_data.std():.4f}")
        print(f"    Median: {method_data.median():.4f}")
        print(f"    Range: [{method_data.min():.4f}, {method_data.max():.4f}]")
        print(f"    N: {len(method_data)}")

# ============================================================================
# PRUEBAS DE NORMALIDAD
# ============================================================================
print("\n" + "="*80)
print("NORMALITY TESTS (Shapiro-Wilk)")
print("="*80)

for condition in ['Normal', '5times']:
    condition_label = "21 Time Points" if condition == 'Normal' else "5 Time Points"
    print(f"\n{condition_label}:")
    for method in ['BNI3', 'BoolNet']:
        data = combined_df[
            (combined_df['Condition'] == condition) & 
            (combined_df['Method_Type'] == method)
        ]['Balanced_Accuracy']
        
        stat, p_value = stats.shapiro(data)
        normal = "YES" if p_value > 0.05 else "NO"
        print(f"  {method}: W={stat:.4f}, p={p_value:.4f} → Normal: {normal}")

# ============================================================================
# PRUEBAS ESTADÍSTICAS - WILCOXON PAREADO
# ============================================================================
print("\n" + "="*80)
print("STATISTICAL TESTS (Wilcoxon signed-rank test)")
print("="*80)
print("(Paired test for dependent samples)")

def perform_paired_test(df, condition_name):
    """Realiza test de Wilcoxon pareado entre BNI3 y BoolNet"""
    condition_data = df[df['Condition'] == condition_name].copy()
    
    # Crear identificador único para cada par Size-Rep
    condition_data['Pair_ID'] = condition_data['Size'] + '_' + condition_data['Rep_Number']
    
    # Pivotar para tener BNI3 y BoolNet en columnas
    pivot_df = condition_data.pivot_table(
        index='Pair_ID',
        columns='Method_Type',
        values='Balanced_Accuracy'
    )
    
    # Verificar que tenemos pares completos
    complete_pairs = pivot_df.dropna()
    n_pairs = len(complete_pairs)
    
    condition_label = "21 Time Points" if condition_name == 'Normal' else "5 Time Points"
    
    if n_pairs == 0:
        print(f"\n{condition_label}: No complete pairs to compare")
        return None
    
    bni3_values = complete_pairs['BNI3'].values
    boolnet_values = complete_pairs['BoolNet'].values
    
    # Test de Wilcoxon pareado
    stat, p_value = stats.wilcoxon(bni3_values, boolnet_values)
    
    # Calcular diferencias
    differences = bni3_values - boolnet_values
    mean_diff = differences.mean()
    median_diff = np.median(differences)
    
    # Determinar significancia
    significance = ""
    if p_value < 0.001:
        significance = "***"
    elif p_value < 0.01:
        significance = "**"
    elif p_value < 0.05:
        significance = "*"
    else:
        significance = "ns"
    
    print(f"\n{condition_label}:")
    print(f"  Complete pairs: {n_pairs}")
    print(f"  BNI3 mean: {bni3_values.mean():.4f}")
    print(f"  BoolNet mean: {boolnet_values.mean():.4f}")
    print(f"  Mean difference (BNI3 - BoolNet): {mean_diff:.4f}")
    print(f"  Median difference: {median_diff:.4f}")
    print(f"  Wilcoxon statistic: {stat:.2f}")
    print(f"  p-value: {p_value:.6f} {significance}")
    
    if p_value < 0.05:
        winner = "BNI3" if mean_diff > 0 else "BoolNet"
        print(f"  → {winner} is significantly better")
    else:
        print(f"  → No significant difference")
    
    return {
        'condition': condition_label,
        'n_pairs': n_pairs,
        'bni3_mean': bni3_values.mean(),
        'boolnet_mean': boolnet_values.mean(),
        'mean_diff': mean_diff,
        'median_diff': median_diff,
        'statistic': stat,
        'p_value': p_value,
        'significance': significance
    }

# Realizar pruebas para ambas condiciones
results_normal = perform_paired_test(combined_df, 'Normal')
results_5times = perform_paired_test(combined_df, '5times')

# ============================================================================
# VISUALIZACIÓN
# ============================================================================
print("\n" + "="*80)
print("GENERATING PLOTS")
print("="*80)

# Configurar estilo
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 11

# Crear figura con dos subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Colores personalizados
colors = {'BNI3': '#2E86AB', 'BoolNet': '#A23B72'}

def add_significance_bar(ax, x1, x2, y, p_value, height=0.02):
    """Añade barra de significancia al gráfico"""
    # Usar siempre un solo asterisco para simplificar visualmente
    if p_value < 0.05:
        sig = "*"
    else:
        sig = "ns"
    
    # Dibujar barra
    ax.plot([x1, x1, x2, x2], [y, y+height, y+height, y], 'k-', lw=1.5)
    ax.text((x1+x2)/2, y+height, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')

# Subplot 1: Datos Normales
ax1 = axes[0]
normal_plot_data = normal_df[['Method_Type', 'Balanced_Accuracy']].copy()
sns.boxplot(
    data=normal_plot_data,
    x='Method_Type',
    y='Balanced_Accuracy',
    hue='Method_Type',
    palette=colors,
    ax=ax1,
    width=0.5,
    legend=False
)
sns.stripplot(
    data=normal_plot_data,
    x='Method_Type',
    y='Balanced_Accuracy',
    color='black',
    alpha=0.3,
    size=4,
    ax=ax1
)

ax1.set_title('21 Time Points', fontsize=14, fontweight='bold')
ax1.set_xlabel('')  # Sin label
ax1.set_ylabel('StrAcc', fontsize=12, fontweight='bold')
ax1.set_ylim(0.45, 1.05)
# Tick labels en negrita
for label in ax1.get_xticklabels():
    label.set_fontweight('bold')
for label in ax1.get_yticklabels():
    label.set_fontweight('bold')

# Añadir estadísticas solo en el primer subplot
if results_normal:
    y_max = normal_plot_data['Balanced_Accuracy'].max()
    add_significance_bar(ax1, 0, 1, y_max + 0.05, results_normal['p_value'])

# Subplot 2: Datos 5times
ax2 = axes[1]
fivex_plot_data = fivex_df[['Method_Type', 'Balanced_Accuracy']].copy()
sns.boxplot(
    data=fivex_plot_data,
    x='Method_Type',
    y='Balanced_Accuracy',
    hue='Method_Type',
    palette=colors,
    ax=ax2,
    width=0.5,
    legend=False
)
sns.stripplot(
    data=fivex_plot_data,
    x='Method_Type',
    y='Balanced_Accuracy',
    color='black',
    alpha=0.3,
    size=4,
    ax=ax2
)

ax2.set_title('5 Time Points', fontsize=14, fontweight='bold')
ax2.set_xlabel('')  # Sin label
ax2.set_ylabel('StrAcc', fontsize=12, fontweight='bold')
ax2.set_ylim(0.45, 1.05)
# Tick labels en negrita
for label in ax2.get_xticklabels():
    label.set_fontweight('bold')
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')

# Añadir barra de significancia
if results_5times:
    y_max = fivex_plot_data['Balanced_Accuracy'].max()
    add_significance_bar(ax2, 0, 1, y_max + 0.05, results_5times['p_value'])

# Ajustar layout
plt.tight_layout()

# Guardar figura
output_path = Path(OUTPUT_DIR) / "BNI3_vs_BoolNet_Balanced_Accuracy_comparison.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nFigure saved at: {output_path}")

# También guardar en formato SVG para publicación
output_path_svg = Path(OUTPUT_DIR) / "BNI3_vs_BoolNet_Balanced_Accuracy_comparison.svg"
plt.savefig(output_path_svg, format='svg', bbox_inches='tight')
print(f"SVG figure saved at: {output_path_svg}")

plt.close()

# ============================================================================
# GUARDAR RESULTADOS ESTADÍSTICOS
# ============================================================================
results_summary = pd.DataFrame([results_normal, results_5times])
# Transponer la tabla
results_summary_transposed = results_summary.set_index('condition').T
results_path = Path(OUTPUT_DIR) / "statistical_comparison_results.tsv"
results_summary_transposed.to_csv(results_path, sep='\t')
print(f"Statistical results saved at: {results_path}")

print("\n" + "="*80)
print("ANALYSIS COMPLETED")
print("="*80)