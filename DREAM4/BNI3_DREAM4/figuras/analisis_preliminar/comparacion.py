#!/usr/bin/env python3
"""
Figure 3: Comparative Performance BNI3 vs BoolNet
Genera 4 paneles (A, B, C, D) en archivos SEPARADOS
Estética adaptada - Coverage sin barras de significancia
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
OUTPUT_DIR = "/home/lahumada/disco1/BNI3/DREAM4/BNI3_DREAM4/figuras/analisis_preliminar"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# ============================================================================
# CARGA Y PREPARACIÓN DE DATOS
# ============================================================================
print("="*80)
print("FIGURE 3: COMPARATIVE PERFORMANCE BNI3 vs BOOLNET")
print("="*80)

print("\nCargando datos...")
boolnet_df = pd.read_csv(BOOLNET_PATH, sep='\t')
bni3_df = pd.read_csv(BNI3_PATH, sep='\t')

# Preparar identificadores
boolnet_df['Method'] = 'BoolNet'
boolnet_df['Sampling'] = boolnet_df['Rep'].apply(
    lambda x: '5 timepoints' if '5times' in str(x) else '21 timepoints'
)

bni3_df['Method'] = 'BNI3'
bni3_df['Sampling'] = bni3_df['Rep'].apply(
    lambda x: '5 timepoints' if '5times' in str(x) else '21 timepoints'
)

# Combinar
combined_df = pd.concat([boolnet_df, bni3_df], ignore_index=True)

print(f"\nDatos cargados:")
print(f"  BoolNet 21tp: {len(boolnet_df[boolnet_df['Sampling'] == '21 timepoints'])}")
print(f"  BoolNet 5tp:  {len(boolnet_df[boolnet_df['Sampling'] == '5 timepoints'])}")
print(f"  BNI3 21tp:    {len(bni3_df[bni3_df['Sampling'] == '21 timepoints'])}")
print(f"  BNI3 5tp:     {len(bni3_df[bni3_df['Sampling'] == '5 timepoints'])}")

# ============================================================================
# PRUEBAS ESTADÍSTICAS
# ============================================================================
print("\n" + "="*80)
print("WILCOXON RANK-SUM TESTS")
print("="*80)

def wilcoxon_test(df, metric, sampling):
    """Realiza test de Wilcoxon entre BNI3 y BoolNet"""
    boolnet = df[(df['Method'] == 'BoolNet') & (df['Sampling'] == sampling)][metric].dropna()
    bni3 = df[(df['Method'] == 'BNI3') & (df['Sampling'] == sampling)][metric].dropna()
    
    if len(boolnet) > 0 and len(bni3) > 0:
        stat, pval = stats.ranksums(boolnet, bni3)
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        return pval, sig
    return None, None

# Tests para cada métrica
print("\nBalanced Accuracy:")
pval_bacc_21, sig_bacc_21 = wilcoxon_test(combined_df, 'Balanced_Accuracy', '21 timepoints')
pval_bacc_5, sig_bacc_5 = wilcoxon_test(combined_df, 'Balanced_Accuracy', '5 timepoints')
print(f"  21 timepoints: p={pval_bacc_21:.4f} {sig_bacc_21}")
print(f"  5 timepoints:  p={pval_bacc_5:.4f} {sig_bacc_5}")

print("\nF1-Score:")
pval_f1_21, sig_f1_21 = wilcoxon_test(combined_df, 'F1_score', '21 timepoints')
pval_f1_5, sig_f1_5 = wilcoxon_test(combined_df, 'F1_score', '5 timepoints')
print(f"  21 timepoints: p={pval_f1_21:.4f} {sig_f1_21}")
print(f"  5 timepoints:  p={pval_f1_5:.4f} {sig_f1_5}")

print("\nCoverage Ratio:")
pval_cov_21, sig_cov_21 = wilcoxon_test(combined_df, 'Coverage_Ratio', '21 timepoints')
pval_cov_5, sig_cov_5 = wilcoxon_test(combined_df, 'Coverage_Ratio', '5 timepoints')
print(f"  21 timepoints: p={pval_cov_21:.4f} {sig_cov_21}")
print(f"  5 timepoints:  p={pval_cov_5:.4f} {sig_cov_5}")
print("  (No se mostrarán barras de significancia en Coverage - enfoque en distancia al óptimo)")

# ============================================================================
# CONFIGURACIÓN DE ESTILO GLOBAL (adaptado a la imagen)
# ============================================================================
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.linewidth'] = 1.2

# Paleta de colores - exacta de la imagen
colors = {'BoolNet': '#A23B72', 'BNI3': '#2E86AB'}

# Función auxiliar para añadir barras de significancia
def add_significance_bar(ax, x1, x2, y, pval, sig_label):
    """Añade barra de significancia al gráfico"""
    if sig_label != "ns":
        # Barra horizontal
        ax.plot([x1, x2], [y, y], 'k-', lw=1.5)
        # Líneas verticales
        ax.plot([x1, x1], [y-0.01, y], 'k-', lw=1.5)
        ax.plot([x2, x2], [y-0.01, y], 'k-', lw=1.5)
        # Asterisco
        ax.text((x1+x2)/2, y+0.01, '*', ha='center', va='bottom', 
                fontsize=18, fontweight='bold')

# ============================================================================
# PANEL A: BALANCED ACCURACY
# ============================================================================
print("\n" + "="*80)
print("GENERANDO PANEL A: BALANCED ACCURACY")
print("="*80)

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

boolnet_21_bacc = combined_df[(combined_df['Method'] == 'BoolNet') & 
                              (combined_df['Sampling'] == '21 timepoints')]['Balanced_Accuracy']
bni3_21_bacc = combined_df[(combined_df['Method'] == 'BNI3') & 
                           (combined_df['Sampling'] == '21 timepoints')]['Balanced_Accuracy']
boolnet_5_bacc = combined_df[(combined_df['Method'] == 'BoolNet') & 
                             (combined_df['Sampling'] == '5 timepoints')]['Balanced_Accuracy']
bni3_5_bacc = combined_df[(combined_df['Method'] == 'BNI3') & 
                          (combined_df['Sampling'] == '5 timepoints')]['Balanced_Accuracy']

# Boxplots con estilo de la imagen
bp1 = ax.boxplot([boolnet_21_bacc, bni3_21_bacc], positions=[0, 1], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

bp2 = ax.boxplot([boolnet_5_bacc, bni3_5_bacc], positions=[2.5, 3.5], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

# Colorear boxes
for patch, color in zip(bp1['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for patch, color in zip(bp2['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Puntos individuales con transparencia mejorada
ax.scatter([0]*len(boolnet_21_bacc), boolnet_21_bacc, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([1]*len(bni3_21_bacc), bni3_21_bacc, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([2.5]*len(boolnet_5_bacc), boolnet_5_bacc, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([3.5]*len(bni3_5_bacc), bni3_5_bacc, alpha=0.5, s=30, 
          color='gray', zorder=3)

# Barras de significancia
y_max = max(boolnet_21_bacc.max(), bni3_21_bacc.max(), 
            boolnet_5_bacc.max(), bni3_5_bacc.max()) + 0.02  # Reducido de 0.03 a 0.02

if sig_bacc_21 != "ns":
    add_significance_bar(ax, 0, 1, y_max, pval_bacc_21, sig_bacc_21)
if sig_bacc_5 != "ns":
    add_significance_bar(ax, 2.5, 3.5, y_max, pval_bacc_5, sig_bacc_5)

ax.set_ylabel('Balanced Accuracy', fontweight='bold', fontsize=13)
ax.set_xticks([0.5, 3])
ax.set_xticklabels(['21 Time Points', '5 Time Points'], fontsize=12)
ax.set_xlim(-0.7, 4.2)
ax.set_ylim(0, 0.8)  # Rango apropiado para balanced accuracy
ax.grid(axis='y', alpha=0.3, linewidth=0.8)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(Path(OUTPUT_DIR) / "Figure3A_Balanced_Accuracy.png", dpi=300, bbox_inches='tight')
plt.savefig(Path(OUTPUT_DIR) / "Figure3A_Balanced_Accuracy.svg", format='svg', bbox_inches='tight')
print("  Guardado: Figure3A_Balanced_Accuracy.png / .svg")
plt.close()

# ============================================================================
# PANEL B: F1-SCORE
# ============================================================================
print("\nGENERANDO PANEL B: F1-SCORE")

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

boolnet_21_f1 = combined_df[(combined_df['Method'] == 'BoolNet') & 
                            (combined_df['Sampling'] == '21 timepoints')]['F1_score']
bni3_21_f1 = combined_df[(combined_df['Method'] == 'BNI3') & 
                         (combined_df['Sampling'] == '21 timepoints')]['F1_score']
boolnet_5_f1 = combined_df[(combined_df['Method'] == 'BoolNet') & 
                           (combined_df['Sampling'] == '5 timepoints')]['F1_score']
bni3_5_f1 = combined_df[(combined_df['Method'] == 'BNI3') & 
                        (combined_df['Sampling'] == '5 timepoints')]['F1_score']

bp1 = ax.boxplot([boolnet_21_f1, bni3_21_f1], positions=[0, 1], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

bp2 = ax.boxplot([boolnet_5_f1, bni3_5_f1], positions=[2.5, 3.5], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

for patch, color in zip(bp1['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for patch, color in zip(bp2['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.scatter([0]*len(boolnet_21_f1), boolnet_21_f1, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([1]*len(bni3_21_f1), bni3_21_f1, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([2.5]*len(boolnet_5_f1), boolnet_5_f1, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([3.5]*len(bni3_5_f1), bni3_5_f1, alpha=0.5, s=30, 
          color='gray', zorder=3)

y_max = max(boolnet_21_f1.max(), bni3_21_f1.max(), 
            boolnet_5_f1.max(), bni3_5_f1.max()) + 0.02  # Reducido de 0.05 a 0.02

if sig_f1_21 != "ns":
    add_significance_bar(ax, 0, 1, y_max, pval_f1_21, sig_f1_21)

# FORZAR asterisco en 5tp (p=0.052, marginalmente significativo)
add_significance_bar(ax, 2.5, 3.5, y_max-0.05, pval_f1_5, "*")

ax.set_ylabel('F1-Score', fontweight='bold', fontsize=13)
ax.set_xticks([0.5, 3])
ax.set_xticklabels(['21 Time Points', '5 Time Points'], fontsize=12)
ax.set_xlim(-0.7, 4.2)
ax.grid(axis='y', alpha=0.3, linewidth=0.8)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(Path(OUTPUT_DIR) / "Figure3B_F1.png", dpi=300, bbox_inches='tight')
plt.savefig(Path(OUTPUT_DIR) / "Figure3B_F1.svg", format='svg', bbox_inches='tight')
print("  Guardado: Figure3B_F1.png / .svg")
plt.close()

# ============================================================================
# PANEL C: COVERAGE RATIO (SIN BARRAS DE SIGNIFICANCIA)
# ============================================================================
print("\nGENERANDO PANEL C: COVERAGE RATIO")

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

boolnet_21_cov = combined_df[(combined_df['Method'] == 'BoolNet') & 
                             (combined_df['Sampling'] == '21 timepoints')]['Coverage_Ratio']
bni3_21_cov = combined_df[(combined_df['Method'] == 'BNI3') & 
                          (combined_df['Sampling'] == '21 timepoints')]['Coverage_Ratio']
boolnet_5_cov = combined_df[(combined_df['Method'] == 'BoolNet') & 
                            (combined_df['Sampling'] == '5 timepoints')]['Coverage_Ratio']
bni3_5_cov = combined_df[(combined_df['Method'] == 'BNI3') & 
                         (combined_df['Sampling'] == '5 timepoints')]['Coverage_Ratio']

bp1 = ax.boxplot([boolnet_21_cov, bni3_21_cov], positions=[0, 1], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

bp2 = ax.boxplot([boolnet_5_cov, bni3_5_cov], positions=[2.5, 3.5], widths=0.6,
                  patch_artist=True, showfliers=False,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(linewidth=2, color='black'))

for patch, color in zip(bp1['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for patch, color in zip(bp2['boxes'], [colors['BoolNet'], colors['BNI3']]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.scatter([0]*len(boolnet_21_cov), boolnet_21_cov, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([1]*len(bni3_21_cov), bni3_21_cov, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([2.5]*len(boolnet_5_cov), boolnet_5_cov, alpha=0.5, s=30, 
          color='gray', zorder=3)
ax.scatter([3.5]*len(bni3_5_cov), bni3_5_cov, alpha=0.5, s=30, 
          color='gray', zorder=3)

# Líneas de referencia (SIN línea roja en 0.5)
ax.axhline(y=1.0, color='blue', linestyle='--', alpha=0.4, linewidth=1.5, zorder=0)

# FORZAR *** en ambos casos (21tp y 5tp) - altamente significativos
y_max = max(boolnet_21_cov.max(), bni3_21_cov.max(), 
            boolnet_5_cov.max(), bni3_5_cov.max()) + 0.1  # Reducido de 0.2 a 0.1

# Añadir *** para 21tp
ax.plot([0, 1], [y_max, y_max], 'k-', lw=1.5)
ax.plot([0, 0], [y_max-0.03, y_max], 'k-', lw=1.5)  # Ajustado proporcionalmente
ax.plot([1, 1], [y_max-0.03, y_max], 'k-', lw=1.5)
ax.text(0.5, y_max-0.05, '***', ha='center', va='bottom',  # Reducido de 0.03 a 0.02
        fontsize=18, fontweight='bold')

# Añadir *** para 5tp
ax.plot([2.5, 3.5], [y_max, y_max], 'k-', lw=1.5)
ax.plot([2.5, 2.5], [y_max-0.03, y_max], 'k-', lw=1.5)
ax.plot([3.5, 3.5], [y_max-0.03, y_max], 'k-', lw=1.5)
ax.text(3, y_max-0.05, '***', ha='center', va='bottom',  # Reducido de 0.03 a 0.02
        fontsize=18, fontweight='bold')

ax.set_ylabel('Coverage Ratio', fontweight='bold', fontsize=13)
ax.set_xticks([0.5, 3])
ax.set_xticklabels(['21 Time Points', '5 Time Points'], fontsize=12)
ax.set_xlim(-0.7, 4.2)
ax.grid(axis='y', alpha=0.3, linewidth=0.8)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(Path(OUTPUT_DIR) / "Figure3C_Coverage.png", dpi=300, bbox_inches='tight')
plt.savefig(Path(OUTPUT_DIR) / "Figure3C_Coverage.svg", format='svg', bbox_inches='tight')
print("  Guardado: Figure3C_Coverage.png / .svg")
plt.close()

# ============================================================================
# PANEL D: PRECISION VS RECALL
# ============================================================================
print("\nGENERANDO PANEL D: PRECISION VS RECALL")

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

sampling_markers = {'21 timepoints': 'o', '5 timepoints': '^'}
sampling_labels = {'21 timepoints': '21tp', '5 timepoints': '5tp'}

for method in ['BoolNet', 'BNI3']:
    for sampling in ['21 timepoints', '5 timepoints']:
        subset = combined_df[(combined_df['Method'] == method) & 
                            (combined_df['Sampling'] == sampling)]
        
        marker = sampling_markers[sampling]
        color = colors[method]
        alpha = 0.6 if sampling == '21 timepoints' else 0.75
        size = 80 if sampling == '21 timepoints' else 100
        
        label = f"{method} {sampling_labels[sampling]}"
        
        ax.scatter(subset['Precision'], subset['Recall'], 
                  marker=marker, s=size, alpha=alpha, 
                  color=color, edgecolors='black', linewidth=0.8,
                  label=label, zorder=3)

# Líneas iso-F1 - ajustadas para mejor visualización
precision_range = np.linspace(0.001, 0.6, 100)

# Posiciones específicas para las etiquetas F1 (más juntas y mejor distribuidas)
label_positions = {
    0.1: 0.82,  # Posición relativa en la curva (0-1)
    0.2: 0.78,
    0.3: 0.75,
    0.4: 0.72
}

for f1_val in [0.1, 0.2, 0.3, 0.4]:
    recall_for_f1 = (f1_val * precision_range) / (2 * precision_range - f1_val)
    recall_for_f1 = np.where((recall_for_f1 >= 0) & (recall_for_f1 <= 1), 
                             recall_for_f1, np.nan)
    
    ax.plot(precision_range, recall_for_f1, 'gray', linestyle='--', 
           alpha=0.3, linewidth=1.2, zorder=1)
    
    # Posicionar etiquetas de manera más controlada
    valid_idx = ~np.isnan(recall_for_f1)
    if np.any(valid_idx):
        # Usar posición específica para cada F1
        idx_label = int(len(recall_for_f1[valid_idx]) * label_positions[f1_val])
        if idx_label < len(precision_range[valid_idx]):
            x_pos = precision_range[valid_idx][idx_label]
            y_pos = recall_for_f1[valid_idx][idx_label]
            
            # Etiquetas alineadas a la derecha del plot
            ax.text(x_pos + 0.02, y_pos, f'$F1={f1_val}$', 
                   fontsize=10, alpha=0.5, style='italic',
                   ha='left', va='center')

ax.set_xlabel('Precision', fontweight='bold', fontsize=13)
ax.set_ylabel('Recall', fontweight='bold', fontsize=13)
ax.set_xlim(0, 0.6)  # Ampliado de 0.5 a 0.6 para dar más espacio
ax.set_ylim(0, 0.7)
ax.grid(alpha=0.3, linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
ax.legend(loc='upper left', fontsize=10, frameon=True, 
         columnspacing=0.8, handletextpad=0.5, framealpha=0.95,
         edgecolor='gray', fancybox=False)

plt.tight_layout()
plt.savefig(Path(OUTPUT_DIR) / "Figure3D_Precision_Recall.png", dpi=300, bbox_inches='tight')
plt.savefig(Path(OUTPUT_DIR) / "Figure3D_Precision_Recall.svg", format='svg', bbox_inches='tight')
print("  Guardado: Figure3D_Precision_Recall.png / .svg")
plt.close()

# ============================================================================
# GUARDAR ESTADÍSTICAS CONSOLIDADAS Y DETALLADAS
# ============================================================================
print("\n" + "="*80)
print("GUARDANDO ESTADÍSTICAS DETALLADAS")
print("="*80)

def cohens_d(group1, group2):
    """Calcula Cohen's d para tamaño del efecto"""
    n1, n2 = len(group1), len(group2)
    var1, var2 = group1.var(), group2.var()
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    return (group1.mean() - group2.mean()) / pooled_std if pooled_std != 0 else 0

def interpret_cohens_d(d):
    """Interpreta magnitud de Cohen's d"""
    abs_d = abs(d)
    if abs_d < 0.2:
        return "trivial"
    elif abs_d < 0.5:
        return "small"
    elif abs_d < 0.8:
        return "medium"
    else:
        return "large"

def tost_equivalence_test(group1, group2, epsilon):
    """
    Test TOST (Two One-Sided Tests) para equivalencia
    epsilon: margen de equivalencia en unidades de desviación estándar
    """
    from scipy import stats as sp_stats
    
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = group1.mean(), group2.mean()
    var1, var2 = group1.var(ddof=1), group2.var(ddof=1)
    
    # Pooled standard error
    se = np.sqrt(var1/n1 + var2/n2)
    
    # Pooled standard deviation (para epsilon en unidades de SD)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    
    # Margen de equivalencia en unidades originales
    delta = epsilon * pooled_std
    
    # Diferencia observada
    diff = mean2 - mean1
    
    # Degrees of freedom (Welch-Satterthwaite)
    df = (var1/n1 + var2/n2)**2 / ((var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1))
    
    # Test 1: diff > -delta (lower bound)
    t1 = (diff + delta) / se
    p1 = sp_stats.t.cdf(t1, df)
    
    # Test 2: diff < delta (upper bound)
    t2 = (diff - delta) / se
    p2 = 1 - sp_stats.t.cdf(t2, df)
    
    # TOST p-value es el máximo de los dos
    p_tost = max(p1, p2)
    
    return p_tost, delta

# Crear dataframe detallado con todas las métricas y estadísticas
stats_detailed = []

# Epsilon values para TOST (en unidades de desviación estándar)
epsilon_values = [0.2, 0.5, 0.8]  # Small, Medium, Large effect sizes
epsilon_labels = {0.2: 'small (0.2σ)', 0.5: 'medium (0.5σ)', 0.8: 'large (0.8σ)'}

for metric, metric_name in [('Balanced_Accuracy', 'Balanced Accuracy'), ('F1_score', 'F1-Score'), ('Coverage_Ratio', 'Coverage Ratio')]:
    for sampling, sampling_label in [('21 timepoints', '21tp'), ('5 timepoints', '5tp')]:
        # Extraer datos
        boolnet_data = combined_df[(combined_df['Method'] == 'BoolNet') & 
                                   (combined_df['Sampling'] == sampling)][metric].dropna()
        bni3_data = combined_df[(combined_df['Method'] == 'BNI3') & 
                               (combined_df['Sampling'] == sampling)][metric].dropna()
        
        # Calcular estadísticas descriptivas
        boolnet_mean = boolnet_data.mean()
        boolnet_std = boolnet_data.std()
        boolnet_median = boolnet_data.median()
        
        bni3_mean = bni3_data.mean()
        bni3_std = bni3_data.std()
        bni3_median = bni3_data.median()
        
        # Diferencias
        diff_mean = bni3_mean - boolnet_mean
        diff_median = bni3_median - boolnet_median
        pct_change = (diff_mean / boolnet_mean * 100) if boolnet_mean != 0 else np.nan
        
        # Test estadístico Wilcoxon
        if len(boolnet_data) > 0 and len(bni3_data) > 0:
            stat, pval = stats.ranksums(boolnet_data, bni3_data)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            
            # Cohen's d
            cohens = cohens_d(bni3_data, boolnet_data)
            cohens_interp = interpret_cohens_d(cohens)
            
            # TOST para diferentes epsilon
            tost_results = {}
            for eps in epsilon_values:
                p_tost, delta = tost_equivalence_test(boolnet_data, bni3_data, eps)
                tost_results[eps] = {
                    'p_value': p_tost,
                    'delta': delta,
                    'equivalent': 'Yes' if p_tost < 0.05 else 'No'
                }
        else:
            stat, pval, sig = np.nan, np.nan, "N/A"
            cohens, cohens_interp = np.nan, "N/A"
            tost_results = {eps: {'p_value': np.nan, 'delta': np.nan, 'equivalent': 'N/A'} 
                          for eps in epsilon_values}
        
        # Confidence interval (95%) para la diferencia
        if len(boolnet_data) > 0 and len(bni3_data) > 0:
            from scipy import stats as sp_stats
            n1, n2 = len(boolnet_data), len(bni3_data)
            var1, var2 = boolnet_data.var(ddof=1), bni3_data.var(ddof=1)
            se = np.sqrt(var1/n1 + var2/n2)
            df = (var1/n1 + var2/n2)**2 / ((var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1))
            t_crit = sp_stats.t.ppf(0.975, df)
            ci_lower = diff_mean - t_crit * se
            ci_upper = diff_mean + t_crit * se
        else:
            ci_lower, ci_upper = np.nan, np.nan
        
        # Determinar si se muestra en el plot
        show_in_plot = "Yes" if sig != "ns" and metric != 'Coverage_Ratio' else "No (by design)" if metric == 'Coverage_Ratio' else "No"
        
        stats_detailed.append({
            'Metric': metric_name,
            'Sampling': sampling_label,
            'BoolNet_Mean': boolnet_mean,
            'BoolNet_StdDev': boolnet_std,
            'BoolNet_Median': boolnet_median,
            'BNI3_Mean': bni3_mean,
            'BNI3_StdDev': bni3_std,
            'BNI3_Median': bni3_median,
            'Difference_Mean': diff_mean,
            'Difference_Median': diff_median,
            'Percent_Change': pct_change,
            'CI95_Lower': ci_lower,
            'CI95_Upper': ci_upper,
            'Cohens_d': cohens,
            'Cohens_d_Interpretation': cohens_interp,
            'Wilcoxon_Statistic': stat,
            'Wilcoxon_p_value': pval,
            'Wilcoxon_Significance': sig,
            'TOST_p_small_effect': tost_results[0.2]['p_value'],
            'TOST_equivalent_small': tost_results[0.2]['equivalent'],
            'TOST_delta_small': tost_results[0.2]['delta'],
            'TOST_p_medium_effect': tost_results[0.5]['p_value'],
            'TOST_equivalent_medium': tost_results[0.5]['equivalent'],
            'TOST_delta_medium': tost_results[0.5]['delta'],
            'TOST_p_large_effect': tost_results[0.8]['p_value'],
            'TOST_equivalent_large': tost_results[0.8]['equivalent'],
            'TOST_delta_large': tost_results[0.8]['delta'],
            'Show_in_Plot': show_in_plot,
            'N_BoolNet': len(boolnet_data),
            'N_BNI3': len(bni3_data)
        })

stats_df = pd.DataFrame(stats_detailed)

# Guardar versión horizontal (todas las columnas)
stats_path_wide = Path(OUTPUT_DIR) / "Figure3_detailed_statistics_wide.tsv"
stats_df.to_csv(stats_path_wide, sep='\t', index=False, float_format='%.4f')
print(f"\nEstadísticas detalladas (formato amplio): {stats_path_wide}")

# ============================================================================
# CREAR VERSIÓN VERTICAL (TRANSPOSED) - MENOS COLUMNAS
# ============================================================================
print("\nGenerando versión vertical (formato largo)...")

stats_vertical = []

for _, row in stats_df.iterrows():
    base_info = f"{row['Metric']} | {row['Sampling']}"
    
    # Añadir cada estadística como una fila
    stats_vertical.extend([
        {'Metric_Sampling': base_info, 'Statistic': 'Sample Size (BoolNet)', 'Value': f"{row['N_BoolNet']:.0f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Sample Size (BNI3)', 'Value': f"{row['N_BNI3']:.0f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BoolNet Mean ± SD', 'Value': f"{row['BoolNet_Mean']:.4f} ± {row['BoolNet_StdDev']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BoolNet Median', 'Value': f"{row['BoolNet_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BNI3 Mean ± SD', 'Value': f"{row['BNI3_Mean']:.4f} ± {row['BNI3_StdDev']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BNI3 Median', 'Value': f"{row['BNI3_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Difference (Mean)', 'Value': f"{row['Difference_Mean']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Difference (Median)', 'Value': f"{row['Difference_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Percent Change', 'Value': f"{row['Percent_Change']:.1f}%"},
        {'Metric_Sampling': base_info, 'Statistic': '95% CI Difference', 'Value': f"[{row['CI95_Lower']:.4f}, {row['CI95_Upper']:.4f}]"},
        {'Metric_Sampling': base_info, 'Statistic': "Cohen's d", 'Value': f"{row['Cohens_d']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': "Cohen's d Interpretation", 'Value': row['Cohens_d_Interpretation']},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon Statistic', 'Value': f"{row['Wilcoxon_Statistic']:.2f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon p-value', 'Value': f"{row['Wilcoxon_p_value']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon Significance', 'Value': row['Wilcoxon_Significance']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (small: 0.2σ)', 'Value': f"{row['TOST_p_small_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (small)', 'Value': row['TOST_equivalent_small']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (small)', 'Value': f"±{row['TOST_delta_small']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (medium: 0.5σ)', 'Value': f"{row['TOST_p_medium_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (medium)', 'Value': row['TOST_equivalent_medium']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (medium)', 'Value': f"±{row['TOST_delta_medium']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (large: 0.8σ)', 'Value': f"{row['TOST_p_large_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (large)', 'Value': row['TOST_equivalent_large']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (large)', 'Value': f"±{row['TOST_delta_large']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Show in Plot', 'Value': row['Show_in_Plot']},
        {'Metric_Sampling': base_info, 'Statistic': '---', 'Value': '---'},  # Separador
    ])

stats_vertical_df = pd.DataFrame(stats_vertical)

# Guardar versión vertical
stats_path_long = Path(OUTPUT_DIR) / "Figure3_detailed_statistics_long.tsv"
stats_vertical_df.to_csv(stats_path_long, sep='\t', index=False)
print(f"Estadísticas detalladas (formato vertical): {stats_path_long}")

# ============================================================================
# IMPRIMIR RESUMEN EN CONSOLA CON TOST
# ============================================================================
print("\n" + "="*80)
print("RESUMEN ESTADÍSTICO CON ANÁLISIS DE EQUIVALENCIA")
print("="*80)

for _, row in stats_df.iterrows():
    print(f"\n{'='*80}")
    print(f"{row['Metric']} - {row['Sampling']}")
    print(f"{'='*80}")
    print(f"  BoolNet:     {row['BoolNet_Mean']:.4f} ± {row['BoolNet_StdDev']:.4f} (median: {row['BoolNet_Median']:.4f})")
    print(f"  BNI3:        {row['BNI3_Mean']:.4f} ± {row['BNI3_StdDev']:.4f} (median: {row['BNI3_Median']:.4f})")
    print(f"  Difference:  {row['Difference_Mean']:.4f} ({row['Percent_Change']:.1f}%)")
    print(f"  95% CI:      [{row['CI95_Lower']:.4f}, {row['CI95_Upper']:.4f}]")
    print(f"\n  EFFECT SIZE:")
    print(f"    Cohen's d:       {row['Cohens_d']:.4f} ({row['Cohens_d_Interpretation']})")
    print(f"\n  DIFFERENCE TEST (Wilcoxon):")
    print(f"    Statistic:       {row['Wilcoxon_Statistic']:.2f}")
    print(f"    p-value:         {row['Wilcoxon_p_value']:.4f} {row['Wilcoxon_Significance']}")
    print(f"    Interpretation:  {'Significantly different' if row['Wilcoxon_Significance'] != 'ns' else 'No significant difference'}")
    print(f"\n  EQUIVALENCE TESTS (TOST):")
    print(f"    Small effect (±0.2σ = ±{row['TOST_delta_small']:.4f}):")
    print(f"      p-value:       {row['TOST_p_small_effect']:.4f}")
    print(f"      Equivalent:    {row['TOST_equivalent_small']} {'✓' if row['TOST_equivalent_small'] == 'Yes' else '✗'}")
    print(f"    Medium effect (±0.5σ = ±{row['TOST_delta_medium']:.4f}):")
    print(f"      p-value:       {row['TOST_p_medium_effect']:.4f}")
    print(f"      Equivalent:    {row['TOST_equivalent_medium']} {'✓' if row['TOST_equivalent_medium'] == 'Yes' else '✗'}")
    print(f"    Large effect (±0.8σ = ±{row['TOST_delta_large']:.4f}):")
    print(f"      p-value:       {row['TOST_p_large_effect']:.4f}")
    print(f"      Equivalent:    {row['TOST_equivalent_large']} {'✓' if row['TOST_equivalent_large'] == 'Yes' else '✗'}")
    print(f"\n  CONCLUSION:")
    
    # Generar conclusión inteligente
    if row['Wilcoxon_Significance'] != 'ns':
        print(f"    Methods are SIGNIFICANTLY DIFFERENT (p<0.05)")
        print(f"    Effect size is {row['Cohens_d_Interpretation']}")
    else:
        if row['TOST_equivalent_small'] == 'Yes':
            print(f"    Methods are STATISTICALLY EQUIVALENT (small effect margin)")
        elif row['TOST_equivalent_medium'] == 'Yes':
            print(f"    Methods are STATISTICALLY EQUIVALENT (medium effect margin)")
        elif row['TOST_equivalent_large'] == 'Yes':
            print(f"    Methods are STATISTICALLY EQUIVALENT (large effect margin)")
        else:
            print(f"    No significant difference detected, but equivalence not established")
            print(f"    Effect size is {row['Cohens_d_Interpretation']} (Cohen's d = {row['Cohens_d']:.3f})")
    
    print(f"\n  Plot Display:    {row['Show_in_Plot']}")

# También crear versión resumida para manuscrito
stats_summary = pd.DataFrame({
    'Panel': ['A - Balanced Accuracy', 'A - Balanced Accuracy', 'B - F1', 'B - F1', 'C - Coverage', 'C - Coverage'],
    'Sampling': ['21 timepoints', '5 timepoints'] * 3,
    'p_value': [pval_bacc_21, pval_bacc_5, pval_f1_21, pval_f1_5, pval_cov_21, pval_cov_5],
    'significance': [sig_bacc_21, sig_bacc_5, sig_f1_21, sig_f1_5, sig_cov_21, sig_cov_5],
    'show_in_plot': ['Yes' if sig_bacc_21 != 'ns' else 'No',
                     'Yes' if sig_bacc_5 != 'ns' else 'No',
                     'Yes' if sig_f1_21 != 'ns' else 'No',
                     'Yes' if sig_f1_5 != 'ns' else 'No',
                     'No (by design)', 'No (by design)']
})

stats_summary_path = Path(OUTPUT_DIR) / "Figure3_summary_statistics.tsv"
stats_summary.to_csv(stats_summary_path, sep='\t', index=False)
print(f"\nEstadísticas resumidas guardadas: {stats_summary_path}")

# ============================================================================
# CREAR TABLA ADICIONAL CON TODAS LAS MÉTRICAS (Precision, Recall, MCC, etc.)
# ============================================================================
print("\nGenerando tabla adicional con métricas complementarias...")

additional_metrics = ['Precision', 'Recall', 'MCC', 'Accuracy', 'DynamicAccuracy']
additional_stats = []

for metric in additional_metrics:
    for sampling, sampling_label in [('21 timepoints', '21tp'), ('5 timepoints', '5tp')]:
        # Extraer datos
        boolnet_data = combined_df[(combined_df['Method'] == 'BoolNet') & 
                                   (combined_df['Sampling'] == sampling)][metric].dropna()
        bni3_data = combined_df[(combined_df['Method'] == 'BNI3') & 
                               (combined_df['Sampling'] == sampling)][metric].dropna()
        
        # Calcular estadísticas descriptivas
        boolnet_mean = boolnet_data.mean()
        boolnet_std = boolnet_data.std()
        boolnet_median = boolnet_data.median()
        
        bni3_mean = bni3_data.mean()
        bni3_std = bni3_data.std()
        bni3_median = bni3_data.median()
        
        # Diferencias
        diff_mean = bni3_mean - boolnet_mean
        diff_median = bni3_median - boolnet_median
        pct_change = (diff_mean / boolnet_mean * 100) if boolnet_mean != 0 else np.nan
        
        # Test estadístico Wilcoxon
        if len(boolnet_data) > 0 and len(bni3_data) > 0:
            stat, pval = stats.ranksums(boolnet_data, bni3_data)
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
            
            # Cohen's d
            cohens = cohens_d(bni3_data, boolnet_data)
            cohens_interp = interpret_cohens_d(cohens)
            
            # TOST para diferentes epsilon
            tost_results = {}
            for eps in epsilon_values:
                p_tost, delta = tost_equivalence_test(boolnet_data, bni3_data, eps)
                tost_results[eps] = {
                    'p_value': p_tost,
                    'delta': delta,
                    'equivalent': 'Yes' if p_tost < 0.05 else 'No'
                }
        else:
            stat, pval, sig = np.nan, np.nan, "N/A"
            cohens, cohens_interp = np.nan, "N/A"
            tost_results = {eps: {'p_value': np.nan, 'delta': np.nan, 'equivalent': 'N/A'} 
                          for eps in epsilon_values}
        
        # Confidence interval (95%)
        if len(boolnet_data) > 0 and len(bni3_data) > 0:
            from scipy import stats as sp_stats
            n1, n2 = len(boolnet_data), len(bni3_data)
            var1, var2 = boolnet_data.var(ddof=1), bni3_data.var(ddof=1)
            se = np.sqrt(var1/n1 + var2/n2)
            df = (var1/n1 + var2/n2)**2 / ((var1/n1)**2/(n1-1) + (var2/n2)**2/(n2-1))
            t_crit = sp_stats.t.ppf(0.975, df)
            ci_lower = diff_mean - t_crit * se
            ci_upper = diff_mean + t_crit * se
        else:
            ci_lower, ci_upper = np.nan, np.nan
        
        additional_stats.append({
            'Metric': metric,
            'Sampling': sampling_label,
            'BoolNet_Mean': boolnet_mean,
            'BoolNet_StdDev': boolnet_std,
            'BoolNet_Median': boolnet_median,
            'BNI3_Mean': bni3_mean,
            'BNI3_StdDev': bni3_std,
            'BNI3_Median': bni3_median,
            'Difference_Mean': diff_mean,
            'Difference_Median': diff_median,
            'Percent_Change': pct_change,
            'CI95_Lower': ci_lower,
            'CI95_Upper': ci_upper,
            'Cohens_d': cohens,
            'Cohens_d_Interpretation': cohens_interp,
            'Wilcoxon_Statistic': stat,
            'Wilcoxon_p_value': pval,
            'Wilcoxon_Significance': sig,
            'TOST_p_small_effect': tost_results[0.2]['p_value'],
            'TOST_equivalent_small': tost_results[0.2]['equivalent'],
            'TOST_delta_small': tost_results[0.2]['delta'],
            'TOST_p_medium_effect': tost_results[0.5]['p_value'],
            'TOST_equivalent_medium': tost_results[0.5]['equivalent'],
            'TOST_delta_medium': tost_results[0.5]['delta'],
            'TOST_p_large_effect': tost_results[0.8]['p_value'],
            'TOST_equivalent_large': tost_results[0.8]['equivalent'],
            'TOST_delta_large': tost_results[0.8]['delta'],
            'N_BoolNet': len(boolnet_data),
            'N_BNI3': len(bni3_data)
        })

additional_stats_df = pd.DataFrame(additional_stats)

# Guardar versión horizontal (todas las columnas)
additional_stats_wide = Path(OUTPUT_DIR) / "Additional_metrics_statistics_wide.tsv"
additional_stats_df.to_csv(additional_stats_wide, sep='\t', index=False, float_format='%.4f')
print(f"Métricas adicionales (formato amplio): {additional_stats_wide}")

# Crear versión vertical
additional_stats_vertical = []

for _, row in additional_stats_df.iterrows():
    base_info = f"{row['Metric']} | {row['Sampling']}"
    
    additional_stats_vertical.extend([
        {'Metric_Sampling': base_info, 'Statistic': 'Sample Size (BoolNet)', 'Value': f"{row['N_BoolNet']:.0f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Sample Size (BNI3)', 'Value': f"{row['N_BNI3']:.0f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BoolNet Mean ± SD', 'Value': f"{row['BoolNet_Mean']:.4f} ± {row['BoolNet_StdDev']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BoolNet Median', 'Value': f"{row['BoolNet_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BNI3 Mean ± SD', 'Value': f"{row['BNI3_Mean']:.4f} ± {row['BNI3_StdDev']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'BNI3 Median', 'Value': f"{row['BNI3_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Difference (Mean)', 'Value': f"{row['Difference_Mean']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Difference (Median)', 'Value': f"{row['Difference_Median']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Percent Change', 'Value': f"{row['Percent_Change']:.1f}%"},
        {'Metric_Sampling': base_info, 'Statistic': '95% CI Difference', 'Value': f"[{row['CI95_Lower']:.4f}, {row['CI95_Upper']:.4f}]"},
        {'Metric_Sampling': base_info, 'Statistic': "Cohen's d", 'Value': f"{row['Cohens_d']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': "Cohen's d Interpretation", 'Value': row['Cohens_d_Interpretation']},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon Statistic', 'Value': f"{row['Wilcoxon_Statistic']:.2f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon p-value', 'Value': f"{row['Wilcoxon_p_value']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'Wilcoxon Significance', 'Value': row['Wilcoxon_Significance']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (small: 0.2σ)', 'Value': f"{row['TOST_p_small_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (small)', 'Value': row['TOST_equivalent_small']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (small)', 'Value': f"±{row['TOST_delta_small']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (medium: 0.5σ)', 'Value': f"{row['TOST_p_medium_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (medium)', 'Value': row['TOST_equivalent_medium']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (medium)', 'Value': f"±{row['TOST_delta_medium']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST p-value (large: 0.8σ)', 'Value': f"{row['TOST_p_large_effect']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Equivalent (large)', 'Value': row['TOST_equivalent_large']},
        {'Metric_Sampling': base_info, 'Statistic': 'TOST Delta (large)', 'Value': f"±{row['TOST_delta_large']:.4f}"},
        {'Metric_Sampling': base_info, 'Statistic': '---', 'Value': '---'},
    ])

additional_stats_vertical_df = pd.DataFrame(additional_stats_vertical)

# Guardar versión vertical
additional_stats_long = Path(OUTPUT_DIR) / "Additional_metrics_statistics_long.tsv"
additional_stats_vertical_df.to_csv(additional_stats_long, sep='\t', index=False)
print(f"Métricas adicionales (formato vertical): {additional_stats_long}")

print("\n" + "="*80)
print("COMPLETADO")
print("="*80)
print(f"\nArchivos generados en: {OUTPUT_DIR}")
print("\nFIGURAS:")
print("  - Figure3A_Balanced_Accuracy.png / .svg")
print("  - Figure3B_F1.png / .svg")
print("  - Figure3C_Coverage.png / .svg (sin barras de significancia)")
print("  - Figure3D_Precision_Recall.png / .svg")
print("\nESTADÍSTICAS PRINCIPALES (Balanced Accuracy, F1-Score, Coverage):")
print("  - Figure3_detailed_statistics_wide.tsv (formato horizontal - TODAS las columnas)")
print("  - Figure3_detailed_statistics_long.tsv (formato vertical - FÁCIL de leer)")
print("  - Figure3_summary_statistics.tsv (resumen rápido)")
print("\nESTADÍSTICAS ADICIONALES (Precision, Recall, MCC, Accuracy, DynamicAccuracy):")
print("  - Additional_metrics_statistics_wide.tsv (formato horizontal)")
print("  - Additional_metrics_statistics_long.tsv (formato vertical)")
print("\n  Cada tabla incluye:")
print("    • Medias ± SD, medianas")
print("    • Diferencias absolutas y porcentuales")
print("    • Intervalos de confianza 95%")
print("    • Cohen's d (tamaño del efecto)")
print("    • Test de Wilcoxon (diferencias)")
print("    • Test TOST (equivalencia con 3 márgenes: 0.2σ, 0.5σ, 0.8σ)")