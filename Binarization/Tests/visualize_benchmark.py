#!/usr/bin/env python3
"""
Benchmark Results Visualization Tool
Generates publication-ready figures from binarization benchmark results.
"""

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication-quality defaults
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9


def load_benchmark_data(input_paths):
    """
    Load and validate benchmark results from one or more paths
    
    Args:
        input_paths: Can be:
            - Single CSV file path
            - Single directory path (looks for benchmark_results.csv)
            - List of CSV files
            - List of directory paths
    
    Returns:
        DataFrame with merged results or None if error
    """
    
    # Convert single path to list
    if isinstance(input_paths, str):
        input_paths = [input_paths]
    
    all_dataframes = []
    
    for input_path in input_paths:
        path = Path(input_path)
        
        # Determine if it's a file or directory
        if path.is_file():
            csv_file = path
        elif path.is_dir():
            # Look for benchmark_results.csv in directory
            csv_file = path / 'benchmark_results.csv'
            if not csv_file.exists():
                print(f"WARNING: No benchmark_results.csv found in {path}")
                continue
        else:
            print(f"WARNING: Path not found: {input_path}")
            continue
        
        # Load CSV
        try:
            df = pd.read_csv(csv_file)
            
            # Validate required columns
            required_cols = ['algorithm', 'n_genes', 'n_timepoints', 'n_processors', 
                            'processing_time', 'success']
            
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f"WARNING: Missing columns in {csv_file}: {missing_cols}")
                continue
            
            # Filter only successful runs
            df = df[df['success'] == True].copy()
            
            if len(df) == 0:
                print(f"WARNING: No successful runs in {csv_file}")
                continue
            
            all_dataframes.append(df)
            print(f"✓ Loaded: {csv_file.name} ({len(df)} successful runs)")
            
        except Exception as e:
            print(f"ERROR loading {csv_file}: {e}")
            continue
    
    if not all_dataframes:
        print("ERROR: No valid benchmark data loaded")
        return None
    
    # Concatenate all dataframes
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    
    # Remove duplicates (same algorithm, genes, timepoints, processors)
    original_len = len(merged_df)
    merged_df = merged_df.drop_duplicates(
        subset=['algorithm', 'n_genes', 'n_timepoints', 'n_processors'],
        keep='last'
    )
    
    if len(merged_df) < original_len:
        print(f"ℹ Removed {original_len - len(merged_df)} duplicate entries")
    
    print(f"\n📊 MERGED DATA SUMMARY:")
    print(f"Total successful runs: {len(merged_df)}")
    print(f"Algorithms: {sorted(merged_df['algorithm'].unique())}")
    print(f"Genes: {sorted(merged_df['n_genes'].unique())}")
    print(f"Timepoints: {sorted(merged_df['n_timepoints'].unique())}")
    print(f"Processors: {sorted(merged_df['n_processors'].unique())}")
    
    return merged_df


def option1_panels_by_processors(df, output_dir, log_scale=True):
    """
    Option 1: Separate panels for each processor count
    Shows how parallelization affects performance
    
    Args:
        log_scale: If True, use logarithmic Y-axis scale (default: True, recommended for wide range)
    """
    processors = sorted(df['n_processors'].unique())
    n_procs = len(processors)
    
    fig, axes = plt.subplots(n_procs, 2, figsize=(14, 4*n_procs))
    
    # If only one processor count, make axes 2D
    if n_procs == 1:
        axes = axes.reshape(1, -1)
    
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(df['n_genes'].unique())))
    
    scale_name = "Log Scale" if log_scale else "Linear Scale"
    
    for proc_idx, n_proc in enumerate(processors):
        data_proc = df[df['n_processors'] == n_proc]
        
        # SSD panel
        ax_ssd = axes[proc_idx, 0]
        ssd_data = data_proc[data_proc['algorithm'] == 'SSD']
        
        for gene_idx, n_genes in enumerate(sorted(ssd_data['n_genes'].unique())):
            gene_data = ssd_data[ssd_data['n_genes'] == n_genes].sort_values('n_timepoints')
            ax_ssd.plot(gene_data['n_timepoints'], gene_data['processing_time'],
                       marker='o', linewidth=2, markersize=6,
                       color=colors[gene_idx], label=f'{n_genes} genes')
        
        ax_ssd.set_xlabel('Timepoints', fontsize=11)
        ax_ssd.set_ylabel('Time (s)', fontsize=11)
        ax_ssd.set_title(f'SSD - {n_proc} Processor(s) ({scale_name})', fontsize=12)
        ax_ssd.legend(loc='best', fontsize=9)
        ax_ssd.grid(True, alpha=0.3)
        if log_scale:
            ax_ssd.set_yscale('log')
        
        # WCSS panel
        ax_wcss = axes[proc_idx, 1]
        wcss_data = data_proc[data_proc['algorithm'] == 'WCSS']
        
        for gene_idx, n_genes in enumerate(sorted(wcss_data['n_genes'].unique())):
            gene_data = wcss_data[wcss_data['n_genes'] == n_genes].sort_values('n_timepoints')
            ax_wcss.plot(gene_data['n_timepoints'], gene_data['processing_time'],
                        marker='s', linewidth=2, markersize=6,
                        color=colors[gene_idx], label=f'{n_genes} genes')
        
        ax_wcss.set_xlabel('Timepoints', fontsize=11)
        ax_wcss.set_ylabel('Time (s)', fontsize=11)
        ax_wcss.set_title(f'WCSS - {n_proc} Processor(s) ({scale_name})', fontsize=12)
        ax_wcss.legend(loc='best', fontsize=9)
        ax_wcss.grid(True, alpha=0.3)
        if log_scale:
            ax_wcss.set_yscale('log')
    
    plt.tight_layout()
    
    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    suffix = "_log" if log_scale else "_linear"
    fig.savefig(output_path / f'option1_panels_by_processors{suffix}.png', bbox_inches='tight')
    fig.savefig(output_path / f'option1_panels_by_processors{suffix}.svg', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Option 1 saved: panels_by_processors ({scale_name})")


def option2_heatmaps(df, output_dir):
    """
    Option 2: Heatmaps showing genes × timepoints for each algorithm and processor count
    """
    processors = sorted(df['n_processors'].unique())
    algorithms = sorted(df['algorithm'].unique())
    
    n_rows = len(algorithms)
    n_cols = len(processors)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4.5*n_rows))
    
    # Make axes 2D if needed
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    for alg_idx, algorithm in enumerate(algorithms):
        for proc_idx, n_proc in enumerate(processors):
            ax = axes[alg_idx, proc_idx]
            
            # Filter data
            data = df[(df['algorithm'] == algorithm) & (df['n_processors'] == n_proc)]
            
            if len(data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax.set_title(f'{algorithm} - {n_proc} proc(s)')
                continue
            
            # Pivot for heatmap
            pivot_table = data.pivot_table(
                values='processing_time',
                index='n_genes',
                columns='n_timepoints',
                aggfunc='mean'
            )
            
            # Smart formatting: fewer decimals for larger numbers
            def smart_format(x):
                if pd.isna(x):
                    return ''
                elif x < 1:
                    return f'{x:.2f}'
                elif x < 10:
                    return f'{x:.1f}'
                else:
                    return f'{x:.0f}'
            
            # Create annotation matrix with smart formatting
            annot_matrix = pivot_table.map(smart_format)
            
            # Create heatmap
            sns.heatmap(pivot_table, annot=annot_matrix, fmt='', cmap='YlOrRd',
                       ax=ax, cbar_kws={'label': 'Time (s)'}, linewidths=0.5,
                       annot_kws={'fontsize': 8})
            
            ax.set_xlabel('Timepoints')
            ax.set_ylabel('Genes')
            ax.set_title(f'{algorithm} - {n_proc} Processor(s)')
    
    plt.tight_layout()
    
    # Save
    output_path = Path(output_dir)
    fig.savefig(output_path / 'option2_heatmaps.png', bbox_inches='tight')
    fig.savefig(output_path / 'option2_heatmaps.svg', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Option 2 saved: heatmaps")


def option3_speedup_analysis(df, output_dir):
    """
    Option 3: Speedup and efficiency analysis for parallelization
    Shows only representative cases for clarity
    """
    # Calculate speedup (time with 1 proc / time with N procs)
    processors = sorted(df['n_processors'].unique())
    
    if len(processors) < 2:
        print("⚠ Option 3 requires multiple processor counts. Skipping.")
        return
    
    # Check for very small processing times that might cause issues
    min_time = df['processing_time'].min()
    if min_time < 0.01:
        print(f"ℹ Note: Some processing times are very small (< 0.01s). Speedup calculations may be less accurate for these cases.")
    
    # Get baseline (1 processor or minimum)
    baseline_proc = min(processors)
    
    # Select representative cases to show
    # Strategy: Show small, medium, and large datasets for each algorithm
    genes_available = sorted(df['n_genes'].unique())
    timepoints_available = sorted(df['n_timepoints'].unique())
    
    # Define representative combinations (small, medium, large)
    representative_cases = []
    
    # Small datasets
    if len(genes_available) >= 1 and len(timepoints_available) >= 1:
        representative_cases.append((genes_available[0], timepoints_available[0], 'Small'))
    
    # Medium datasets (middle values)
    if len(genes_available) >= 2 and len(timepoints_available) >= 2:
        mid_g = len(genes_available) // 2
        mid_t = len(timepoints_available) // 2
        representative_cases.append((genes_available[mid_g], timepoints_available[mid_t], 'Medium'))
    
    # Large datasets
    if len(genes_available) >= 1 and len(timepoints_available) >= 1:
        representative_cases.append((genes_available[-1], timepoints_available[-1], 'Large'))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax1, ax2 = axes
    
    # Color scheme: different colors for different sizes, different line styles for algorithms
    colors = {'Small': '#1f77b4', 'Medium': '#ff7f0e', 'Large': '#2ca02c'}
    markers = {'SSD': 'o', 'WCSS': 's'}
    linestyles = {'SSD': '-', 'WCSS': '--'}
    
    # For each algorithm
    for algorithm in sorted(df['algorithm'].unique()):
        alg_data = df[df['algorithm'] == algorithm]
        
        # For each representative case
        for n_genes, n_timepoints, size_label in representative_cases:
            
            subset = alg_data[(alg_data['n_genes'] == n_genes) & 
                             (alg_data['n_timepoints'] == n_timepoints)]
            
            if len(subset) < 2:
                continue
            
            # Get baseline time
            baseline = subset[subset['n_processors'] == baseline_proc]['processing_time'].values
            if len(baseline) == 0 or baseline[0] <= 0:
                continue
            baseline_time = baseline[0]
            
            # Calculate speedup for each processor count
            speedups = []
            efficiencies = []
            proc_counts = []
            
            for n_proc in processors:
                proc_time = subset[subset['n_processors'] == n_proc]['processing_time'].values
                if len(proc_time) > 0 and proc_time[0] > 0:
                    speedup = baseline_time / proc_time[0]
                    efficiency = (speedup / n_proc) * 100
                    
                    speedups.append(speedup)
                    efficiencies.append(efficiency)
                    proc_counts.append(n_proc)
            
            if not speedups:
                continue
            
            # Plot speedup
            label = f'{algorithm} {size_label} ({n_genes}g, {n_timepoints}t)'
            ax1.plot(proc_counts, speedups, 
                    marker=markers[algorithm], 
                    linestyle=linestyles[algorithm],
                    linewidth=2.5,
                    markersize=8,
                    color=colors[size_label], 
                    label=label,
                    alpha=0.8)
            
            # Plot efficiency
            ax2.plot(proc_counts, efficiencies, 
                    marker=markers[algorithm], 
                    linestyle=linestyles[algorithm],
                    linewidth=2.5,
                    markersize=8,
                    color=colors[size_label], 
                    label=label,
                    alpha=0.8)
    
    # Ideal speedup line
    ax1.plot(processors, processors, 'k--', linewidth=2, label='Ideal (linear)', alpha=0.7)
    
    ax1.set_xlabel('Number of Processors', fontsize=12)
    ax1.set_ylabel('Speedup', fontsize=12)
    ax1.set_title('Speedup vs Processors\n(Representative Cases)', fontsize=13, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)
    
    ax2.axhline(y=100, color='k', linestyle='--', linewidth=2, label='Ideal (100%)', alpha=0.7)
    ax2.set_xlabel('Number of Processors', fontsize=12)
    ax2.set_ylabel('Efficiency (%)', fontsize=12)
    ax2.set_title('Parallel Efficiency\n(Representative Cases)', fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 110)
    
    # Save
    output_path = Path(output_dir)
    fig.savefig(output_path / 'option3_speedup_analysis.png', bbox_inches='tight')
    fig.savefig(output_path / 'option3_speedup_analysis.svg', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Option 3 saved: speedup_analysis")


def option4_combined_figure(df, output_dir):
    """
    Option 4: Simplified figure focusing on scaling behavior
    A) Time per gene vs Timepoints (SSD)
    B) Time per gene vs Timepoints (WCSS)
    C) Scaling with gene count (log-log)
    D) Scaling with timepoint count (log-log)
    """
    processors = sorted(df['n_processors'].unique())
    gene_counts = sorted(df['n_genes'].unique())
    timepoint_counts = sorted(df['n_timepoints'].unique())
    
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)
    
    # Colors for different processor counts
    proc_colors = {processors[i]: plt.cm.viridis(i / max(1, len(processors)-1)) 
                   for i in range(len(processors))}
    
    # Panel A: SSD - Time per gene (normalized performance)
    ax1 = fig.add_subplot(gs[0, 0])
    
    for n_proc in processors:
        subset = df[(df['algorithm'] == 'SSD') & (df['n_processors'] == n_proc)].copy()
        if len(subset) == 0:
            continue
        
        subset['time_per_gene'] = subset['processing_time'] / subset['n_genes']
        grouped = subset.groupby('n_timepoints')['time_per_gene'].mean().reset_index()
        
        ax1.plot(grouped['n_timepoints'], grouped['time_per_gene'],
                marker='o', linewidth=2.5, markersize=7,
                color=proc_colors[n_proc], label=f'{n_proc}p', alpha=0.8)
    
    ax1.set_xlabel('Timepoints', fontsize=12)
    ax1.set_ylabel('Time per Gene (s/gene)', fontsize=12)
    ax1.set_title('A) SSD: Normalized Performance', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # Panel B: WCSS - Time per gene (normalized performance)
    ax2 = fig.add_subplot(gs[0, 1])
    
    for n_proc in processors:
        subset = df[(df['algorithm'] == 'WCSS') & (df['n_processors'] == n_proc)].copy()
        if len(subset) == 0:
            continue
        
        subset['time_per_gene'] = subset['processing_time'] / subset['n_genes']
        grouped = subset.groupby('n_timepoints')['time_per_gene'].mean().reset_index()
        
        ax2.plot(grouped['n_timepoints'], grouped['time_per_gene'],
                marker='s', linewidth=2.5, markersize=7,
                color=proc_colors[n_proc], label=f'{n_proc}p', alpha=0.8)
    
    ax2.set_xlabel('Timepoints', fontsize=12)
    ax2.set_ylabel('Time per Gene (s/gene)', fontsize=12)
    ax2.set_title('B) WCSS: Normalized Performance', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    # Panel C: Scaling with gene count (both algorithms, fixed timepoint)
    ax3 = fig.add_subplot(gs[1, 0])
    
    # Use middle or largest timepoint for comparison
    comparison_timepoint = timepoint_counts[-1] if len(timepoint_counts) > 2 else timepoint_counts[len(timepoint_counts)//2]
    
    for algorithm in ['SSD', 'WCSS']:
        # Use only one processor count for clarity (the minimum)
        n_proc = min(processors)
        subset = df[(df['algorithm'] == algorithm) & 
                   (df['n_processors'] == n_proc) &
                   (df['n_timepoints'] == comparison_timepoint)].sort_values('n_genes')
        
        if len(subset) == 0:
            continue
        
        marker = 'o' if algorithm == 'SSD' else 's'
        linestyle = '-' if algorithm == 'SSD' else '--'
        color = '#1f77b4' if algorithm == 'SSD' else '#2ca02c'
        
        ax3.plot(subset['n_genes'], subset['processing_time'],
                marker=marker, linestyle=linestyle, linewidth=2.5, markersize=7,
                color=color, label=algorithm, alpha=0.8)
    
    ax3.set_xlabel('Number of Genes', fontsize=12)
    ax3.set_ylabel('Processing Time (s)', fontsize=12)
    ax3.set_title(f'C) Scaling with Gene Count\n({comparison_timepoint} timepoints, {min(processors)}p)', fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    
    # Panel D: Scaling with timepoint count (both algorithms, fixed genes)
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Use middle or largest gene count for comparison
    comparison_genes = gene_counts[-1] if len(gene_counts) > 2 else gene_counts[len(gene_counts)//2]
    
    for algorithm in ['SSD', 'WCSS']:
        # Use only one processor count for clarity
        n_proc = min(processors)
        subset = df[(df['algorithm'] == algorithm) & 
                   (df['n_processors'] == n_proc) &
                   (df['n_genes'] == comparison_genes)].sort_values('n_timepoints')
        
        if len(subset) == 0:
            continue
        
        marker = 'o' if algorithm == 'SSD' else 's'
        linestyle = '-' if algorithm == 'SSD' else '--'
        color = '#1f77b4' if algorithm == 'SSD' else '#2ca02c'
        
        ax4.plot(subset['n_timepoints'], subset['processing_time'],
                marker=marker, linestyle=linestyle, linewidth=2.5, markersize=7,
                color=color, label=algorithm, alpha=0.8)
    
    ax4.set_xlabel('Number of Timepoints', fontsize=12)
    ax4.set_ylabel('Processing Time (s)', fontsize=12)
    ax4.set_title(f'D) Scaling with Timepoint Count\n({comparison_genes} genes, {min(processors)}p)', fontsize=13, fontweight='bold')
    ax4.legend(loc='best', fontsize=11)
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')
    
    # Save
    output_path = Path(output_dir)
    fig.savefig(output_path / 'option4_combined_figure.png', bbox_inches='tight', dpi=300)
    fig.savefig(output_path / 'option4_combined_figure.svg', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Option 4 saved: combined_figure")



def main():
    parser = argparse.ArgumentParser(
        description="Visualize binarization benchmark results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single CSV file
  python3 visualize_benchmark.py -i benchmark_results.csv -o figures
  
  # Single directory (looks for benchmark_results.csv inside)
  python3 visualize_benchmark.py -i results_p1/ -o figures
  
  # Multiple directories (merges results automatically)
  python3 visualize_benchmark.py -i results_p1/ results_p8/ -o figures
  
  # Multiple CSV files
  python3 visualize_benchmark.py -i results_p1/benchmark_results.csv results_p8/benchmark_results.csv -o figures
  
  # Specific options only
  python3 visualize_benchmark.py -i results_p1/ results_p8/ -o plots --options 2,4
        """
    )
    
    parser.add_argument('-i', '--input',
                        nargs='+',
                        required=True,
                        help='Input CSV file(s) or directory(ies) with benchmark results (space-separated for multiple)')
    
    parser.add_argument('-o', '--output',
                        default='./figures',
                        help='Output directory for figures (default: ./figures)')
    
    parser.add_argument('--options',
                        default='1,2,3,4',
                        help='Comma-separated list of figure options to generate (default: 1,2,3,4)')
    
    parser.add_argument('--version',
                        action='version',
                        version='Benchmark Visualizer v1.0')
    
    args = parser.parse_args()
    
    # Parse options
    try:
        options = [int(x.strip()) for x in args.options.split(',')]
    except ValueError:
        print("ERROR: Invalid options format. Use comma-separated integers (e.g., 1,2,3)")
        return 1
    
    print("="*60)
    print("BENCHMARK VISUALIZATION")
    print("="*60)
    print(f"Input paths: {len(args.input)}")
    for inp in args.input:
        print(f"  - {inp}")
    print(f"Output dir:    {args.output}")
    print(f"Options:       {options}")
    print("="*60 + "\n")
    
    # Load data (handles merging automatically)
    df = load_benchmark_data(args.input)
    if df is None:
        return 1
    
    print("\nGenerating figures...\n")
    
    # Generate requested figures
    if 1 in options:
        try:
            option1_panels_by_processors(df, args.output, log_scale=True)
            option1_panels_by_processors(df, args.output, log_scale=False)
        except Exception as e:
            print(f"✗ Option 1 failed: {e}")
    
    if 2 in options:
        try:
            option2_heatmaps(df, args.output)
        except Exception as e:
            print(f"✗ Option 2 failed: {e}")
    
    if 3 in options:
        try:
            option3_speedup_analysis(df, args.output)
        except Exception as e:
            print(f"✗ Option 3 failed: {e}")
    
    if 4 in options:
        try:
            option4_combined_figure(df, args.output)
        except Exception as e:
            print(f"✗ Option 4 failed: {e}")
    
    print("\n" + "="*60)
    print(f"Figures saved to: {Path(args.output).absolute()}")
    print("="*60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())