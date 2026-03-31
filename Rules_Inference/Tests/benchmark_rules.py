#!/usr/bin/env python3
"""
GEP Rule Inference Benchmark Tool
Runs GEP-based Boolean network inference on test datasets and generates performance metrics.
"""

import argparse
import sys
import os
import subprocess
import time
import re
from pathlib import Path
import pandas as pd
import numpy as np


def parse_filename(filename):
    """
    Parse binarized dataset filename to extract genes and timepoints
    Expected format: bin_ssd_g{n_genes}_t{n_timepoints}_p{n_proc}_test_data_g{n_genes}_t{n_timepoints}.tsv
    """
    pattern = r'bin_ssd_g(\d+)_t(\d+)_p\d+_test_data_g\d+_t\d+\.tsv'
    match = re.match(pattern, filename)
    
    if match:
        n_genes = int(match.group(1))
        n_timepoints = int(match.group(2))
        return n_genes, n_timepoints
    
    return None, None


def find_raw_data_file(binary_file, raw_data_dir):
    """
    Find the corresponding raw data file for a binary file
    
    Args:
        binary_file: Path to binary file
        raw_data_dir: Directory containing raw data files
    
    Returns:
        Path to raw data file or None if not found
    """
    # Extract g and t from binary filename
    n_genes, n_timepoints = parse_filename(binary_file.name)
    
    if n_genes is None or n_timepoints is None:
        return None
    
    # Look for raw data file
    raw_filename = f"test_data_g{n_genes}_t{n_timepoints}.tsv"
    raw_path = Path(raw_data_dir) / raw_filename
    
    if raw_path.exists():
        return raw_path
    
    return None


def filter_datasets(binary_files, gene_counts=None, timepoint_counts=None):
    """
    Filter dataset files based on gene and timepoint counts
    
    Args:
        binary_files: List of binary file paths
        gene_counts: List of gene counts to include (None = all)
        timepoint_counts: List of timepoint counts to include (None = all)
    
    Returns:
        Filtered list of binary file paths
    """
    filtered = []
    
    for binary_file in binary_files:
        n_genes, n_timepoints = parse_filename(binary_file.name)
        
        if n_genes is None or n_timepoints is None:
            continue
        
        # Filter by gene count if specified
        if gene_counts is not None and n_genes not in gene_counts:
            continue
        
        # Filter by timepoint count if specified
        if timepoint_counts is not None and n_timepoints not in timepoint_counts:
            continue
        
        filtered.append(binary_file)
    
    return filtered


def run_gep_inference(raw_file, binary_file, output_dir, n_processors=1, 
                     population=50, generations=50, n_repetitions=3,
                     target_genes=None, gep_script=None, verbose=False):
    """
    Run GEP inference and measure execution time
    
    Args:
        raw_file: Path to raw data file
        binary_file: Path to binary data file
        output_dir: Output directory for results
        n_processors: Number of processors to use
        population: Population size for GEP
        generations: Number of generations
        n_repetitions: Number of repetitions per gene
        target_genes: Specific genes to target (None = all genes)
        gep_script: Custom path to GEP script
        verbose: Enable verbose output
    
    Returns:
        dict: Results including time, success status, and metrics
    """
    
    # Determine which script to use
    if gep_script:
        script = gep_script
    else:
        script = 'BNI3_GEP.py'
    
    # Validate script exists
    if not os.path.exists(script):
        raise FileNotFoundError(f"Script not found: {script}")
    
    # Build command
    cmd = [
        'python3',
        script,
        '-i', str(raw_file),
        '-i_binary', str(binary_file),
        '-o', str(output_dir),
        '-p', str(n_processors),
        '-pop', str(population),
        '-gen', str(generations),
        '-n_rep', str(n_repetitions),
        '--no_analysis'  # Skip automatic analysis for benchmarking
    ]
    
    if target_genes:
        cmd.extend(['-targets', target_genes])
    
    if verbose:
        cmd.append('-v')
    
    print(f"  Running GEP (p={n_processors})...", end=' ', flush=True)
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=3600  # 1 hour timeout
        )
        
        elapsed_time = time.time() - start_time
        
        # Parse metrics from output
        total_experiments = None
        successful_experiments = None
        avg_time_per_exp = None
        
        for line in result.stdout.split('\n'):
            if 'Total experiments:' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    total_experiments = int(match.group(1))
            elif 'Successful:' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    successful_experiments = int(match.group(1))
            elif 'Average time/exp:' in line:
                match = re.search(r'(\d+\.?\d*)\s+seconds', line)
                if match:
                    avg_time_per_exp = float(match.group(1))
        
        print(f"✓ ({elapsed_time:.2f}s)")
        
        return {
            'success': True,
            'total_time': elapsed_time,
            'total_experiments': total_experiments,
            'successful_experiments': successful_experiments,
            'avg_time_per_exp': avg_time_per_exp,
            'error': None
        }
        
    except subprocess.TimeoutExpired:
        elapsed_time = time.time() - start_time
        print(f"✗ (Timeout after {elapsed_time:.2f}s)")
        
        return {
            'success': False,
            'total_time': elapsed_time,
            'total_experiments': None,
            'successful_experiments': None,
            'avg_time_per_exp': None,
            'error': 'Timeout after 1 hour'
        }
        
    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print(f"✗ (Failed after {elapsed_time:.2f}s)")
        
        return {
            'success': False,
            'total_time': elapsed_time,
            'total_experiments': None,
            'successful_experiments': None,
            'avg_time_per_exp': None,
            'error': str(e)
        }


def benchmark_datasets(binary_dir, raw_dir, output_base_dir, processor_counts,
                      gene_counts, timepoint_counts, population, generations, 
                      n_repetitions, target_genes, gep_script, verbose):
    """
    Run GEP inference on filtered test datasets with different processor counts
    
    Args:
        binary_dir: Directory containing binary datasets
        raw_dir: Directory containing raw datasets
        output_base_dir: Base directory for output files
        processor_counts: List of processor counts to test
        gene_counts: List of gene counts to test (None = all)
        timepoint_counts: List of timepoint counts to test (None = all)
        population: Population size for GEP
        generations: Number of generations
        n_repetitions: Number of repetitions per gene
        target_genes: Specific genes to target (None = all genes)
        gep_script: Custom path to GEP script
        verbose: Enable verbose output
    """
    
    binary_path = Path(binary_dir)
    output_base = Path(output_base_dir)
    
    # Find all SSD binary datasets
    all_binary_files = sorted(binary_path.glob('bin_ssd_g*_t*_p*_test_data_g*_t*.tsv'))
    
    if not all_binary_files:
        print(f"ERROR: No SSD binary datasets found in {binary_dir}")
        print("Expected format: bin_ssd_g{n_genes}_t{n_timepoints}_p{n_proc}_test_data_g{n_genes}_t{n_timepoints}.tsv")
        return None
    
    # Filter datasets based on gene and timepoint counts
    binary_files = filter_datasets(all_binary_files, gene_counts, timepoint_counts)
    
    if not binary_files:
        print(f"ERROR: No datasets match the specified criteria")
        if gene_counts:
            print(f"  Gene counts: {gene_counts}")
        if timepoint_counts:
            print(f"  Timepoint counts: {timepoint_counts}")
        return None
    
    print(f"\nFound {len(all_binary_files)} total SSD binary datasets")
    print(f"Selected {len(binary_files)} datasets matching criteria")
    print("="*80)
    
    results = []
    
    for binary_file in binary_files:
        # Parse filename
        n_genes, n_timepoints = parse_filename(binary_file.name)
        
        if n_genes is None or n_timepoints is None:
            print(f"WARNING: Skipping file with unexpected name: {binary_file.name}")
            continue
        
        # Find corresponding raw data file
        raw_file = find_raw_data_file(binary_file, raw_dir)
        
        if raw_file is None:
            print(f"WARNING: Raw data file not found for {binary_file.name}")
            continue
        
        print(f"\n{binary_file.name}")
        print(f"  Dataset: {n_genes} genes × {n_timepoints} timepoints")
        
        # Run with different processor counts
        for n_proc in processor_counts:
            # Create output directory for this run
            output_dir = output_base / f"gep_g{n_genes}_t{n_timepoints}_p{n_proc}_pop{population}_gen{generations}_rep{n_repetitions}"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run GEP inference
            result = run_gep_inference(
                raw_file=raw_file,
                binary_file=binary_file,
                output_dir=output_dir,
                n_processors=n_proc,
                population=population,
                generations=generations,
                n_repetitions=n_repetitions,
                target_genes=target_genes,
                gep_script=gep_script,
                verbose=verbose
            )
            
            # Add metadata
            result.update({
                'binary_file': binary_file.name,
                'raw_file': raw_file.name,
                'n_genes': n_genes,
                'n_timepoints': n_timepoints,
                'n_processors': n_proc,
                'population': population,
                'generations': generations,
                'n_repetitions': n_repetitions
            })
            
            # Calculate per-gene metrics
            if result['success'] and result['avg_time_per_exp'] is not None:
                result['total_processing_time'] = result['avg_time_per_exp'] * result['total_experiments']
                result['time_per_gene'] = result['total_processing_time'] / n_genes if n_genes > 0 else None
                result['time_per_gene_per_rep'] = result['time_per_gene'] / n_repetitions if n_repetitions > 0 else None
            else:
                result['total_processing_time'] = None
                result['time_per_gene'] = None
                result['time_per_gene_per_rep'] = None
            
            results.append(result)
    
    return results


def generate_summary_report(results, output_dir):
    """
    Generate comprehensive summary report with performance metrics
    """
    
    if not results:
        print("\nNo results to summarize.")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Create output directory for reports
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save detailed results to CSV
    results_file = output_path / 'inference_benchmark_results.csv'
    df.to_csv(results_file, index=False, float_format='%.6f')
    
    print("\n" + "="*80)
    print("GEP INFERENCE BENCHMARK SUMMARY")
    print("="*80)
    
    # Overall statistics
    print("\n📊 OVERALL STATISTICS:")
    print("-" * 80)
    
    successful = df['success'].sum()
    total = len(df)
    
    print(f"Total runs:        {total}")
    print(f"Successful:        {successful}")
    print(f"Failed:            {total - successful}")
    print(f"Success rate:      {100*successful/total:.1f}%")
    
    if successful > 0:
        success_data = df[df['success']]
        print(f"\nTotal time:        {success_data['total_time'].sum():.2f} seconds")
        print(f"Avg time/run:      {success_data['total_time'].mean():.2f} seconds")
        print(f"Min time/run:      {success_data['total_time'].min():.2f} seconds")
        print(f"Max time/run:      {success_data['total_time'].max():.2f} seconds")
    
    # Performance by processor count
    print("\n⚡ PERFORMANCE BY PROCESSOR COUNT:")
    print("-" * 80)
    
    successful_data = df[df['success']]
    
    if len(successful_data) > 0:
        for n_proc in sorted(successful_data['n_processors'].unique()):
            proc_data = successful_data[successful_data['n_processors'] == n_proc]
            
            print(f"\n{n_proc} Processor(s):")
            print(f"  Runs:                        {len(proc_data)}")
            print(f"  Avg total time:              {proc_data['total_time'].mean():.2f} seconds")
            print(f"  Avg time per gene:           {proc_data['time_per_gene'].mean():.4f} seconds")
            print(f"  Avg time per gene per rep:   {proc_data['time_per_gene_per_rep'].mean():.4f} seconds")
            print(f"  Throughput (genes/sec):      {(1/proc_data['time_per_gene'].mean()):.4f}")
    
    # Detailed results table
    print("\n📋 DETAILED RESULTS:")
    print("-" * 80)
    print(f"{'Genes':>6} {'Times':>6} {'Procs':>6} {'Pop':>5} {'Gen':>5} {'Rep':>5} "
          f"{'Total(s)':>10} {'Time/Gene':>12} {'Status':<10}")
    print("-" * 80)
    
    for _, row in df.iterrows():
        status = "✓ Success" if row['success'] else "✗ Failed"
        time_str = f"{row['total_time']:.2f}" if row['success'] else "N/A"
        time_per_gene_str = f"{row['time_per_gene']:.6f}" if row['success'] and row['time_per_gene'] is not None else "N/A"
        
        print(f"{row['n_genes']:>6} {row['n_timepoints']:>6} {row['n_processors']:>6} "
              f"{row['population']:>5} {row['generations']:>5} {row['n_repetitions']:>5} "
              f"{time_str:>10} {time_per_gene_str:>12} {status:<10}")
    
    # Scaling analysis
    print("\n📈 SCALING ANALYSIS:")
    print("-" * 80)
    
    scaling_data = []
    
    # Group by dataset (genes × timepoints) and compare processor counts
    for n_genes in sorted(df['n_genes'].unique()):
        for n_timepoints in sorted(df['n_timepoints'].unique()):
            dataset_data = successful_data[
                (successful_data['n_genes'] == n_genes) & 
                (successful_data['n_timepoints'] == n_timepoints)
            ]
            
            if len(dataset_data) > 1:
                # Sort by processor count
                dataset_data = dataset_data.sort_values('n_processors')
                
                baseline_time = None
                baseline_procs = None
                
                for _, row in dataset_data.iterrows():
                    if baseline_time is None:
                        baseline_time = row['total_time']
                        baseline_procs = row['n_processors']
                        speedup = 1.0
                        efficiency = 1.0
                    else:
                        speedup = baseline_time / row['total_time']
                        efficiency = speedup / (row['n_processors'] / baseline_procs)
                    
                    scaling_data.append({
                        'n_genes': n_genes,
                        'n_timepoints': n_timepoints,
                        'n_processors': row['n_processors'],
                        'total_time': row['total_time'],
                        'speedup': speedup,
                        'efficiency': efficiency * 100
                    })
    
    if scaling_data:
        scaling_df = pd.DataFrame(scaling_data)
        
        print(f"\n{'Genes':>6} {'Times':>6} {'Procs':>6} {'Time(s)':>10} {'Speedup':>10} {'Efficiency(%)':>15}")
        print("-" * 80)
        
        for _, row in scaling_df.iterrows():
            print(f"{row['n_genes']:>6} {row['n_timepoints']:>6} {row['n_processors']:>6} "
                  f"{row['total_time']:>10.2f} {row['speedup']:>10.2f}x {row['efficiency']:>14.1f}%")
        
        # Save scaling analysis to CSV
        scaling_file = output_path / 'inference_scaling_analysis.csv'
        scaling_df.to_csv(scaling_file, index=False, float_format='%.6f')
        
        print(f"\n📁 Scaling analysis saved to: {scaling_file}")
    
    print("\n" + "="*80)
    print(f"📁 Detailed results saved to: {results_file}")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark GEP-based Boolean network inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic benchmark with 1,2,4 processors (all datasets)
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/
  
  # Custom processor counts
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -p 1,2,4,8
  
  # Test specific gene counts only
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -g 10,20,50
  
  # Test specific timepoint counts only
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -t 10,25,50
  
  # Test specific combinations
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -g 10,20 -t 10,25 -p 1,8,16
  
  # Quick test: small datasets with few processors
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -g 10 -t 10,25 -p 1,8
  
  # Small-scale test with reduced GEP parameters
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -pop 30 -gen 30 -n_rep 2
  
  # Test specific genes within datasets
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ -targets Gene1,Gene2,Gene3
  
  # Use custom GEP script
  python3 benchmark_inference.py -i_binary binarized/ -i_raw raw_data/ -o results/ --gep_script /path/to/BNI3_GEP.py
        """
    )
    
    parser.add_argument('-i_binary', '--input_binary',
                        required=True,
                        help='Input directory containing binarized datasets (SSD only)')
    
    parser.add_argument('-i_raw', '--input_raw',
                        required=True,
                        help='Input directory containing raw datasets')
    
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory for results and reports')
    
    parser.add_argument('-p', '--processors',
                        default='1,2,4',
                        help='Comma-separated list of processor counts to test (default: 1,2,4)')
    
    parser.add_argument('-g', '--genes',
                        default=None,
                        help='Comma-separated list of gene counts to test (default: all available)')
    
    parser.add_argument('-t', '--timepoints',
                        default=None,
                        help='Comma-separated list of timepoint counts to test (default: all available)')
    
    parser.add_argument('-pop', '--population',
                        type=int,
                        default=50,
                        help='Population size for GEP (default: 50)')
    
    parser.add_argument('-gen', '--generations',
                        type=int,
                        default=50,
                        help='Number of generations (default: 50)')
    
    parser.add_argument('-n_rep', '--n_repetitions',
                        type=int,
                        default=3,
                        help='Number of repetitions per gene (default: 3)')
    
    parser.add_argument('-targets', '--target_genes',
                        default=None,
                        help='Target specific genes (comma-separated). Default: all genes')
    
    parser.add_argument('--gep_script',
                        default=None,
                        help='Path to GEP script (default: BNI3_GEP.py in current directory)')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Enable verbose output')
    
    parser.add_argument('--version',
                        action='version',
                        version='GEP Inference Benchmark v2.0')
    
    args = parser.parse_args()
    
    # Validate directories
    if not os.path.exists(args.input_binary):
        print(f"ERROR: Binary data directory '{args.input_binary}' does not exist.")
        return 1
    
    if not os.path.exists(args.input_raw):
        print(f"ERROR: Raw data directory '{args.input_raw}' does not exist.")
        return 1
    
    # Parse processor counts
    try:
        processor_counts = [int(p.strip()) for p in args.processors.split(',')]
        processor_counts = sorted(set(processor_counts))  # Remove duplicates and sort
    except ValueError:
        print(f"ERROR: Invalid processor counts '{args.processors}'. Use comma-separated integers.")
        return 1
    
    # Parse gene counts (if specified)
    gene_counts = None
    if args.genes:
        try:
            gene_counts = [int(g.strip()) for g in args.genes.split(',')]
            gene_counts = sorted(set(gene_counts))  # Remove duplicates and sort
        except ValueError:
            print(f"ERROR: Invalid gene counts '{args.genes}'. Use comma-separated integers.")
            return 1
    
    # Parse timepoint counts (if specified)
    timepoint_counts = None
    if args.timepoints:
        try:
            timepoint_counts = [int(t.strip()) for t in args.timepoints.split(',')]
            timepoint_counts = sorted(set(timepoint_counts))  # Remove duplicates and sort
        except ValueError:
            print(f"ERROR: Invalid timepoint counts '{args.timepoints}'. Use comma-separated integers.")
            return 1
    
    # Validate GEP script if provided
    if args.gep_script and not os.path.exists(args.gep_script):
        print(f"ERROR: GEP script not found: {args.gep_script}")
        return 1
    
    print("="*80)
    print("GEP BOOLEAN NETWORK INFERENCE BENCHMARK")
    print("="*80)
    print(f"Binary data dir:  {args.input_binary}")
    print(f"Raw data dir:     {args.input_raw}")
    print(f"Output dir:       {args.output}")
    print(f"Processors:       {', '.join(map(str, processor_counts))}")
    if gene_counts:
        print(f"Gene counts:      {', '.join(map(str, gene_counts))}")
    else:
        print(f"Gene counts:      all available")
    if timepoint_counts:
        print(f"Timepoint counts: {', '.join(map(str, timepoint_counts))}")
    else:
        print(f"Timepoint counts: all available")
    print(f"Population:       {args.population}")
    print(f"Generations:      {args.generations}")
    print(f"Repetitions:      {args.n_repetitions}")
    if args.target_genes:
        print(f"Target genes:     {args.target_genes}")
    if args.gep_script:
        print(f"GEP script:       {args.gep_script}")
    print("="*80)
    
    # Run benchmark
    try:
        results = benchmark_datasets(
            binary_dir=args.input_binary,
            raw_dir=args.input_raw,
            output_base_dir=args.output,
            processor_counts=processor_counts,
            gene_counts=gene_counts,
            timepoint_counts=timepoint_counts,
            population=args.population,
            generations=args.generations,
            n_repetitions=args.n_repetitions,
            target_genes=args.target_genes,
            gep_script=args.gep_script,
            verbose=args.verbose
        )
        
        if results:
            # Generate summary report
            generate_summary_report(results, args.output)
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\nBenchmark interrupted by user.")
        return 1
    except Exception as e:
        print(f"\n\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())