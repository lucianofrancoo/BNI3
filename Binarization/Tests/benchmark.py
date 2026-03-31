#!/usr/bin/env python3
"""
Binarization Algorithm Benchmark Tool
Runs SSD and WCSS binarization algorithms on test datasets and generates performance metrics.
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
    Parse test dataset filename to extract genes and timepoints
    Expected format: test_data_g{n_genes}_t{n_timepoints}.tsv
    """
    pattern = r'test_data_g(\d+)_t(\d+)\.tsv'
    match = re.match(pattern, filename)
    
    if match:
        n_genes = int(match.group(1))
        n_timepoints = int(match.group(2))
        return n_genes, n_timepoints
    
    return None, None


def run_binarization(algorithm, input_file, output_file, n_processors=1, script_path=None):
    """
    Run a binarization algorithm and measure execution time
    
    Args:
        algorithm: Name of the algorithm ('SSD' or 'WCSS')
        input_file: Path to input file
        output_file: Path to output file
        n_processors: Number of processors to use
        script_path: Custom path to the binarization script
    
    Returns:
        dict: Results including time, success status, and algorithm name
    """
    
    # Determine which script to use
    if script_path:
        script = script_path
    else:
        if algorithm.upper() == 'SSD':
            script = 'BNI3_bin_ssd.py'
        elif algorithm.upper() == 'WCSS':
            script = 'BNI3_bin_wcss.py'
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")
    
    # Validate script exists
    if not os.path.exists(script):
        raise FileNotFoundError(f"Script not found: {script}")
    
    # Build command
    cmd = [
        'python3',
        script,
        '-i', str(input_file),
        '-o', str(output_file),
        '-p', str(n_processors)
    ]
    
    print(f"  Running {algorithm}...", end=' ', flush=True)
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        elapsed_time = time.time() - start_time
        
        # Parse processing time from output
        processing_time = None
        for line in result.stdout.split('\n'):
            if 'Processing time:' in line:
                match = re.search(r'(\d+\.?\d*)\s+seconds', line)
                if match:
                    processing_time = float(match.group(1))
                    break
        
        if processing_time is None:
            processing_time = elapsed_time
        
        print(f"✓ ({elapsed_time:.2f}s)")
        
        return {
            'algorithm': algorithm,
            'success': True,
            'total_time': elapsed_time,
            'processing_time': processing_time,
            'error': None
        }
        
    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print(f"✗ (Failed after {elapsed_time:.2f}s)")
        
        return {
            'algorithm': algorithm,
            'success': False,
            'total_time': elapsed_time,
            'processing_time': None,
            'error': str(e)
        }


def benchmark_datasets(input_dir, output_dir, algorithms, n_processors, ssd_script=None, wcss_script=None):
    """
    Run binarization algorithms on all test datasets and collect metrics
    
    Args:
        input_dir: Directory containing test datasets
        output_dir: Directory for output files
        algorithms: List of algorithms to run
        n_processors: Number of processors to use
        ssd_script: Custom path to SSD script (optional)
        wcss_script: Custom path to WCSS script (optional)
    """
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Find all test datasets
    test_files = sorted(input_path.glob('test_data_g*_t*.tsv'))
    
    if not test_files:
        print(f"ERROR: No test datasets found in {input_dir}")
        print("Expected format: test_data_g{n_genes}_t{n_timepoints}.tsv")
        return None
    
    print(f"\nFound {len(test_files)} test datasets")
    print("="*80)
    
    results = []
    
    for test_file in test_files:
        # Parse filename
        n_genes, n_timepoints = parse_filename(test_file.name)
        
        if n_genes is None or n_timepoints is None:
            print(f"WARNING: Skipping file with unexpected name: {test_file.name}")
            continue
        
        print(f"\n{test_file.name} ({n_genes} genes × {n_timepoints} timepoints)")
        
        # Run each algorithm
        for algorithm in algorithms:
            # Create output filename with detailed information
            output_file = output_path / f"bin_{algorithm.lower()}_g{n_genes}_t{n_timepoints}_p{n_processors}_{test_file.name}"
            
            # Determine script path
            script_path = None
            if algorithm.upper() == 'SSD' and ssd_script:
                script_path = ssd_script
            elif algorithm.upper() == 'WCSS' and wcss_script:
                script_path = wcss_script
            
            # Run binarization
            result = run_binarization(
                algorithm=algorithm,
                input_file=test_file,
                output_file=output_file,
                n_processors=n_processors,
                script_path=script_path
            )
            
            # Add metadata
            result.update({
                'input_file': test_file.name,
                'output_file': output_file.name,
                'n_genes': n_genes,
                'n_timepoints': n_timepoints,
                'n_processors': n_processors
            })
            
            # Calculate per-gene and per-timepoint metrics
            if result['success'] and result['processing_time'] is not None:
                result['time_per_gene'] = result['processing_time'] / n_genes
                result['time_per_timepoint'] = result['processing_time'] / n_timepoints
                result['time_per_gene_timepoint'] = result['processing_time'] / (n_genes * n_timepoints)
            else:
                result['time_per_gene'] = None
                result['time_per_timepoint'] = None
                result['time_per_gene_timepoint'] = None
            
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
    results_file = output_path / 'benchmark_results.csv'
    df.to_csv(results_file, index=False, float_format='%.6f')
    
    print("\n" + "="*80)
    print("BENCHMARK SUMMARY")
    print("="*80)
    
    # Overall statistics
    print("\n📊 OVERALL STATISTICS:")
    print("-" * 80)
    
    for algorithm in df['algorithm'].unique():
        alg_data = df[df['algorithm'] == algorithm]
        successful = alg_data['success'].sum()
        total = len(alg_data)
        
        print(f"\n{algorithm} Algorithm:")
        print(f"  Total runs:        {total}")
        print(f"  Successful:        {successful}")
        print(f"  Failed:            {total - successful}")
        
        if successful > 0:
            success_data = alg_data[alg_data['success']]
            print(f"  Total time:        {success_data['total_time'].sum():.2f} seconds")
            print(f"  Avg time/dataset:  {success_data['processing_time'].mean():.4f} seconds")
            print(f"  Min time/dataset:  {success_data['processing_time'].min():.4f} seconds")
            print(f"  Max time/dataset:  {success_data['processing_time'].max():.4f} seconds")
    
    # Performance metrics
    print("\n⚡ PERFORMANCE METRICS:")
    print("-" * 80)
    
    successful_data = df[df['success']]
    
    if len(successful_data) > 0:
        # Group by algorithm
        for algorithm in successful_data['algorithm'].unique():
            alg_data = successful_data[successful_data['algorithm'] == algorithm]
            
            print(f"\n{algorithm} Algorithm:")
            print(f"  Avg time per gene:              {alg_data['time_per_gene'].mean():.6f} seconds")
            print(f"  Avg time per timepoint:         {alg_data['time_per_timepoint'].mean():.6f} seconds")
            print(f"  Avg time per gene×timepoint:    {alg_data['time_per_gene_timepoint'].mean():.8f} seconds")
            print(f"  Throughput (genes/sec):         {(1/alg_data['time_per_gene'].mean()):.2f}")
    
    # Detailed results table
    print("\n📋 DETAILED RESULTS:")
    print("-" * 80)
    print(f"{'Dataset':<25} {'Algorithm':<8} {'Genes':>6} {'Times':>6} {'Time(s)':>10} {'Time/Gene':>12} {'Status':<10}")
    print("-" * 80)
    
    for _, row in df.iterrows():
        status = "✓ Success" if row['success'] else "✗ Failed"
        time_str = f"{row['processing_time']:.4f}" if row['success'] else "N/A"
        time_per_gene_str = f"{row['time_per_gene']:.6f}" if row['success'] else "N/A"
        
        print(f"{row['input_file']:<25} {row['algorithm']:<8} {row['n_genes']:>6} "
              f"{row['n_timepoints']:>6} {time_str:>10} {time_per_gene_str:>12} {status:<10}")
    
    # Algorithm comparison
    print("\n🔄 ALGORITHM COMPARISON:")
    print("-" * 80)
    
    # Group by dataset and compare
    comparison_data = []
    
    for input_file in df['input_file'].unique():
        file_data = successful_data[successful_data['input_file'] == input_file]
        
        if len(file_data) == 2:  # Both algorithms succeeded
            ssd_time = file_data[file_data['algorithm'] == 'SSD']['processing_time'].values[0]
            wcss_time = file_data[file_data['algorithm'] == 'WCSS']['processing_time'].values[0]
            n_genes = file_data['n_genes'].values[0]
            n_timepoints = file_data['n_timepoints'].values[0]
            
            ratio = wcss_time / ssd_time
            faster = "SSD" if ratio > 1 else "WCSS"
            
            comparison_data.append({
                'dataset': input_file,
                'n_genes': n_genes,
                'n_timepoints': n_timepoints,
                'ssd_time': ssd_time,
                'wcss_time': wcss_time,
                'ratio': ratio,
                'faster': faster
            })
    
    if comparison_data:
        comp_df = pd.DataFrame(comparison_data)
        
        print(f"\n{'Dataset':<25} {'SSD(s)':>10} {'WCSS(s)':>10} {'Ratio':>10} {'Faster':<10}")
        print("-" * 80)
        
        for _, row in comp_df.iterrows():
            print(f"{row['dataset']:<25} {row['ssd_time']:>10.4f} {row['wcss_time']:>10.4f} "
                  f"{row['ratio']:>10.2f}x {row['faster']:<10}")
        
        # Save comparison to CSV
        comparison_file = output_path / 'algorithm_comparison.csv'
        comp_df.to_csv(comparison_file, index=False, float_format='%.6f')
        
        print(f"\n📈 Average speed ratio (WCSS/SSD): {comp_df['ratio'].mean():.2f}x")
    
    print("\n" + "="*80)
    print(f"📁 Detailed results saved to: {results_file}")
    if comparison_data:
        print(f"📁 Comparison saved to: {comparison_file}")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark binarization algorithms (SSD and WCSS)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 run_benchmark.py -i test_data -o results
  python3 run_benchmark.py -i ./datasets -o ./binarized -p 4
  python3 run_benchmark.py -i test_data -o results -a SSD
  python3 run_benchmark.py -i test_data -o results -a WCSS
  python3 run_benchmark.py -i test_data -o results --ssd_script /path/to/BNI3_bin_ssd.py
  python3 run_benchmark.py -i test_data -o results --ssd_script ./ssd.py --wcss_script ./wcss.py
        """
    )
    
    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input directory containing test datasets')
    
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory for binarized results and reports')
    
    parser.add_argument('-a', '--algorithms',
                        default='SSD,WCSS',
                        help='Comma-separated list of algorithms to run (default: SSD,WCSS)')
    
    parser.add_argument('-p', '--processors',
                        type=int,
                        default=1,
                        help='Number of processors to use for each run (default: 1)')
    
    parser.add_argument('--ssd_script',
                        default=None,
                        help='Path to SSD binarization script (default: BNI3_bin_ssd.py in current directory)')
    
    parser.add_argument('--wcss_script',
                        default=None,
                        help='Path to WCSS binarization script (default: BNI3_bin_wcss.py in current directory)')
    
    parser.add_argument('--version',
                        action='version',
                        version='Binarization Benchmark v1.0')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.input):
        print(f"ERROR: Input directory '{args.input}' does not exist.")
        return 1
    
    # Parse algorithms
    algorithms = [alg.strip().upper() for alg in args.algorithms.split(',')]
    valid_algorithms = ['SSD', 'WCSS']
    
    for alg in algorithms:
        if alg not in valid_algorithms:
            print(f"ERROR: Unknown algorithm '{alg}'. Valid options: {', '.join(valid_algorithms)}")
            return 1
    
    # Validate script paths if provided
    if args.ssd_script and not os.path.exists(args.ssd_script):
        print(f"ERROR: SSD script not found: {args.ssd_script}")
        return 1
    
    if args.wcss_script and not os.path.exists(args.wcss_script):
        print(f"ERROR: WCSS script not found: {args.wcss_script}")
        return 1
    
    print("="*80)
    print("BINARIZATION ALGORITHM BENCHMARK")
    print("="*80)
    print(f"Input directory:  {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Algorithms:       {', '.join(algorithms)}")
    print(f"Processors:       {args.processors}")
    if args.ssd_script:
        print(f"SSD script:       {args.ssd_script}")
    if args.wcss_script:
        print(f"WCSS script:      {args.wcss_script}")
    print("="*80)
    
    # Run benchmark
    try:
        results = benchmark_datasets(
            input_dir=args.input,
            output_dir=args.output,
            algorithms=algorithms,
            n_processors=args.processors,
            ssd_script=args.ssd_script,
            wcss_script=args.wcss_script
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
        return 1


if __name__ == "__main__":
    sys.exit(main())