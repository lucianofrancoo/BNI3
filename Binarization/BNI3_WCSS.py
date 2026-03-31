#!/usr/bin/env python3
"""
WCSS Binarization Tool
Implementation of the WCSS (Within-Cluster Sum of Squares) algorithm for gene expression matrix binarization.
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
import time


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def discretize_row_wcss(row, gene_name, verbose):
    """
    Apply WCSS algorithm to a specific row
    
    The WCSS algorithm finds the optimal split of values into two clusters
    by minimizing the Within-Cluster Sum of Squares.
    """
    log_message(f"\n### PROCESSING {gene_name} ###", verbose)
    log_message(f"Original values: {', '.join(map(str, row))}", verbose)
    
    # Sort the row values
    sorted_row = np.sort(row)
    num_cols = len(sorted_row)
    
    log_message(f"Sorted values: {', '.join(map(str, sorted_row))}", verbose)
    
    # Initialize best WCSS as infinity
    best_wcss = np.inf
    best_cluster = None
    best_k = None
    
    # Find the best division into two clusters
    for k in range(1, num_cols):
        # Divide the sorted row into two clusters
        cluster_1 = sorted_row[:k]
        cluster_2 = sorted_row[k:]
        
        # Calculate the mean of each cluster
        mean_1 = np.mean(cluster_1)
        mean_2 = np.mean(cluster_2)
        
        # Calculate the WCSS as the sum of squared Euclidean distances to the mean
        # WCSS = sum of (value - mean)^2 for all values in the cluster
        wcss = 0
        wcss_details_1 = []
        wcss_details_2 = []
        
        for value in cluster_1:
            dist = (value - mean_1) ** 2
            wcss += dist
            wcss_details_1.append(f"{value:.4f}→{dist:.4f}")
        
        for value in cluster_2:
            dist = (value - mean_2) ** 2
            wcss += dist
            wcss_details_2.append(f"{value:.4f}→{dist:.4f}")
        
        log_message(f"  Split at k={k}: cluster_1={list(cluster_1)}, cluster_2={list(cluster_2)}", verbose)
        log_message(f"    Mean_1={mean_1:.4f}, Mean_2={mean_2:.4f}", verbose)
        log_message(f"    Euclidean distances cluster_1: [{', '.join(wcss_details_1)}]", verbose)
        log_message(f"    Euclidean distances cluster_2: [{', '.join(wcss_details_2)}]", verbose)
        log_message(f"    WCSS={wcss:.4f}", verbose)
        
        # Update best WCSS and best cluster if a better division is found
        if wcss < best_wcss:
            best_wcss = wcss
            best_cluster = {'cluster_1': cluster_1.copy(), 'cluster_2': cluster_2.copy()}
            best_k = k
            log_message(f"    → NEW BEST! (WCSS={best_wcss:.4f})", verbose)
    
    log_message(f"\n=== BEST SPLIT FOR {gene_name} ===", verbose)
    log_message(f"Best k: {best_k}", verbose)
    log_message(f"Best WCSS: {best_wcss:.4f}", verbose)
    log_message(f"Cluster 1 (→0): {list(best_cluster['cluster_1'])}", verbose)
    log_message(f"Cluster 2 (→1): {list(best_cluster['cluster_2'])}", verbose)
    
    # Replace values in the original row with 0 and 1 according to the best cluster
    result = np.where(np.isin(row, best_cluster['cluster_1']), 0, 1)
    
    log_message(f"Result: {result}", verbose)
    
    return result


def process_individual_gene(args):
    """
    Helper function to process an individual gene (used in parallelization)
    """
    gene_data, gene_name, verbose = args
    result = discretize_row_wcss(gene_data, gene_name, verbose)
    return gene_name, result


def binarize_matrix(input_file, output_file, verbose, n_processors):
    """
    Process complete file applying WCSS to each gene
    """
    start_time = time.time()
    
    try:
        # Read file
        log_message(f"Reading file: {input_file}", verbose)
        df = pd.read_csv(input_file, sep='\t', index_col='ID')
        
        log_message(f"Matrix dimensions: {df.shape}", verbose)
        
        # Determine number of processors to use
        if n_processors is None:
            n_processors = 1  # Default to 1 processor
        
        # Validate number of processors
        n_processors = max(1, min(n_processors, cpu_count(), df.shape[0]))
        
        log_message(f"Using {n_processors} processors", verbose)
        
        # Create result matrix
        discretized_df = df.copy()
        
        # Processing start time
        processing_start_time = time.time()
        
        if n_processors == 1:
            # Sequential processing
            log_message("Sequential processing", verbose)
            for gene_name, row in df.iterrows():
                result = discretize_row_wcss(row.values, gene_name, verbose)
                discretized_df.loc[gene_name] = result
        else:
            # Parallel processing
            log_message("Parallel processing", verbose)
            
            # Prepare data for parallelization
            # In parallel mode, only disable detailed verbose with many processors
            # For traceability, keep verbose with few processors
            verbose_parallel = verbose if n_processors <= 4 else False
            
            if verbose and n_processors > 4:
                log_message("Note: Detailed verbose disabled in parallel mode to avoid mixed outputs", verbose)
            
            args_list = [(row.values, gene_name, verbose_parallel) 
                        for gene_name, row in df.iterrows()]
            
            # Process in parallel
            with Pool(processes=n_processors) as pool:
                results = pool.map(process_individual_gene, args_list)
            
            # Assign results
            for gene_name, result in results:
                discretized_df.loc[gene_name] = result
        
        processing_time = time.time() - processing_start_time
        
        # Transpose for final format
        final_df = discretized_df.T
        
        # Convert to integers to avoid decimal values (1.0 -> 1, 0.0 -> 0)
        final_df = final_df.astype(int)
        
        # Create output directory if it doesn't exist
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save result
        log_message(f"Saving result to: {output_file}", verbose)
        final_df.to_csv(output_file, sep='\t', index=False)
        
        total_time = time.time() - start_time
        
        # FINAL SUMMARY - ALWAYS DISPLAYED
        print("\n" + "="*60)
        print("WCSS BINARIZATION SUMMARY")
        print("="*60)
        print(f"Input file:        {os.path.basename(input_file)}")
        print(f"Output file:       {os.path.basename(output_file)}")
        print(f"Genes processed:   {df.shape[0]}")
        print(f"Samples:           {df.shape[1]}")
        print(f"Processors used:   {n_processors}")
        print(f"Processing time:   {processing_time:.2f} seconds")
        print(f"Total time:        {total_time:.2f} seconds")
        print(f"Speed:             {df.shape[0]/processing_time:.1f} genes/sec")
        print("Status:            SUCCESS")
        print("="*60)
        
        return final_df
        
    except Exception as e:
        total_time = time.time() - start_time
        print(f"\nERROR: Failed to process file after {total_time:.2f} seconds: {str(e)}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Gene expression matrix binarization using WCSS algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 BNI3_bin_wcss.py -i input.tsv -o output.tsv
  python3 BNI3_bin_wcss.py -i data.tsv -o binary_data.tsv --verbose
  python3 BNI3_bin_wcss.py -i data.tsv -o binary_data.tsv -p 4
  python3 BNI3_bin_wcss.py -i data.tsv -o binary_data.tsv -v -p 8
        """
    )
    
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Input file (expression matrix in TSV format)')
    
    parser.add_argument('-o', '--output',
                        required=True, 
                        help='Output file (binarized matrix)')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Show detailed processing information')
    
    parser.add_argument('-p', '--processors',
                        type=int,
                        default=1,
                        help='Number of processors to use (default: 1)')
    
    parser.add_argument('--version',
                        action='version',
                        version='WCSS Binarizer v1.0')
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"ERROR: Input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Process file
    try:
        binarize_matrix(args.input, args.output, args.verbose, args.processors)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()