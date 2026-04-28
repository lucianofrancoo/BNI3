#!/usr/bin/env python3
"""
SSD Binarization Tool
Implementation of the SSD (Short Series Discretization) algorithm for gene expression matrix binarization.
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from pathlib import Path
from multiprocessing import Pool, cpu_count
import functools
import time


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def detect_isolated_elements(distance_matrix, verbose):
    """
    Detect elements that are completely isolated
    (all their non-diagonal connections are NaN)
    """
    isolated_elements = []
    
    for j in range(len(distance_matrix)):
        # Get all distances excluding diagonal
        element_distances = np.concatenate([
            distance_matrix[j, :j], 
            distance_matrix[j, j+1:]
        ])
        
        if np.all(np.isnan(element_distances)):
            isolated_elements.append(j)
    
    if verbose and isolated_elements:
        log_message(f"Isolated elements detected at positions: {isolated_elements}", verbose)
        log_message("Individual verification:", verbose)
        for elem in isolated_elements:
            elem_distances = np.concatenate([
                distance_matrix[elem, :elem], 
                distance_matrix[elem, elem+1:]
            ])
            log_message(f"  Element {elem} - All NaN? {np.all(np.isnan(elem_distances))}", verbose)
    
    return isolated_elements


def determine_extreme_groups(isolated_elements, connected_elements, ordered_row, verbose):
    """
    Determine which groups go to 0 and which go to 1 based on extreme values
    """
    log_message("=== DETERMINING EXTREME GROUPS ===", verbose)
    log_message(f"Isolated elements: {isolated_elements}", verbose)
    log_message(f"Connected elements: {connected_elements}", verbose)
    
    # Identify extreme values
    min_isolated_idx = min(isolated_elements)
    max_isolated_idx = max(isolated_elements)
    min_connected_idx = min(connected_elements)
    max_connected_idx = max(connected_elements)
    
    min_isolated_value = ordered_row[min_isolated_idx]
    min_connected_value = ordered_row[min_connected_idx]
    max_isolated_value = ordered_row[max_isolated_idx]
    max_connected_value = ordered_row[max_connected_idx]
    
    # DECISION: Which group goes to 0 and which to 1?
    if min_isolated_value < min_connected_value:
        # Case 1: Isolated elements are smaller than connected ones
        # The smallest isolated goes to 0, connected ones go to 1
        group_0 = [min_isolated_idx]
        group_1 = connected_elements.copy()
        intermediates = [x for x in isolated_elements if x != min_isolated_idx]
        
        log_message("CASE 1: Isolated smaller → Smallest isolated goes to 0, connected go to 1", verbose)
        log_message(f"Decision boundary: {min_isolated_value} (smallest isolated) vs {min_connected_value} (smallest connected)", verbose)
    else:
        # Case 2: Connected elements are smaller than isolated ones
        # Connected ones go to 0, largest isolated goes to 1
        group_0 = connected_elements.copy()
        group_1 = [max_isolated_idx]
        intermediates = [x for x in isolated_elements if x != max_isolated_idx]
        
        log_message("CASE 2: Connected smaller → Connected go to 0, largest isolated goes to 1", verbose)
        log_message(f"Decision boundary: {max_connected_value} (largest connected) vs {min_isolated_value} (smallest isolated)", verbose)
    
    log_message(f"Initial group 0: {[ordered_row[i] for i in group_0]}", verbose)
    log_message(f"Initial group 1: {[ordered_row[i] for i in group_1]}", verbose)
    log_message(f"Intermediate elements to evaluate: {[ordered_row[i] for i in intermediates]}", verbose)
    
    return {
        'group_0': group_0,
        'group_1': group_1,
        'intermediates': intermediates
    }


def assign_intermediates_by_proximity(initial_groups, ordered_row, verbose):
    """
    Assign intermediate elements by proximity using dynamic evaluation
    """
    group_0 = initial_groups['group_0']
    group_1 = initial_groups['group_1']
    intermediates = initial_groups['intermediates']
    
    if not intermediates:
        log_message("No intermediate elements to evaluate", verbose)
        return {'group_0': group_0, 'group_1': group_1}
    
    log_message("\n=== ASSIGNING INTERMEDIATE ELEMENTS BY PROXIMITY ===", verbose)
    
    # Evaluate each intermediate element in order
    for elem_idx in intermediates:
        value = ordered_row[elem_idx]
        
        # Calculate current limits of each group
        group_0_limit = max([ordered_row[i] for i in group_0])  # Largest in group 0
        group_1_limit = min([ordered_row[i] for i in group_1])  # Smallest in group 1
        
        # Calculate distances to nearest limits
        dist_to_group_0 = abs(value - group_0_limit)
        dist_to_group_1 = abs(value - group_1_limit)
        
        log_message(f"Evaluating element: {value} (position {elem_idx})", verbose)
        log_message(f"  Current group 0 limit: {group_0_limit} → Distance: {dist_to_group_0}", verbose)
        log_message(f"  Current group 1 limit: {group_1_limit} → Distance: {dist_to_group_1}", verbose)
        
        # DECISION: Assign to nearest group
        if dist_to_group_0 <= dist_to_group_1:
            group_0.append(elem_idx)
            log_message(f"  → Assigned to GROUP 0 (distance {dist_to_group_0} ≤ {dist_to_group_1})", verbose)
            log_message(f"  → Updated group 0: {[ordered_row[i] for i in group_0]}", verbose)
        else:
            group_1.append(elem_idx)
            log_message(f"  → Assigned to GROUP 1 (distance {dist_to_group_1} < {dist_to_group_0})", verbose)
            log_message(f"  → Updated group 1: {[ordered_row[i] for i in group_1]}", verbose)
    
    return {'group_0': group_0, 'group_1': group_1}


def determine_separator(final_groups, ordered_row, verbose):
    """
    Determine final separator based on groups
    """
    group_1 = final_groups['group_1']
    separator = min([ordered_row[i] for i in group_1])
    
    log_message("\n=== FINAL SEPARATOR ===", verbose)
    log_message(f"Final group 0: {[ordered_row[i] for i in final_groups['group_0']]}", verbose)
    log_message(f"Final group 1: {[ordered_row[i] for i in final_groups['group_1']]}", verbose)
    log_message(f"Selected separator: {separator}", verbose)
    
    return separator


def discretize_row_ssd(row, gene_name, verbose):
    """
    Apply SSD algorithm to a specific row
    """
    ordered_row = np.sort(row)
    num_cols = len(ordered_row)
    
    log_message(f"\n\n### PROCESSING {gene_name} ###", verbose)
    log_message(f"Ordered values: {', '.join(map(str, ordered_row))}", verbose)
    
    # Create distance matrix
    distance_matrix = squareform(pdist(ordered_row.reshape(-1, 1)))
    
    separator = None
    isolated_elements = []

    while separator is None:
        # Find maximum distance
        valid_distances = distance_matrix[~np.isnan(distance_matrix)]
        valid_distances = valid_distances[valid_distances > 0]  # Exclude diagonal (0s)
        
        if len(valid_distances) == 0:
            separator = ordered_row[0]
            log_message("No valid distances remaining - using first element as separator", verbose)
            break
        
        max_dist = np.max(valid_distances)
        
        # Find ALL positions with maximum distance
        max_indices = np.where(distance_matrix == max_dist)
        num_max = len(max_indices[0])
        
        log_message(f"Removing maximum distance: {max_dist} at {num_max} positions", verbose)
        
        # Remove ALL maximum distances
        distance_matrix[max_indices] = np.nan
        
        # Detect isolated elements
        isolated_elements = detect_isolated_elements(distance_matrix, verbose)
        
        if isolated_elements:
            if verbose:
                log_message("\n=== FINAL DISTANCE MATRIX ===", verbose)
                # Show rounded matrix only for visualization (calculations use original values)
                display_matrix = np.round(distance_matrix, 3)
                log_message(str(display_matrix), verbose)
                log_message("===============================", verbose)
            
            if len(isolated_elements) == 1:
                # SIMPLE CASE: Only one isolated element
                separator = ordered_row[isolated_elements[0]]
                log_message(f"SIMPLE CASE: Single isolated element → {separator}", verbose)
                
                # Show groups for consistency with complex cases
                if verbose:
                    # Calculate final groups the same way as the main algorithm
                    separator_pos_display = np.where(ordered_row == separator)[0][0]
                    
                    # Apply same distance-based logic for display
                    if separator_pos_display > 0 and separator_pos_display < num_cols - 1:
                        lower_neighbor = ordered_row[separator_pos_display - 1]
                        upper_neighbor = ordered_row[separator_pos_display + 1]
                        dist_to_lower = abs(separator - lower_neighbor)
                        dist_to_upper = abs(separator - upper_neighbor)
                        
                        if dist_to_lower <= dist_to_upper:
                            group_0_values = ordered_row[:separator_pos_display + 1].tolist()
                            group_1_values = ordered_row[separator_pos_display + 1:].tolist()
                        else:
                            group_0_values = ordered_row[:separator_pos_display].tolist()
                            group_1_values = ordered_row[separator_pos_display:].tolist()
                    else:
                        # Handle edge cases
                        if separator_pos_display == 0:
                            group_0_values = [separator]
                            group_1_values = ordered_row[1:].tolist()
                        else:
                            group_0_values = ordered_row[:num_cols-1].tolist()
                            group_1_values = [separator]
                    
                    log_message("\n=== FINAL SEPARATOR ===", verbose)
                    log_message(f"Final group 0: {group_0_values}", verbose)
                    log_message(f"Final group 1: {group_1_values}", verbose)
                    log_message(f"Selected separator: {separator}", verbose)
            
            else:
                # COMPLEX CASE: Multiple isolated elements
                connected_elements = [i for i in range(num_cols) if i not in isolated_elements]
                
                if not connected_elements:
                    # SUBCASE: All elements are isolated
                    separator_idx = isolated_elements[len(isolated_elements) // 2]
                    separator = ordered_row[separator_idx]
                    log_message(f"SUBCASE: All isolated → using middle element: {separator}", verbose)
                    
                    # Show groups for consistency
                    if verbose:
                        separator_pos = np.where(ordered_row == separator)[0][0]
                        if separator_pos == 0:
                            group_0_values = [separator]
                            group_1_values = ordered_row[1:].tolist()
                        elif separator_pos == num_cols - 1:
                            group_0_values = ordered_row[:num_cols-1].tolist()
                            group_1_values = [separator]
                        else:
                            group_0_values = ordered_row[:separator_pos].tolist()
                            group_1_values = ordered_row[separator_pos:].tolist()
                        
                        log_message("\n=== FINAL SEPARATOR ===", verbose)
                        log_message(f"Final group 0: {group_0_values}", verbose)
                        log_message(f"Final group 1: {group_1_values}", verbose)
                        log_message(f"Selected separator: {separator}", verbose)
                
                else:
                    # SUBCASE: There are isolated AND connected elements
                    initial_groups = determine_extreme_groups(isolated_elements, connected_elements, ordered_row, verbose)
                    final_groups = assign_intermediates_by_proximity(initial_groups, ordered_row, verbose)
                    separator = determine_separator(final_groups, ordered_row, verbose)
    
    # Assign 0s and 1s according to separator
    separator_pos = np.where(ordered_row == separator)[0][0]
    
    # For isolated elements, calculate distance to decide group assignment
    if len(isolated_elements) == 1 and separator_pos > 0 and separator_pos < num_cols - 1:
        # Calculate distances to nearest neighbors
        lower_neighbor = ordered_row[separator_pos - 1]
        upper_neighbor = ordered_row[separator_pos + 1]
        
        dist_to_lower = abs(separator - lower_neighbor)
        dist_to_upper = abs(separator - upper_neighbor)
        
        log_message(f"Distance-based assignment for isolated element {separator}:", verbose)
        log_message(f"  Distance to lower neighbor ({lower_neighbor}): {dist_to_lower}", verbose)
        log_message(f"  Distance to upper neighbor ({upper_neighbor}): {dist_to_upper}", verbose)
        
        if dist_to_lower <= dist_to_upper:
            # Closer to lower group - separator goes to group 0
            values_0 = ordered_row[:separator_pos + 1].tolist()
            values_1 = ordered_row[separator_pos + 1:].tolist()
            log_message(f"  → Assigned to lower group (distance {dist_to_lower} ≤ {dist_to_upper})", verbose)
        else:
            # Closer to upper group - separator goes to group 1
            values_0 = ordered_row[:separator_pos].tolist()
            values_1 = ordered_row[separator_pos:].tolist()
            log_message(f"  → Assigned to upper group (distance {dist_to_upper} < {dist_to_lower})", verbose)
    else:
        # Handle extreme cases or use fixed rule for complex cases
        if separator_pos == 0:
            # Separator is the first element
            values_0 = [separator]
            values_1 = ordered_row[1:].tolist()
        elif separator_pos == num_cols - 1:
            # Separator is the last element
            values_0 = ordered_row[:num_cols-1].tolist()
            values_1 = [separator]
        else:
            # Normal case (including complex cases with multiple isolated elements)
            values_0 = ordered_row[:separator_pos].tolist()
            values_1 = ordered_row[separator_pos:].tolist()
    
    # Create binary vector
    result = np.where(np.isin(row, values_0), 0, 1)
    
    return result


def process_individual_gene(args):
    """
    Helper function to process an individual gene (used in parallelization)
    """
    gene_data, gene_name, verbose = args
    result = discretize_row_ssd(gene_data, gene_name, verbose)
    return gene_name, result


def binarize_matrix(input_file, output_file, verbose, n_processors):
    """
    Process complete file applying SSD to each gene
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
                result = discretize_row_ssd(row.values, gene_name, verbose)
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
        print("SSD BINARIZATION SUMMARY")
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
        description="Gene expression matrix binarization using SSD algorithm",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 BNI3_bin.py -i input.tsv -o output.tsv
  python3 BNI3_bin.py -i data.tsv -o binary_data.tsv --verbose
  python3 BNI3_bin.py -i data.tsv -o binary_data.tsv -p 4
  python3 BNI3_bin.py -i data.tsv -o binary_data.tsv -v -p 8
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
                        version='SSD Binarizer v1.0')
    
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