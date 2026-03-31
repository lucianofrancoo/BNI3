#!/usr/bin/env python3
"""
BoolNet Rules Evaluator with Attractor Metrics
Evaluates BoolNet-generated Boolean rules using attractor analysis and parsimony metrics.

This script processes multiple BoolNet rule files (one rule per gene, no combinatorics)
and generates a consolidated report with all evaluation metrics.

Author: Luciano
"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import sys
import os
import glob
from datetime import datetime


def convert_boolnet_to_python(rule: str) -> str:
    """
    Convert BoolNet syntax to Python Boolean syntax
    
    BoolNet uses: ! (NOT), & (AND), | (OR)
    Python uses: ~ (NOT), & (AND), | (OR)
    
    Args:
        rule: BoolNet Boolean rule string
        
    Returns:
        Python-compatible Boolean rule string
    """
    # Replace ! with ~ for NOT operator
    python_rule = rule.replace('!', '~')
    return python_rule


def evaluate_rule(rule: str, gene_state: Dict[str, bool]) -> bool:
    """
    Evaluate a Boolean rule given a gene state
    
    Args:
        rule: Boolean rule string (Python syntax)
        gene_state: Dictionary mapping gene names to boolean values
        
    Returns:
        Result of rule evaluation
    """
    rule_eval = rule.replace('&', ' and ').replace('|', ' or ').replace('~', ' not ')
    
    for gene, value in gene_state.items():
        rule_eval = rule_eval.replace(gene, str(value))
    
    try:
        result = eval(rule_eval)
        return bool(result)
    except:
        return False


def calculate_next_state(current_state: List[bool], 
                         gene_rules: Dict[str, str], 
                         genes: List[str]) -> List[bool]:
    """
    Calculate next state using synchronous Boolean dynamics
    
    Args:
        current_state: Current state as list of booleans
        gene_rules: Dictionary of gene rules
        genes: List of gene names in order
        
    Returns:
        Next state as list of booleans
    """
    gene_state = {gene: current_state[i] for i, gene in enumerate(genes)}
    next_state = []
    
    for gene in genes:
        if gene in gene_rules:
            rule = gene_rules[gene]
            next_value = evaluate_rule(rule, gene_state)
        else:
            next_value = gene_state[gene]
        next_state.append(next_value)
    
    return next_state


def state_to_tuple(state: List[bool]) -> Tuple[bool, ...]:
    """Convert state list to hashable tuple"""
    return tuple(state)


def count_literals_in_rule(rule: str) -> int:
    """
    Count number of gene literals in a Boolean rule
    
    Args:
        rule: Boolean rule string (e.g., "G1 & G2 | ~G3")
        
    Returns:
        Number of gene mentions (literals)
    """
    # Remove operators and parentheses
    cleaned = rule.replace('&', ' ').replace('|', ' ').replace('~', ' ')
    cleaned = cleaned.replace('(', ' ').replace(')', ' ')
    cleaned = cleaned.replace('!', ' ')  # In case original has !
    
    # Split and count non-empty tokens that look like genes
    tokens = [t.strip() for t in cleaned.split() if t.strip()]
    # Filter out boolean constants
    gene_tokens = [t for t in tokens if t not in ['True', 'False', 'true', 'false']]
    
    return len(gene_tokens)


def count_not_operators(rule: str) -> int:
    """
    Count number of NOT operators in a Boolean rule
    Handles both ~ (Python) and ! (BoolNet) syntax
    
    Args:
        rule: Boolean rule string
        
    Returns:
        Number of NOT operators
    """
    return rule.count('~') + rule.count('!')


def extract_regulators(rule: str) -> Set[str]:
    """
    Extract unique gene names that appear in a rule
    
    Args:
        rule: Boolean rule string
        
    Returns:
        Set of unique gene names
    """
    # Remove operators and parentheses
    cleaned = rule.replace('&', ' ').replace('|', ' ').replace('~', ' ')
    cleaned = cleaned.replace('(', ' ').replace(')', ' ')
    cleaned = cleaned.replace('!', ' ')  # In case original has !
    
    # Split and get unique genes
    tokens = [t.strip() for t in cleaned.split() if t.strip()]
    # Filter out boolean constants
    genes = {t for t in tokens if t not in ['True', 'False', 'true', 'false']}
    
    return genes


def calculate_parsimony_metrics(gene_rules: Dict[str, str]) -> Dict[str, float]:
    """
    Calculate parsimony-based metrics for rule complexity
    
    Metrics:
    1. Total literals: Sum of all gene mentions across all rules (prefer fewer)
    2. Total NOT operators: Sum of all negations (prefer fewer - less repression)
    3. Average K (connectivity): Average number of regulators per gene (prefer K≈2)
    
    Args:
        gene_rules: Dictionary mapping gene names to their Boolean rules
        
    Returns:
        Dictionary with parsimony metrics
    """
    total_literals = 0
    total_nots = 0
    k_values = []
    
    for gene, rule in gene_rules.items():
        # Skip if rule is a constant
        if rule in ['True', 'False', 'true', 'false']:
            total_literals += 0
            total_nots += 0
            k_values.append(0)
            continue
        
        # Count literals
        total_literals += count_literals_in_rule(rule)
        
        # Count NOT operators
        total_nots += count_not_operators(rule)
        
        # Count unique regulators for K
        regulators = extract_regulators(rule)
        k_values.append(len(regulators))
    
    # Calculate statistics
    avg_k = np.mean(k_values) if k_values else 0
    std_k = np.std(k_values) if k_values else 0
    
    # K distance from optimal (Kauffman's K=2)
    k_distance_from_2 = abs(avg_k - 2.0)
    
    return {
        'total_literals': total_literals,
        'total_nots': total_nots,
        'avg_k': avg_k,
        'std_k': std_k,
        'k_distance_from_2': k_distance_from_2
    }


def find_attractors_for_ruleset(gene_rules: Dict[str, str], 
                                 genes: List[str], 
                                 max_iterations: int = 1000) -> Tuple[List[List[List[bool]]], Dict]:
    """
    Find all attractors for a given ruleset
    
    Args:
        gene_rules: Dictionary of gene rules
        genes: List of gene names
        max_iterations: Maximum iterations for trajectory simulation
        
    Returns:
        Tuple of (list of attractors, basin sizes dictionary)
    """
    n_genes = len(genes)
    total_states = 2 ** n_genes
    
    visited_states = set()
    attractors = []
    basins = defaultdict(int)
    
    # Test all possible initial states
    for i in range(total_states):
        # Generate binary state
        binary = format(i, f'0{n_genes}b')
        initial_state = [bool(int(bit)) for bit in binary]
        initial_tuple = state_to_tuple(initial_state)
        
        if initial_tuple in visited_states:
            continue
        
        # Simulate trajectory
        trajectory = [initial_state]
        trajectory_set = {initial_tuple}
        current_state = initial_state
        
        for _ in range(max_iterations):
            next_state = calculate_next_state(current_state, gene_rules, genes)
            next_tuple = state_to_tuple(next_state)
            
            if next_tuple in trajectory_set:
                # Found a cycle - extract attractor
                cycle_start = next((i for i, s in enumerate(trajectory) 
                                   if state_to_tuple(s) == next_tuple))
                attractor = trajectory[cycle_start:]
                
                # Check if this attractor is new
                attractor_tuples = [state_to_tuple(s) for s in attractor]
                is_new = True
                attractor_id = -1
                
                for idx, existing_attractor in enumerate(attractors):
                    existing_tuples = [state_to_tuple(s) for s in existing_attractor]
                    if set(attractor_tuples) == set(existing_tuples):
                        is_new = False
                        attractor_id = idx
                        break
                
                if is_new:
                    attractors.append(attractor)
                    attractor_id = len(attractors) - 1
                
                # Mark all states in trajectory as visited and assign to basin
                for state in trajectory:
                    visited_states.add(state_to_tuple(state))
                basins[attractor_id] += 1
                
                break
            
            trajectory.append(next_state)
            trajectory_set.add(next_tuple)
            current_state = next_state
    
    return attractors, basins


def calculate_attractor_metrics(attractors: List[List[List[bool]]], 
                                basins: Dict[int, int],
                                genes: List[str],
                                binarized_matrix: pd.DataFrame = None) -> Dict[str, float]:
    """
    Calculate metrics for attractor analysis
    
    Args:
        attractors: List of attractors
        basins: Basin sizes for each attractor
        genes: List of gene names
        binarized_matrix: Optional binarized expression matrix
        
    Returns:
        Dictionary with metric scores
    """
    metrics = {}
    
    # Metric 1: Number of basins (fewer is better, more interpretable)
    n_basins = len(attractors)
    metrics['n_basins'] = n_basins
    metrics['basin_penalty'] = n_basins
    
    # Metric 2: Attractor complexity (cycle length)
    cycle_lengths = [len(att) for att in attractors]
    avg_cycle_length = np.mean(cycle_lengths) if cycle_lengths else 0
    metrics['avg_cycle_length'] = avg_cycle_length
    metrics['complexity_penalty'] = avg_cycle_length
    
    # Metric 3: Basin distribution (entropy - uniformity)
    total_states = sum(basins.values())
    basin_probs = [basins[i] / total_states for i in range(len(attractors))]
    if basin_probs:
        entropy = -sum(p * np.log2(p) if p > 0 else 0 for p in basin_probs)
        metrics['basin_entropy'] = entropy
    else:
        metrics['basin_entropy'] = 0
    
    # Metric 4: Fixed point concordance with final timepoint (if matrix provided)
    metrics['final_state_concordance'] = 0
    if binarized_matrix is not None:
        # Check if any fixed point matches the last timepoint
        last_timepoint = binarized_matrix.iloc[-1]
        
        for attractor in attractors:
            if len(attractor) == 1:  # Fixed point
                fixed_point = attractor[0]
                # Compare with last timepoint
                matches = sum(1 for i, gene in enumerate(genes) 
                             if gene in last_timepoint.index and 
                             fixed_point[i] == bool(last_timepoint[gene]))
                concordance = matches / len(genes)
                metrics['final_state_concordance'] = max(
                    metrics['final_state_concordance'], 
                    concordance
                )
    
    # Metric 5: Proportion of fixed points (higher is better)
    fixed_points = sum(1 for att in attractors if len(att) == 1)
    metrics['fixed_point_ratio'] = fixed_points / n_basins if n_basins > 0 else 0
    
    return metrics


def load_boolnet_rules_table(rules_file: str) -> pd.DataFrame:
    """Load BoolNet rules table from TSV file"""
    df = pd.read_csv(rules_file, sep='\t', encoding='utf-8')
    required_cols = ['Gen', 'Regla']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}. Found: {df.columns.tolist()}")
    return df


def load_binarized_matrix(matrix_file: str) -> pd.DataFrame:
    """Load binarized expression matrix"""
    df = pd.read_csv(matrix_file, sep='\t')
    
    # Check if first column looks like gene names or timepoints
    first_col = df.columns[0]
    
    # If first column doesn't look like a gene name, set it as index
    if not any(df.columns[0].startswith(prefix) for prefix in ['G', 'AT', 'LOC', 'ENSG']):
        df = df.set_index(df.columns[0])
    
    return df


def evaluate_boolnet_ruleset(rules_file: str,
                             binarized_matrix_file: str = None,
                             max_iterations: int = 1000,
                             verbose: bool = False) -> Dict:
    """
    Evaluate a single BoolNet ruleset
    
    Args:
        rules_file: Path to BoolNet rules TSV file
        binarized_matrix_file: Optional path to binarized expression matrix
        max_iterations: Maximum iterations for attractor search
        verbose: Print progress information
        
    Returns:
        Dictionary with evaluation results
    """
    # Load rules
    df = load_boolnet_rules_table(rules_file)
    
    # Load binarized matrix if provided
    binarized_matrix = None
    if binarized_matrix_file and os.path.exists(binarized_matrix_file):
        binarized_matrix = load_binarized_matrix(binarized_matrix_file)
    
    # Build gene_rules dictionary (convert BoolNet syntax to Python)
    gene_rules = {}
    genes = []
    
    for _, row in df.iterrows():
        gene = str(row['Gen'])
        boolnet_rule = str(row['Regla'])
        python_rule = convert_boolnet_to_python(boolnet_rule)
        
        gene_rules[gene] = python_rule
        genes.append(gene)
    
    genes = sorted(genes)
    
    # Find attractors
    attractors, basins = find_attractors_for_ruleset(gene_rules, genes, max_iterations)
    
    # Calculate attractor metrics
    metrics = calculate_attractor_metrics(attractors, basins, genes, binarized_matrix)
    
    # Calculate parsimony metrics
    parsimony_metrics = calculate_parsimony_metrics(gene_rules)
    metrics.update(parsimony_metrics)
    
    # Calculate final integrated score (same formula as BNI3)
    def safe_normalize_single(value, min_val, max_val):
        """Normalize single value to [0,1] range"""
        if max_val == min_val:
            return 0.0
        return (value - min_val) / (max_val - min_val)
    
    # For single evaluation, we can't normalize across multiple samples
    # So we'll just store the raw metrics and calculate the score components
    # The normalization will happen when we consolidate all results
    
    # Store all metrics
    result = {
        'rules_file': os.path.basename(rules_file),
        'n_genes': len(genes),
        'n_basins': metrics['n_basins'],
        'avg_cycle_length': metrics['avg_cycle_length'],
        'basin_entropy': metrics['basin_entropy'],
        'final_state_concordance': metrics['final_state_concordance'],
        'fixed_point_ratio': metrics['fixed_point_ratio'],
        'total_literals': metrics['total_literals'],
        'total_nots': metrics['total_nots'],
        'avg_k': metrics['avg_k'],
        'std_k': metrics['std_k'],
        'k_distance_from_2': metrics['k_distance_from_2']
    }
    
    return result


def find_all_boolnet_rules_files(base_dir: str) -> List[Tuple[str, str, str]]:
    """
    Find all BoolNet rules files in subdirectories
    
    Args:
        base_dir: Base directory to search (e.g., /path/to/size10)
        
    Returns:
        List of tuples: (rules_file_path, binarized_matrix_path, subdirectory_name)
    """
    files_info = []
    
    # Find all subdirectories (size10_1, size10_2, etc.)
    subdirs = sorted(glob.glob(os.path.join(base_dir, 'size10_*')))
    
    for subdir in subdirs:
        subdir_name = os.path.basename(subdir)
        
        # Find all reglas_Boolnet_*.tsv files
        rules_files = sorted(glob.glob(os.path.join(subdir, 'reglas_Boolnet_*.tsv')))
        
        for rules_file in rules_files:
            # Extract the base name to find corresponding binarized matrix
            # reglas_Boolnet_bestfit_size10_1_1.tsv -> bin_size10_1_1.tsv
            rules_basename = os.path.basename(rules_file)
            
            # Remove 'reglas_Boolnet_bestfit_' prefix and '.tsv' suffix
            if rules_basename.startswith('reglas_Boolnet_bestfit_'):
                matrix_suffix = rules_basename.replace('reglas_Boolnet_bestfit_', '')
                matrix_file = os.path.join(subdir, f'bin_{matrix_suffix}')
            else:
                # Fallback: try to guess
                matrix_file = None
            
            files_info.append((rules_file, matrix_file, subdir_name))
    
    return files_info


def evaluate_all_boolnet_rules(base_dir: str,
                               output_file: str,
                               max_iterations: int = 1000,
                               verbose: bool = True):
    """
    Main function to evaluate all BoolNet rules in subdirectories
    
    Args:
        base_dir: Base directory containing size10_* subdirectories
        output_file: Path to consolidated output TSV file
        max_iterations: Maximum iterations for attractor search
        verbose: Print progress information
    """
    # Create log file
    log_file = output_file.replace('.tsv', '_log.txt')
    log_lines = []
    
    def log(message):
        """Helper to log both to console and file"""
        if verbose:
            print(message)
        log_lines.append(message)
    
    log("="*80)
    log("BOOLNET RULES EVALUATION WITH ATTRACTOR METRICS")
    log("="*80)
    log(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Base directory: {base_dir}")
    log(f"Output file: {output_file}")
    log("")
    
    # Find all rules files
    log("1. Scanning for BoolNet rules files...")
    files_info = find_all_boolnet_rules_files(base_dir)
    
    if not files_info:
        log("ERROR: No BoolNet rules files found!")
        log(f"Searched in: {base_dir}/size10_*/reglas_Boolnet_*.tsv")
        return
    
    log(f"   Found {len(files_info)} rules files to evaluate")
    
    # Count by subdirectory
    subdirs_count = {}
    for _, _, subdir in files_info:
        subdirs_count[subdir] = subdirs_count.get(subdir, 0) + 1
    
    log(f"   Distribution:")
    for subdir in sorted(subdirs_count.keys()):
        log(f"      {subdir}: {subdirs_count[subdir]} files")
    log("")
    
    # Evaluate all rulesets
    log("2. Evaluating rulesets...")
    results = []
    
    for i, (rules_file, matrix_file, subdir) in enumerate(files_info, 1):
        rules_basename = os.path.basename(rules_file)
        
        if verbose and i % 5 == 0:
            log(f"   Progress: {i}/{len(files_info)} files evaluated")
        
        try:
            result = evaluate_boolnet_ruleset(
                rules_file=rules_file,
                binarized_matrix_file=matrix_file,
                max_iterations=max_iterations,
                verbose=False
            )
            
            # Add subdirectory info
            result['subdirectory'] = subdir
            result['full_path'] = rules_file
            
            # Extract additional metadata from filename
            # reglas_Boolnet_bestfit_size10_1_1.tsv -> network: size10_1, replicate: 1
            basename = rules_basename.replace('reglas_Boolnet_bestfit_', '').replace('.tsv', '')
            
            # Check if it's a _5times variant
            if '_5times' in basename:
                is_5times = True
                basename = basename.replace('_5times', '')
            else:
                is_5times = False
            
            result['is_5times'] = is_5times
            result['variant'] = '5times' if is_5times else 'normal'
            
            # Parse network and replicate
            # Expected format: size10_X_Y where X is network, Y is replicate
            parts = basename.split('_')
            if len(parts) >= 3:
                result['network_id'] = f"size10_{parts[1]}"
                result['replicate'] = parts[2] if len(parts) > 2 else 'unknown'
            else:
                result['network_id'] = basename
                result['replicate'] = 'unknown'
            
            results.append(result)
            
        except Exception as e:
            log(f"   ERROR evaluating {rules_basename}: {str(e)}")
            continue
    
    log(f"   Completed: {len(results)}/{len(files_info)} files successfully evaluated")
    log("")
    
    if not results:
        log("ERROR: No results to save!")
        return
    
    # Create results dataframe
    log("3. Processing and calculating final scores...")
    results_df = pd.DataFrame(results)
    
    # Calculate normalized final score across ALL results
    def safe_normalize(series):
        """Normalize series to 0-1, handling edge cases"""
        max_val = series.max()
        min_val = series.min()
        if max_val == min_val:
            return pd.Series([0.0] * len(series), index=series.index)
        return (series - min_val) / (max_val - min_val)
    
    # Normalize each metric
    norm_n_basins = safe_normalize(results_df['n_basins'])
    norm_cycle_length = safe_normalize(results_df['avg_cycle_length'])
    norm_basin_entropy = safe_normalize(results_df['basin_entropy'])
    # For these, higher is better, so invert
    norm_concordance = 1 - results_df['final_state_concordance']
    norm_fixed_ratio = 1 - results_df['fixed_point_ratio']
    
    # Parsimony metrics
    norm_literals = safe_normalize(results_df['total_literals'])
    norm_nots = safe_normalize(results_df['total_nots'])
    norm_k_distance = safe_normalize(results_df['k_distance_from_2'])
    
    # Calculate FINAL UNIFIED SCORE (same weights as BNI3)
    results_df['final_score'] = (
        # Attractor metrics (60% total weight)
        0.20 * norm_n_basins +
        0.15 * norm_cycle_length +
        0.10 * norm_basin_entropy +
        0.10 * norm_concordance +
        0.05 * norm_fixed_ratio +
        
        # Parsimony metrics (40% total weight)
        0.20 * norm_k_distance +
        0.12 * norm_literals +
        0.08 * norm_nots
    )
    
    # Sort by final_score
    results_df = results_df.sort_values('final_score', ascending=True)
    
    # Reorder columns for better readability
    column_order = [
        'subdirectory',
        'rules_file',
        'network_id',
        'replicate',
        'variant',
        'is_5times',
        'n_genes',
        'final_score',
        'n_basins',
        'avg_cycle_length',
        'basin_entropy',
        'final_state_concordance',
        'fixed_point_ratio',
        'total_literals',
        'total_nots',
        'avg_k',
        'std_k',
        'k_distance_from_2',
        'full_path'
    ]
    
    # Only include columns that exist
    column_order = [col for col in column_order if col in results_df.columns]
    results_df = results_df[column_order]
    
    # Save results
    log("4. Saving results...")
    results_df.to_csv(output_file, sep='\t', index=False)
    log(f"   Saved consolidated results to: {output_file}")
    log("")
    
    # Generate summary statistics
    log("="*80)
    log("EVALUATION SUMMARY")
    log("="*80)
    log(f"Total files evaluated: {len(results_df)}")
    log(f"Total networks: {results_df['network_id'].nunique()}")
    log(f"Variants: normal={len(results_df[~results_df['is_5times']])}, "
        f"5times={len(results_df[results_df['is_5times']])}")
    log("")
    
    log("Best performing ruleset:")
    best = results_df.iloc[0]
    log(f"  File: {best['rules_file']}")
    log(f"  Network: {best['network_id']}, Replicate: {best['replicate']}, Variant: {best['variant']}")
    log(f"  Final score: {best['final_score']:.4f} (LOWER IS BETTER)")
    log(f"  Attractor metrics:")
    log(f"    - Number of basins: {int(best['n_basins'])}")
    log(f"    - Avg cycle length: {best['avg_cycle_length']:.2f}")
    log(f"    - Basin entropy: {best['basin_entropy']:.4f}")
    log(f"    - Fixed point ratio: {best['fixed_point_ratio']:.2f}")
    log(f"    - Final state concordance: {best['final_state_concordance']:.4f}")
    log(f"  Parsimony metrics:")
    log(f"    - Total literals: {int(best['total_literals'])}")
    log(f"    - Total NOTs: {int(best['total_nots'])}")
    log(f"    - Avg K: {best['avg_k']:.2f} (optimal ≈ 2.0)")
    log(f"    - K distance from 2.0: {best['k_distance_from_2']:.3f}")
    log("")
    
    log("Score statistics:")
    log(f"  Mean final_score: {results_df['final_score'].mean():.4f}")
    log(f"  Std final_score: {results_df['final_score'].std():.4f}")
    log(f"  Min final_score: {results_df['final_score'].min():.4f}")
    log(f"  Max final_score: {results_df['final_score'].max():.4f}")
    log("")
    
    log("Comparison by variant:")
    for variant in ['normal', '5times']:
        variant_df = results_df[results_df['variant'] == variant]
        if len(variant_df) > 0:
            log(f"  {variant}:")
            log(f"    Count: {len(variant_df)}")
            log(f"    Mean final_score: {variant_df['final_score'].mean():.4f}")
            log(f"    Best final_score: {variant_df['final_score'].min():.4f}")
            log(f"    Worst final_score: {variant_df['final_score'].max():.4f}")
    log("")
    
    log(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("="*80)
    
    # Save log file
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    
    if verbose:
        print(f"\nLog saved to: {log_file}")


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Evaluate BoolNet-generated Boolean rules using attractor metrics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Evaluate all BoolNet rules in size10 subdirectories
  python3 BNI3_Evaluate_BoolNet_rules.py \\
      -i /home/lahumada/disco1/BNI3/BoolNet_comparison/size10 \\
      -o /home/lahumada/disco1/BNI3/BoolNet_comparison/size10/boolnet_evaluation_results.tsv \\
      -v

Input structure expected:
  base_dir/
    size10_1/
      reglas_Boolnet_bestfit_size10_1_1.tsv
      reglas_Boolnet_bestfit_size10_1_1_5times.tsv
      bin_size10_1_1.tsv
      bin_size10_1_1_5times.tsv
      ...
    size10_2/
      ...

Output:
  - Consolidated TSV with all evaluations (sorted by final_score)
  - Log file with execution summary and statistics

Metrics calculated:
  - final_score: Unified score combining ALL metrics (LOWER IS BETTER)
  - Attractor metrics: n_basins, avg_cycle_length, basin_entropy, etc.
  - Parsimony metrics: total_literals, total_nots, avg_k, k_distance_from_2
        """
    )
    
    parser.add_argument('-i', '--input', type=str, required=True,
                       help='Base directory containing size10_* subdirectories')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Output consolidated TSV file')
    parser.add_argument('--max-iter', type=int, default=1000,
                       help='Maximum iterations for attractor search (default: 1000)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print detailed progress information')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Run evaluation
    try:
        evaluate_all_boolnet_rules(
            base_dir=args.input,
            output_file=args.output,
            max_iterations=args.max_iter,
            verbose=args.verbose
        )
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()