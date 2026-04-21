#!/usr/bin/env python3
"""
Boolean Rules Evaluator with Attractor Metrics
Evaluates different combinations of Boolean rules using attractor analysis
without generating intermediate files.

Author: Luciano
"""

import argparse
import pandas as pd
import numpy as np
import itertools
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import sys
import os
import multiprocessing as mp
from functools import partial


def evaluate_rule(rule: str, gene_state: Dict[str, bool]) -> bool:
    """
    Evaluate a Boolean rule given a gene state
    
    Args:
        rule: Boolean rule string
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
    
    # Split and count non-empty tokens that look like genes
    tokens = [t.strip() for t in cleaned.split() if t.strip()]
    # Filter out boolean constants
    gene_tokens = [t for t in tokens if t not in ['True', 'False', 'true', 'false']]
    
    return len(gene_tokens)


def count_not_operators(rule: str) -> int:
    """
    Count number of NOT operators (~) in a Boolean rule
    
    Args:
        rule: Boolean rule string
        
    Returns:
        Number of NOT operators
    """
    return rule.count('~')


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


def state_to_tuple(state: List[bool]) -> Tuple[bool, ...]:
    """Convert state list to hashable tuple"""
    return tuple(state)


def find_attractors_for_ruleset(gene_rules: Dict[str, str], 
                                 genes: List[str], 
                                 max_iterations: int = 1000) -> Tuple[List[List[List[bool]]], Dict]:
    """
    Find all attractors for a given ruleset without parallelization
    
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
    metrics['basin_penalty'] = n_basins  # Penalize more basins
    
    # Metric 2: Attractor complexity (cycle length)
    cycle_lengths = [len(att) for att in attractors]
    avg_cycle_length = np.mean(cycle_lengths) if cycle_lengths else 0
    metrics['avg_cycle_length'] = avg_cycle_length
    metrics['complexity_penalty'] = avg_cycle_length  # Fixed points (length=1) are better
    
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
    
    # Calculate composite score (lower is better)
    # Weight different factors
    composite_score = (
        0.3 * metrics['basin_penalty'] +           # Prefer fewer basins
        0.3 * metrics['complexity_penalty'] +      # Prefer simpler attractors
        0.2 * (1 - metrics['final_state_concordance']) +  # Prefer concordance with data
        0.2 * (1 - metrics['fixed_point_ratio'])   # Prefer fixed points
    )
    
    metrics['composite_score'] = composite_score
    
    return metrics


def load_rules_table(rules_file: str) -> pd.DataFrame:
    """Load rules table from TSV file"""
    df = pd.read_csv(rules_file, sep='\t', encoding='utf-8')
    required_cols = ['Gene', 'Position', 'Rule']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    return df


def load_binarized_matrix(matrix_file: str) -> pd.DataFrame:
    """Load binarized expression matrix"""
    # Load matrix and handle with or without index column
    df = pd.read_csv(matrix_file, sep='\t')
    
    # Check if first column looks like gene names or timepoints
    first_col = df.columns[0]
    
    # If first column doesn't look like a gene name (e.g., starts with 'G', 'AT', etc)
    # assume it's a time/sample identifier and set it as index
    if not any(df.columns[0].startswith(prefix) for prefix in ['G', 'AT', 'LOC', 'ENSG']):
        df = df.set_index(df.columns[0])
    
    return df


def calculate_intelligent_defaults(rules_by_gene: Dict[str, List], 
                                  target_time_minutes: int = 10,
                                  n_processes: int = 8) -> Dict[str, int]:
    """
    Calculate intelligent default values for top_n and max_combinations
    based on the problem size and desired execution time.
    
    Args:
        rules_by_gene: Dictionary with rules per gene
        target_time_minutes: Target execution time in minutes (default: 10)
        n_processes: Number of processes (default: 8, conservative)
        
    Returns:
        Dictionary with recommended 'top_n' and 'max_combinations'
    """
    # Get number of rules per gene
    rules_counts = [len(rules) for rules in rules_by_gene.values()]
    n_genes = len(rules_by_gene)
    
    # ===== LÍMITES ABSOLUTOS CONSERVADORES =====
    MAX_COMBINATIONS_ABSOLUTE = 50_000  # Límite máximo absoluto
    MAX_RULES_PER_GENE_ABSOLUTE = 30    # Máximo de reglas por gen
    TARGET_COMBINATIONS_FOR_SAMPLING = 20_000  # Target preferido para muestreo
    
    # Calculate total possible combinations (with overflow protection)
    try:
        total_combinations = np.prod(rules_counts)
        # Check for overflow (negative numbers or unreasonably large)
        if total_combinations < 0 or total_combinations > 1e18:
            total_combinations = float('inf')
    except (OverflowError, ValueError):
        total_combinations = float('inf')
    
    # Estimate throughput (combinations per second)
    # Conservative estimate: 3 combinations/second/core
    throughput_per_core = 3
    total_throughput = throughput_per_core * n_processes
    
    # Calculate how many combinations we can evaluate in target time
    max_feasible = int(total_throughput * target_time_minutes * 60)
    
    # Apply absolute limit (no more than 50K regardless of time)
    max_feasible = min(max_feasible, MAX_COMBINATIONS_ABSOLUTE)
    
    # Strategy 1: If total combinations is small, evaluate all
    if total_combinations <= max_feasible:
        return {
            'top_n': None,  # Use all rules with max score
            'max_combinations': None,  # Evaluate all
            'reason': 'small_problem',
            'estimated_time_minutes': total_combinations / total_throughput / 60 if total_combinations != float('inf') else target_time_minutes
        }
    
    # Strategy 2: If limiting combinations is enough (within 10x of feasible)
    if total_combinations <= max_feasible * 10:
        # We can sample a reasonable fraction
        # Prefer our target of 20K if possible
        recommended_combos = min(TARGET_COMBINATIONS_FOR_SAMPLING, total_combinations)
        return {
            'top_n': None,  # Use all rules with max score
            'max_combinations': int(recommended_combos),
            'reason': 'limit_combinations_only',
            'estimated_time_minutes': target_time_minutes
        }
    
    # Strategy 3: Need to reduce rules per gene AND limit combinations
    # This is the most common case for large problems
    
    # Calculate recommended top_n to get ~20K combinations
    # Formula: top_n^n_genes ≈ TARGET_COMBINATIONS_FOR_SAMPLING
    if n_genes > 0:
        recommended_top_n = int(TARGET_COMBINATIONS_FOR_SAMPLING ** (1 / n_genes))
    else:
        recommended_top_n = 5
    
    # Apply constraints on top_n
    recommended_top_n = max(3, recommended_top_n)  # Minimum 3 rules per gene
    recommended_top_n = min(MAX_RULES_PER_GENE_ABSOLUTE, recommended_top_n)  # Maximum 30
    
    # Calculate how many combinations this would generate
    estimated_combos = recommended_top_n ** n_genes
    
    # If still too many, apply max_combinations limit
    if estimated_combos > TARGET_COMBINATIONS_FOR_SAMPLING:
        recommended_max_combos = TARGET_COMBINATIONS_FOR_SAMPLING
    else:
        recommended_max_combos = None  # Can evaluate all with this top_n
    
    return {
        'top_n': recommended_top_n,
        'max_combinations': recommended_max_combos,
        'reason': 'reduce_both',
        'estimated_time_minutes': target_time_minutes
    }


def get_top_rules_per_gene(df: pd.DataFrame, 
                           top_n: int = None,
                           score_column: str = 'Score') -> Dict[str, List[Tuple[int, str]]]:
    """
    Get top rules for each gene based on score
    
    By default (top_n=None), returns ALL rules with the maximum score (tied for best).
    If top_n is specified, returns only the first top_n rules.
    
    Args:
        df: Rules dataframe
        top_n: Number of top rules to consider per gene. 
               If None (default), returns all rules with maximum score.
        score_column: Column name to use for ranking
        
    Returns:
        Dictionary mapping gene names to list of (position, rule) tuples
    """
    rules_by_gene = {}
    
    for gene in sorted(df['Gene'].unique()):
        gene_rules = df[df['Gene'] == gene].copy()
        
        # Sort by score (descending)
        if score_column in gene_rules.columns:
            gene_rules = gene_rules.sort_values(score_column, ascending=False)
            
            if top_n is None:
                # Get all rules with the maximum score (tied for best)
                max_score = gene_rules[score_column].max()
                top_rules = gene_rules[gene_rules[score_column] == max_score]
            else:
                # Get top_n rules
                top_rules = gene_rules.head(top_n)
        else:
            # Fall back to position sorting if no score column
            gene_rules = gene_rules.sort_values('Position', ascending=True)
            if top_n is None:
                top_rules = gene_rules.head(10)  # Default to 10 if no score and no top_n
            else:
                top_rules = gene_rules.head(top_n)
        
        rules_list = [(int(row['Position']), str(row['Rule'])) 
                     for _, row in top_rules.iterrows()]
        
        rules_by_gene[gene] = rules_list
    
    return rules_by_gene


def evaluate_single_combination(combo_indices: List[int],
                                genes: List[str],
                                rules_by_gene: Dict[str, List[Tuple[int, str]]],
                                binarized_matrix: pd.DataFrame = None,
                                max_iterations: int = 1000) -> Tuple[List[int], Dict[str, float]]:
    """
    Evaluate a single combination of rules
    
    Args:
        combo_indices: List of rule indices for each gene
        genes: List of gene names
        rules_by_gene: Dictionary of available rules per gene
        binarized_matrix: Optional binarized expression matrix
        max_iterations: Maximum iterations for attractor search
        
    Returns:
        Tuple of (combo_indices, metrics_dict)
    """
    # Build gene_rules dictionary for this combination
    gene_rules = {}
    for i, gene in enumerate(genes):
        rule_idx = combo_indices[i]
        position, rule = rules_by_gene[gene][rule_idx]
        gene_rules[gene] = rule
    
    # Find attractors
    attractors, basins = find_attractors_for_ruleset(gene_rules, genes, max_iterations)
    
    # Calculate attractor-based metrics
    metrics = calculate_attractor_metrics(attractors, basins, genes, binarized_matrix)
    
    # Calculate parsimony-based metrics (NEW)
    parsimony_metrics = calculate_parsimony_metrics(gene_rules)
    metrics.update(parsimony_metrics)
    
    # Add position information to metrics
    metrics['positions'] = [rules_by_gene[gene][combo_indices[i]][0] 
                           for i, gene in enumerate(genes)]
    
    return combo_indices, metrics


def evaluate_rule_combinations(rules_file: str,
                               output_file: str = None,
                               binarized_matrix_file: str = None,
                               top_n: int = None,
                               max_combinations: int = None,
                               max_iterations: int = 1000,
                               n_processes: int = 8,
                               score_column: str = 'Score',
                               verbose: bool = True):
    """
    Main function to evaluate rule combinations
    
    Args:
        rules_file: Path to rules TSV file
        output_file: Path to output results file (default: same dir as rules_file)
        binarized_matrix_file: Optional path to binarized expression matrix
        top_n: Number of top rules to consider per gene. 
               If None (default), uses all rules with maximum score (tied for best).
        max_combinations: Maximum number of combinations to evaluate.
                         If None (default), evaluates ALL possible combinations.
        max_iterations: Maximum iterations for attractor search
        n_processes: Number of processes for parallelization (default: 8)
        score_column: Column to use for rule ranking
        verbose: Print progress information
    """
    if verbose:
        print("="*80)
        print("BOOLEAN RULES EVALUATION WITH INTEGRATED SCORE")
        print("="*80)
    
    # Load data
    if verbose:
        print(f"\n1. Loading rules from: {rules_file}")
    df = load_rules_table(rules_file)
    
    # Generate output file path if not provided
    if output_file is None:
        rules_dir = os.path.dirname(rules_file)
        output_file = os.path.join(rules_dir, "evaluation_results.tsv")
        if verbose:
            print(f"   Output will be saved to: {output_file}")
    
    binarized_matrix = None
    if binarized_matrix_file:
        if verbose:
            print(f"2. Loading binarized matrix from: {binarized_matrix_file}")
        binarized_matrix = load_binarized_matrix(binarized_matrix_file)
    
    # Get top rules per gene
    # If both top_n and max_combinations are None, calculate intelligent defaults
    if top_n is None and max_combinations is None:
        if verbose:
            print(f"\n3. Calculating intelligent defaults based on problem size...")
        
        # First pass: get rules with max score to assess problem size
        temp_rules = get_top_rules_per_gene(df, top_n=None, score_column=score_column)
        
        # Calculate intelligent defaults
        defaults = calculate_intelligent_defaults(
            temp_rules, 
            target_time_minutes=10,
            n_processes=n_processes
        )
        
        # Apply recommended defaults
        top_n = defaults['top_n']
        max_combinations = defaults['max_combinations']
        
        if verbose:
            print(f"   Problem size analysis:")
            print(f"   - Strategy: {defaults['reason']}")
            if top_n is not None:
                print(f"   - Recommended top_n: {top_n} rules per gene")
            else:
                print(f"   - Using ALL rules with maximum score per gene")
            if max_combinations is not None:
                print(f"   - Recommended max_combinations: {max_combinations:,}")
            else:
                print(f"   - Will evaluate ALL possible combinations")
            print(f"   - Estimated time: ~{defaults['estimated_time_minutes']:.1f} minutes")
    
    if verbose:
        if top_n is None:
            print(f"\n4. Extracting rules with MAXIMUM score per gene (all tied for best)...")
        else:
            print(f"\n4. Extracting top {top_n} rules per gene...")
    
    rules_by_gene = get_top_rules_per_gene(df, top_n, score_column)
    genes = sorted(rules_by_gene.keys())
    
    if verbose:
        print(f"   Found {len(genes)} genes")
        for gene in genes:
            n_rules = len(rules_by_gene[gene])
            if n_rules > 1:
                print(f"   - {gene}: {n_rules} candidate rules")
            else:
                print(f"   - {gene}: {n_rules} rule")
    
    # Calculate total possible combinations
    total_combinations = np.prod([len(rules_by_gene[gene]) for gene in genes])
    
    if verbose:
        print(f"\n5. Total possible combinations: {total_combinations:,}")
    
    # Decide how many combinations to evaluate
    if max_combinations is None:
        # Evaluate ALL combinations (default behavior)
        n_combos_to_eval = total_combinations
        if verbose:
            print(f"   Evaluating ALL {total_combinations:,} combinations")
    else:
        # User specified a limit
        n_combos_to_eval = min(max_combinations, total_combinations)
        if total_combinations > max_combinations:
            if verbose:
                print(f"   Limiting evaluation to {max_combinations:,} combinations (random sampling)")
        else:
            if verbose:
                print(f"   Evaluating all {total_combinations:,} combinations (less than max_combinations)")
    
    # Generate combinations to evaluate
    if n_combos_to_eval == total_combinations:
        # Evaluate all combinations
        if verbose:
            print(f"   Generating all possible combinations...")
        combinations = list(itertools.product(*[range(len(rules_by_gene[gene])) 
                                                for gene in genes]))
    else:
        # Random sampling of combinations
        if verbose:
            print(f"   Random sampling {n_combos_to_eval:,} combinations...")
        np.random.seed(42)
        combinations = []
        for _ in range(n_combos_to_eval):
            combo = [np.random.randint(0, len(rules_by_gene[gene])) 
                    for gene in genes]
            combinations.append(tuple(combo))
        combinations = list(set(combinations))  # Remove duplicates
    
    if verbose:
        print(f"   Ready to evaluate {len(combinations):,} unique combinations")
    
    # Set up parallelization
    if n_processes is None:
        n_processes = mp.cpu_count()
    
    if verbose:
        print(f"\n6. Starting evaluation with {n_processes} processes...")
    
    # Create partial function with fixed parameters
    eval_func = partial(evaluate_single_combination,
                       genes=genes,
                       rules_by_gene=rules_by_gene,
                       binarized_matrix=binarized_matrix,
                       max_iterations=max_iterations)
    
    # Evaluate combinations in parallel
    results = []
    combinations_list = list(combinations)  # Convert to list for iteration
    
    if n_processes > 1:
        with mp.Pool(processes=n_processes) as pool:
            for i, result in enumerate(pool.imap_unordered(eval_func, combinations_list), 1):
                results.append(result)
                if verbose and i % 100 == 0:
                    print(f"   Progress: {i}/{len(combinations_list)} combinations evaluated", 
                          end='\r')
    else:
        for i, combo in enumerate(combinations_list, 1):
            result = eval_func(combo)
            results.append(result)
            if verbose and i % 100 == 0:
                print(f"   Progress: {i}/{len(combinations_list)} combinations evaluated", 
                      end='\r')
    
    if verbose:
        print(f"\n   Completed: {len(results)} combinations evaluated")
    
    # Process results
    if verbose:
        print(f"\n7. Processing and ranking results...")
    
    results_data = []
    for combo_indices, metrics in results:
        row = {
            'combination_id': '_'.join(map(str, combo_indices)),
            'composite_score': metrics['composite_score'],
            # tiebreak_score will be calculated after creating DataFrame
            'n_basins': metrics['n_basins'],
            'avg_cycle_length': metrics['avg_cycle_length'],
            'basin_entropy': metrics['basin_entropy'],
            'final_state_concordance': metrics['final_state_concordance'],
            'fixed_point_ratio': metrics['fixed_point_ratio'],
            # Parsimony metrics
            'total_literals': metrics['total_literals'],
            'total_nots': metrics['total_nots'],
            'avg_k': metrics['avg_k'],
            'std_k': metrics['std_k'],
            'k_distance_from_2': metrics['k_distance_from_2']
        }
        
        # Add individual gene positions and rules
        for i, gene in enumerate(genes):
            position = metrics['positions'][i]
            rule = rules_by_gene[gene][combo_indices[i]][1]
            row[f'{gene}_position'] = position
            row[f'{gene}_rule'] = rule
        
        results_data.append(row)
    
    # Create results dataframe
    results_df = pd.DataFrame(results_data)
    
    # Calculate FINAL INTEGRATED SCORE that combines ALL metrics
    # Lower is better for all components
    # This replaces the two-tier system (composite + tiebreak) with a single unified score
    
    # Normalize each metric to [0, 1] range for fair comparison
    def safe_normalize(series):
        """Normalize series to 0-1, handling edge cases"""
        max_val = series.max()
        min_val = series.min()
        if max_val == min_val:
            return pd.Series([0.0] * len(series), index=series.index)
        return (series - min_val) / (max_val - min_val)
    
    # Attractor-based penalties (from composite_score)
    norm_n_basins = safe_normalize(results_df['n_basins'])
    norm_cycle_length = safe_normalize(results_df['avg_cycle_length'])
    norm_basin_entropy = safe_normalize(results_df['basin_entropy'])
    # For these, higher is better, so invert
    norm_concordance = 1 - results_df['final_state_concordance']  # Already 0-1
    norm_fixed_ratio = 1 - results_df['fixed_point_ratio']  # Already 0-1
    
    # Parsimony-based penalties
    norm_literals = safe_normalize(results_df['total_literals'])
    norm_nots = safe_normalize(results_df['total_nots'])
    norm_k_distance = safe_normalize(results_df['k_distance_from_2'])
    
    # Calculate FINAL UNIFIED SCORE
    # Weights are carefully chosen to balance attractor dynamics vs parsimony
    results_df['final_score'] = (
        # Attractor metrics (60% total weight)
        0.20 * norm_n_basins +              # Fewer basins = better
        0.15 * norm_cycle_length +          # Shorter cycles = better  
        0.10 * norm_basin_entropy +         # Dominant basin = better
        0.10 * norm_concordance +           # Match data = better
        0.05 * norm_fixed_ratio +           # Fixed points = better
        
        # Parsimony metrics (40% total weight)
        0.20 * norm_k_distance +            # K≈2 = better (most important parsimony)
        0.12 * norm_literals +              # Fewer literals = simpler
        0.08 * norm_nots                    # Fewer NOTs = less repression
    )
    
    # Keep old scores for reference/transparency
    results_df['composite_score_old'] = results_df['composite_score']
    
    # Sort by final_score only (single unified ranking)
    results_df = results_df.sort_values('final_score', ascending=True)
    
    # Reorder columns to put final_score first
    cols = list(results_df.columns)
    cols.remove('final_score')
    cols.insert(1, 'final_score')  # After combination_id
    results_df = results_df[cols]
    
    # Save results
    if verbose:
        print(f"\n8. Saving results to: {output_file}")
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Generate rules_by_gene_evaluated.tsv with winning rules
    if verbose:
        print(f"\n9. Generating rules_by_gene_evaluated.tsv with winning rules...")
    
    # Get the best combination(s)
    best_final_score = results_df['final_score'].min()
    best_combinations = results_df[results_df['final_score'] == best_final_score]
    
    # Create a new dataframe with the winning rules
    winning_rules_data = []
    
    for gene in genes:
        # Get the position(s) and rule(s) from best combination(s)
        positions_in_best = set()
        rules_in_best = set()
        
        for _, row in best_combinations.iterrows():
            pos = row[f'{gene}_position']
            rule = row[f'{gene}_rule']
            positions_in_best.add(pos)
            rules_in_best.add((pos, rule))
        
        # Look up original scores for these positions in the original dataframe
        for pos, rule in rules_in_best:
            # Find this rule in the original dataframe
            original_row = df[(df['Gene'] == gene) & 
                            (df['Position'] == pos) & 
                            (df['Rule'] == rule)]
            
            if not original_row.empty:
                # Get the first match (should be unique)
                orig = original_row.iloc[0]
                
                # Create entry for winning rules file maintaining original column order
                # Original order: Gene, Position, Rule, Correct, N_Regulators, MSE, Score
                entry = {
                    'Gene': gene,
                    'Position': int(pos),
                    'Rule': rule,
                }
                
                # Add columns in the original order if they exist
                for col in df.columns:
                    if col not in ['Gene', 'Position', 'Rule']:
                        entry[col] = orig[col]
                
                winning_rules_data.append(entry)
    
    # Create dataframe and sort
    winning_rules_df = pd.DataFrame(winning_rules_data)
    winning_rules_df = winning_rules_df.sort_values(['Gene', 'Position'])
    
    # Reorder columns to match original file format
    # Standard order: Gene, Position, Rule, Correct, N_Regulators, MSE, Score
    # But preserve whatever columns exist in the original file
    original_columns = ['Gene', 'Position', 'Rule']
    
    # Add remaining columns in the order they appear in the original dataframe
    for col in df.columns:
        if col not in original_columns and col in winning_rules_df.columns:
            original_columns.append(col)
    
    # Reorder dataframe columns
    winning_rules_df = winning_rules_df[original_columns]
    
    # Save to file in same directory as output_file
    output_dir = os.path.dirname(output_file)
    winning_rules_file = os.path.join(output_dir, "rules_by_gene_evaluated.tsv")
    winning_rules_df.to_csv(winning_rules_file, sep='\t', index=False)
    
    if verbose:
        print(f"   Saved {len(winning_rules_df)} winning rules to: {winning_rules_file}")
        n_best_combos = len(best_combinations)
        if n_best_combos > 1:
            print(f"   Note: {n_best_combos} combinations tied for best final_score")
            print(f"         The file contains all rules from these tied combinations")
    
    # Print summary
    if verbose:
        print("\n" + "="*80)
        print("EVALUATION SUMMARY")
        print("="*80)
        print(f"Best combination (sorted by final integrated score):")
        best = results_df.iloc[0]
        print(f"  Final score: {best['final_score']:.4f} (LOWER IS BETTER)")
        print(f"  [Combines ALL metrics: attractor dynamics + parsimony]")
        print(f"\n  Attractor metrics:")
        print(f"    Number of basins: {int(best['n_basins'])}")
        print(f"    Avg cycle length: {best['avg_cycle_length']:.2f}")
        print(f"    Basin entropy: {best['basin_entropy']:.4f}")
        print(f"    Fixed point ratio: {best['fixed_point_ratio']:.2f}")
        if binarized_matrix is not None:
            print(f"    Final state concordance: {best['final_state_concordance']:.4f}")
        print(f"\n  Parsimony metrics:")
        print(f"    Total literals: {int(best['total_literals'])} (fewer = simpler)")
        print(f"    Total NOT operators: {int(best['total_nots'])} (fewer = less repression)")
        print(f"    Avg K (connectivity): {best['avg_k']:.2f} (optimal ≈ 2.0)")
        print(f"    K distance from 2.0: {best['k_distance_from_2']:.3f} (closer to 0 = better)")
        
        # Show comparison with old composite_score
        if 'composite_score_old' in best:
            print(f"\n  Reference (old composite_score): {best['composite_score_old']:.4f}")
        
        # Show score range
        print(f"\n  Score range in dataset:")
        print(f"    Best (this): {results_df['final_score'].min():.4f}")
        print(f"    Worst: {results_df['final_score'].max():.4f}")
        
        print(f"\nAll {len(results_df)} combinations saved to output file")
        print(f"Sorted by final_score (single unified metric)")
        print("="*80)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Evaluate Boolean rule combinations using attractor metrics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: Evaluate ALL rules with max score, ALL combinations
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv -m binarized_matrix.tsv -v
  
  # Limit to top 5 rules per gene (instead of all tied for best)
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --top-n 5 -v
  
  # Limit to 500 combinations max (instead of all)
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --max-combos 500 -v
  
  # Both limits combined
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --top-n 3 --max-combos 200 -v
  
  # Custom output directory
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv -o /path/to/results.tsv -v
  
  # Use custom score column for ranking
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --score-col MSE -v

Metrics calculated:
  - final_score: UNIFIED score combining ALL metrics (LOWER IS BETTER)
      Weighted combination of normalized metrics:
      Attractor dynamics (60% weight):
        - 20% Number of basins (fewer = better)
        - 15% Cycle length (shorter = better)
        - 10% Basin entropy (dominant basin = better)
        - 10% Data concordance (match = better)
        - 5%  Fixed point ratio (more = better)
      Parsimony (40% weight):
        - 20% K distance from 2.0 (K≈2 optimal, Kauffman)
        - 12% Total literals (fewer = simpler)
        - 8%  Total NOTs (fewer = less repression)
  
  Individual metrics (for reference):
  - n_basins: Number of attractor basins
  - avg_cycle_length: Average attractor cycle length
  - basin_entropy: Basin size distribution uniformity
  - final_state_concordance: Match with last timepoint
  - fixed_point_ratio: Proportion of fixed point attractors
  - total_literals: Total gene mentions (parsimony)
  - total_nots: Total NOT operators (repression level)
  - avg_k: Average regulators per gene (connectivity)
  - k_distance_from_2: Distance from optimal K=2
  - composite_score_old: Original composite score (reference)

Output format:
  Two TSV files are generated:
  
  1. evaluation_results.tsv:
     - Sorted by final_score (ascending)
     - First row is the overall best combination
     - All metrics normalized and integrated into single ranking
     - Contains all evaluated combinations
  
  2. rules_by_gene_evaluated.tsv:
     - Contains the winning rules (one per gene, or more if tied)
     - Same format as input rules_by_gene.tsv
     - Can be used directly with BNI3_Attractors.py
     - If multiple combinations tied for best, includes all their rules
        """
    )
    
    # Required arguments
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input', type=str, required=True,
                         help='Input rules TSV file')
    
    # Optional arguments
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-o', '--output', type=str, default=None,
                         help='Output results TSV file (default: evaluation_results.tsv in same directory as input)')
    optional.add_argument('-m', '--matrix', type=str, default=None,
                         help='Binarized expression matrix TSV file (for concordance metric)')
    optional.add_argument('--top-n', type=int, default=None,
                         help='Number of top rules to consider per gene. '
                              'Default: None (AUTOMATIC - calculates intelligent default based on problem size, typically 3-30)')
    optional.add_argument('--max-combos', type=int, default=None,
                         help='Maximum number of combinations to evaluate. '
                              'Default: None (AUTOMATIC - calculates intelligent default, typically ~20,000)')
    optional.add_argument('--max-iter', type=int, default=1000,
                         help='Maximum iterations for attractor search (default: 1000)')
    optional.add_argument('-n', '--processes', type=int, default=8,
                         help='Number of parallel processes (default: 8, conservative and optimal)')
    optional.add_argument('--score-col', type=str, default='Score',
                         help='Column name for rule ranking (default: Score)')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Print detailed progress information')
    
    args = parser.parse_args()
    
    # Run evaluation
    try:
        evaluate_rule_combinations(
            rules_file=args.input,
            output_file=args.output,
            binarized_matrix_file=args.matrix,
            top_n=args.top_n,
            max_combinations=args.max_combos,
            max_iterations=args.max_iter,
            n_processes=args.processes,
            score_column=args.score_col,
            verbose=args.verbose
        )
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

print(f"   WARNING: Combination overflow detected. Total would be > 1e18")