#!/usr/bin/env python3
"""
Boolean Network Attractor Analysis Tool
Analyzes Boolean network attractors and basins of attraction from logical rules.
Part of the Boolean Network Inference (BNI) pipeline.
"""

import argparse
import os
import sys
import pandas as pd
import itertools
import textwrap
import re
from collections import defaultdict
import multiprocessing as mp
from pathlib import Path


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def read_rules_table(rules_file, verbose):
    """
    Read the rules table from TSV file
    
    Args:
        rules_file (str): Path to the rules table file
        verbose (bool): Enable verbose output
        
    Returns:
        pd.DataFrame: Rules dataframe or None if error
    """
    try:
        log_message(f"Reading rules table from: {rules_file}", verbose)
        
        if rules_file.endswith('.tsv'):
            df = pd.read_csv(rules_file, sep='\t', encoding='utf-8')
        else:
            df = pd.read_csv(rules_file, encoding='utf-8')
        
        log_message(f"Table loaded successfully with {len(df)} rules", verbose)
        
        # Validate required columns
        required_cols = ['Gene', 'Rule']
        # Check if Position column exists (optional)
        if 'Position' not in df.columns:
            # If no Position column, add it with value 1 for all rows
            df['Position'] = 1
            log_message("No 'Position' column found, using Position=1 for all rules", verbose)
        
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
            
        return df
        
    except Exception as e:
        raise ValueError(f"Error reading rules file: {str(e)}")


def select_rules_by_criteria(df, criteria, custom_positions, verbose):
    """
    Select rules based on Position values and criteria
    
    Args:
        df (pd.DataFrame): Rules dataframe
        criteria (str): Selection criteria ('top_position' or 'custom_positions')
        custom_positions (str): Comma-separated list of positions (only used with 'custom_positions')
        verbose (bool): Enable verbose output
        
    Returns:
        dict: Dictionary with selected rules per gene
    """
    log_message(f"Selecting rules using criteria: {criteria}", verbose)
    
    selected_rules = {}
    genes = sorted(df['Gene'].unique())
    
    # Parse custom positions if provided
    position_list = None
    if criteria == "custom_positions":
        if not custom_positions:
            raise ValueError("custom_positions criteria requires --positions parameter")
        
        try:
            position_list = [int(x.strip()) for x in custom_positions.split(',')]
            log_message(f"Custom positions: {position_list}", verbose)
            
            if len(position_list) != len(genes):
                raise ValueError(f"Number of positions ({len(position_list)}) must match number of genes ({len(genes)}). Genes found: {', '.join(genes)}")
                
        except ValueError as e:
            raise ValueError(f"Invalid positions format. Use comma-separated integers: {e}")
    
    for i, gene in enumerate(genes):
        gene_rules = df[df['Gene'] == gene].copy()
        
        # Select rule based ONLY on Position values
        if criteria == "top_position":
            # Select rule with Position = 1
            top_rules = gene_rules[gene_rules['Position'] == 1]
            if top_rules.empty:
                # If no Position=1, take the minimum position available
                selected_rule = gene_rules.loc[gene_rules['Position'].idxmin()]
                log_message(f"Warning: No Position=1 found for {gene}, using Position={selected_rule['Position']}", verbose)
            else:
                selected_rule = top_rules.iloc[0]  # Take first if multiple Position=1
                
        elif criteria == "custom_positions":
            # Select rule with specified position
            target_position = position_list[i]
            target_rules = gene_rules[gene_rules['Position'] == target_position]
            
            if target_rules.empty:
                # If target position not found, take the minimum position available
                selected_rule = gene_rules.loc[gene_rules['Position'].idxmin()]
                log_message(f"Warning: Position={target_position} not found for {gene}, using Position={selected_rule['Position']}", verbose)
            else:
                selected_rule = target_rules.iloc[0]  # Take first if multiple with same position
        
        # Store only the rule string - convert pandas types to native Python types
        selected_rules[gene] = str(selected_rule['Rule'])
        log_message(f"{gene}: Position {int(selected_rule['Position'])} - {selected_rules[gene]}", verbose)
    
    log_message(f"Selected {len(selected_rules)} rules total", verbose)
    return selected_rules


def is_constant_rule(rule):
    """
    Check if a rule is a boolean constant.
    
    Args:
        rule (str): Boolean rule string
        
    Returns:
        tuple: (is_constant, value) where value is True/False if is_constant is True
    """
    rule_stripped = rule.strip().upper()
    
    if rule_stripped in ['1', 'TRUE']:
        return True, True
    elif rule_stripped in ['0', 'FALSE']:
        return True, False
    else:
        return False, None


def extract_genes_from_rules(gene_rules):
    """
    Extract all unique genes that appear in the rules.
    Excludes boolean constants (0, 1, TRUE, FALSE).
    
    Args:
        gene_rules (dict): Dictionary of gene names to rules
        
    Returns:
        list: Sorted list of all unique gene names
    """
    all_genes = set(gene_rules.keys())
    
    # Boolean constants to exclude
    bool_constants = {'0', '1', 'TRUE', 'FALSE', 'True', 'False', 'true', 'false'}
    
    # Extract genes mentioned in rules (right side)
    for rule in gene_rules.values():
        # Skip if rule is a boolean constant
        rule_stripped = rule.strip()
        if rule_stripped in bool_constants:
            continue
            
        # Parse rule to find mentioned genes
        tokens = rule.replace('&', ' ').replace('|', ' ').replace('~', ' ').replace('!', ' ').replace('(', ' ').replace(')', ' ').split()
        for token in tokens:
            if token and token not in ['AND', 'OR', 'NOT'] and token not in bool_constants:
                all_genes.add(token)
    
    return sorted(list(all_genes))


def parse_mutations(mutations_str, genes):
    """
    Parse mutations string and validate gene names.
    
    Args:
        mutations_str (str): Mutations in format "GENE1:1,GENE2:0"
        genes (list): List of valid gene names
        
    Returns:
        dict: Dictionary mapping gene names to fixed values (0 or 1)
    """
    mutations = {}
    if not mutations_str:
        return mutations
    
    try:
        for mutation in mutations_str.split(','):
            gene, value = mutation.strip().split(':')
            gene = gene.strip()
            value = int(value.strip())
            
            if gene not in genes:
                raise ValueError(f"Gene '{gene}' not found in network. Available genes: {', '.join(genes)}")
            
            if value not in [0, 1]:
                raise ValueError(f"Mutation value must be 0 or 1, got {value} for gene {gene}")
            
            mutations[gene] = value
            
    except ValueError as e:
        if "not enough values to unpack" in str(e):
            raise ValueError(f"Invalid mutations format. Use 'GENE1:1,GENE2:0'. Got: {mutations_str}")
        else:
            raise e
    
    return mutations


def generate_mutation_suffix(mutations):
    """
    Generate filename suffix based on mutations.
    
    Args:
        mutations (dict): Dictionary of gene mutations
        
    Returns:
        str: Suffix for filenames (e.g., "_ABF3_1_MLP7_0")
    """
    if not mutations:
        return ""
    
    mutation_parts = []
    for gene, value in sorted(mutations.items()):
        mutation_parts.append(f"{gene}_{value}")
    
    return "_" + "_".join(mutation_parts)


def evaluate_rule(rule, gene_state):
    """
    Evaluate a Boolean rule given a gene state.
    Handles constant rules (0, 1, TRUE, FALSE).
    
    Args:
        rule (str): Boolean rule string
        gene_state (dict): Dictionary mapping gene names to boolean values
        
    Returns:
        bool: Result of rule evaluation
    """
    # Check if rule is a constant
    is_const, const_value = is_constant_rule(rule)
    if is_const:
        return const_value
    
    # Replace operators - handle both ! and ~ for NOT
    rule_eval = rule.replace('&', ' and ').replace('|', ' or ')
    rule_eval = rule_eval.replace('!', ' not ').replace('~', ' not ')
    
    # Replace gene names with their values
    for gene, value in gene_state.items():
        # Use word boundaries to avoid partial replacements
        pattern = r'\b' + re.escape(gene) + r'\b'
        rule_eval = re.sub(pattern, str(value), rule_eval)
    
    try:
        result = eval(rule_eval)
        return bool(result)
    except Exception as e:
        # If evaluation error, return False by default
        print(f"Warning: Error evaluating rule '{rule}': {e}")
        return False


def calculate_next_state(current_state, gene_rules, genes, mutations=None):
    """
    Calculate the next state of the network given the current state.
    
    Args:
        current_state (list): Current state as list of booleans
        gene_rules (dict): Dictionary of gene rules
        genes (list): List of gene names in order
        mutations (dict): Dictionary of gene mutations {gene: fixed_value}
        
    Returns:
        list: Next state as list of booleans
    """
    next_state = {}
    mutations = mutations or {}
    
    for gene in genes:
        # Check if gene is mutated (fixed value)
        if gene in mutations:
            next_state[gene] = bool(mutations[gene])
        elif gene in gene_rules:
            # Create state dictionary for evaluation
            state_dict = {genes[i]: current_state[i] for i in range(len(genes))}
            next_state[gene] = evaluate_rule(gene_rules[gene], state_dict)
        else:
            # If no rule for gene, maintain current state
            idx = genes.index(gene)
            next_state[gene] = current_state[idx]
    
    return [next_state[gene] for gene in genes]


def state_to_string(state):
    """
    Convert a state (list of booleans) to string for easy comparison.
    
    Args:
        state (list): List of boolean values
        
    Returns:
        str: Binary string representation
    """
    return ''.join(['1' if x else '0' for x in state])


def _process_state_chunk(chunk_states, gene_rules, genes, max_iterations, mutations=None):
    """
    Process a chunk of states and return mapping of states to attractors
    and found attractors.
    
    Args:
        chunk_states (list): Chunk of initial states to process
        gene_rules (dict): Gene rules dictionary
        genes (list): List of gene names
        max_iterations (int): Maximum iterations for trajectory simulation
        mutations (dict): Dictionary of gene mutations
        
    Returns:
        tuple: (state_to_cycle_mapping, found_cycles)
    """
    # Local mapping of states to their attractors (as cycles)
    local_state_to_cycle = {}
    # List of unique cycles found
    found_cycles = []
    
    for initial_state in chunk_states:
        initial_state = list(initial_state)
        initial_state_str = state_to_string(initial_state)
        
        # If already processed this state, skip
        if initial_state_str in local_state_to_cycle:
            continue
        
        # Simulate trajectory
        trajectory = []
        current_state = initial_state[:]
        visited_trajectory = set()
        
        for _ in range(max_iterations):
            state_str = state_to_string(current_state)
            
            # If we found a cycle
            if state_str in visited_trajectory:
                # Find the cycle
                cycle_idx = next(i for i, s in enumerate(trajectory) 
                               if state_to_string(s) == state_str)
                cycle = trajectory[cycle_idx:]
                
                # Convert cycle to canonical format (starting with smallest state)
                cycle_strings = [state_to_string(s) for s in cycle]
                min_idx = cycle_strings.index(min(cycle_strings))
                canonical_cycle = cycle[min_idx:] + cycle[:min_idx]
                canonical_cycle_str = tuple(state_to_string(s) for s in canonical_cycle)
                
                # Map entire trajectory to this cycle
                for traj_state in trajectory:
                    local_state_to_cycle[state_to_string(traj_state)] = canonical_cycle_str
                
                # Add cycle if new
                if canonical_cycle_str not in [tuple(state_to_string(s) for s in c) for c in found_cycles]:
                    found_cycles.append(canonical_cycle)
                
                break
            
            trajectory.append(current_state[:])
            visited_trajectory.add(state_str)
            current_state = calculate_next_state(current_state, gene_rules, genes, mutations)
    
    return local_state_to_cycle, found_cycles


def find_attractors_parallel(gene_rules, genes, max_iterations=1000, n_procs=None, verbose=False, mutations=None):
    """
    Find all attractors using parallel processing and correctly count basin sizes.
    
    Args:
        gene_rules (dict): Gene rules dictionary
        genes (list): List of gene names
        max_iterations (int): Maximum iterations for simulation
        n_procs (int): Number of processes (None for auto-detect)
        verbose (bool): Enable verbose output
        mutations (dict): Dictionary of gene mutations
        
    Returns:
        tuple: (attractors_list, basin_sizes_list)
    """
    n_genes = len(genes)
    n_mutated = len(mutations) if mutations else 0
    n_free_genes = n_genes - n_mutated

    # Generate only valid states (non-mutated genes)
    free_genes_indices = [i for i, gene in enumerate(genes) if gene not in (mutations or {})]
    if free_genes_indices:
        free_states = list(itertools.product([False, True], repeat=len(free_genes_indices)))
        
        # Expand to full states by adding mutated gene values
        all_states = []
        for free_state in free_states:
            full_state = [False] * n_genes
            # Set free genes
            for i, free_idx in enumerate(free_genes_indices):
                full_state[free_idx] = free_state[i]
            # Set mutated genes
            for gene, value in (mutations or {}).items():
                gene_idx = genes.index(gene)
                full_state[gene_idx] = bool(value)
            all_states.append(tuple(full_state))
    else:
        # All genes are mutated - only one possible state
        full_state = [bool(mutations[gene]) for gene in genes]
        all_states = [tuple(full_state)]

    total_states = len(all_states)

    log_message(f"Analyzing {total_states} possible states (2^{n_free_genes}) with mutations on {n_mutated} genes...", verbose)
    
    if n_procs is None:
        n_procs = mp.cpu_count()
    
    log_message(f"Using {n_procs} processes...", verbose)
    
    # Divide states into chunks for each process
    chunk_size = len(all_states) // n_procs
    chunks = []
    for i in range(n_procs):
        start = i * chunk_size
        if i == n_procs - 1:
            # Last chunk takes all remaining
            end = len(all_states)
        else:
            end = start + chunk_size
        chunks.append(all_states[start:end])
    
    # Process chunks in parallel
    with mp.Pool(processes=n_procs) as pool:
        results = pool.starmap(
            _process_state_chunk, 
            [(chunk, gene_rules, genes, max_iterations, mutations) for chunk in chunks]
        )
    
    log_message("Consolidating results from all processes...", verbose)
    
    # Consolidate results
    # First, identify all unique attractors
    unique_attractors = []
    attractor_strings_to_idx = {}  # Mapping from canonical representation to index
    
    for _, chunk_cycles in results:
        for cycle in chunk_cycles:
            cycle_str = tuple(state_to_string(s) for s in cycle)
            if cycle_str not in attractor_strings_to_idx:
                idx = len(unique_attractors)
                unique_attractors.append(cycle)
                attractor_strings_to_idx[cycle_str] = idx
    
    # Now count states in each basin
    basins = [0] * len(unique_attractors)
    processed_states = set()
    
    for local_state_to_cycle, _ in results:
        for state_str, cycle_str in local_state_to_cycle.items():
            if state_str not in processed_states:
                processed_states.add(state_str)
                if cycle_str in attractor_strings_to_idx:
                    idx = attractor_strings_to_idx[cycle_str]
                    basins[idx] += 1
    
    # Verify we processed all states
    log_message(f"Unique states processed: {len(processed_states)}/{total_states}", verbose)
    
    # Process any missing states (shouldn't happen)
    missing_states = []
    for state in all_states:
        if state_to_string(list(state)) not in processed_states:
            missing_states.append(state)
    
    if missing_states:
        log_message(f"Processing {len(missing_states)} missing states...", verbose)
        extra_state_to_cycle, extra_cycles = _process_state_chunk(
            missing_states, gene_rules, genes, max_iterations, mutations
        )
        
        for cycle in extra_cycles:
            cycle_str = tuple(state_to_string(s) for s in cycle)
            if cycle_str not in attractor_strings_to_idx:
                idx = len(unique_attractors)
                unique_attractors.append(cycle)
                attractor_strings_to_idx[cycle_str] = idx
                basins.append(0)
        
        for state_str, cycle_str in extra_state_to_cycle.items():
            if state_str not in processed_states:
                processed_states.add(state_str)
                if cycle_str in attractor_strings_to_idx:
                    idx = attractor_strings_to_idx[cycle_str]
                    basins[idx] += 1
    
    log_message(f"Total states mapped: {len(processed_states)}", verbose)
    log_message(f"Sum of all basins: {sum(basins)}", verbose)
    
    return unique_attractors, basins


def save_attractors_tsv(attractors, genes, basins, output_path, filename, verbose=False):
    """
    Save found attractors to TSV format including basin size information.
    
    Args:
        attractors (list): List of attractor cycles
        genes (list): List of gene names
        basins (list): Basin sizes for each attractor
        output_path (str): Output directory path
        filename (str): Output filename
        verbose (bool): Enable verbose output
        
    Returns:
        pd.DataFrame: DataFrame with attractor data
    """
    attractor_data = []
    total_basin_size = sum(basins) if basins else 0
    
    for i, attractor in enumerate(attractors, 1):
        att_type = "fixed_point" if len(attractor) == 1 else "cycle"
        length = len(attractor)
        basin_size = basins[i-1] if basins else 0
        basin_percentage = (basin_size / total_basin_size) * 100 if total_basin_size > 0 else 0
        
        for j, state in enumerate(attractor):
            row = {
                'attractor_id': i,
                'type': att_type,
                'cycle_length': length,
                'step_in_cycle': j + 1,
                'basin_size': basin_size,
                'basin_percentage': round(basin_percentage, 2)
            }
            
            # Add state of each gene
            for k, gene in enumerate(genes):
                row[gene] = int(state[k])
            
            # Add complete binary representation
            row['binary_state'] = state_to_string(state)
            
            attractor_data.append(row)
    
    # Create DataFrame and save
    df_attractors = pd.DataFrame(attractor_data)
    full_path = os.path.join(output_path, filename)
    df_attractors.to_csv(full_path, sep='\t', index=False)
    
    log_message(f"Attractors saved to: {filename}", verbose)
    return df_attractors


def print_attractors_summary(attractors, genes, basins=None, verbose=False, mutations=None):
    """
    Print detailed information about attractors and their basins.
    
    Args:
        attractors (list): List of attractor cycles
        genes (list): List of gene names
        basins (list): List of basin sizes
        verbose (bool): Enable verbose output
        mutations (dict): Dictionary of gene mutations
    """
    print(f"\n{'='*60}")
    print(f"Found {len(attractors)} attractor(s)")
    print(f"{'='*60}\n")
    
    if basins:
        total_states = sum(basins)
        mutations = mutations or {}
        n_mutated = len(mutations)
        n_free_genes = len(genes) - n_mutated
        expected = 2**n_free_genes
        
        print(f"Total states in system: {total_states}")
        if n_mutated > 0:
            print(f"Expected states: 2^{n_free_genes} = {expected} (with {n_mutated} gene(s) mutated)")
            print(f"Mutated genes: {', '.join([f'{gene}={value}' for gene, value in mutations.items()])}")
        else:
            print(f"Expected states: 2^{len(genes)} = {expected}")
        
        if total_states == expected:
            print(f"✓ Sum matches expected total\n")
        else:
            print(f"✗ WARNING: Missing {expected - total_states} states\n")
    
    for i, attractor in enumerate(attractors, 1):
        print(f"{'─'*50}")
        if basins:
            percentage = (basins[i-1] / sum(basins)) * 100 if sum(basins) > 0 else 0
            print(f"Attractor {i}:")
            print(f"  Basin size: {basins[i-1]} states ({percentage:.2f}%)")
        else:
            print(f"Attractor {i}:")
        
        if len(attractor) == 1:
            print(f"  Type: Fixed point")
        else:
            print(f"  Type: Cycle of length {len(attractor)}")
        
        # Show all states of the attractor with active/inactive genes
        for j, state in enumerate(attractor):
            state_dict = {genes[k]: int(state[k]) for k in range(len(genes))}
            active_genes = [g for g, v in state_dict.items() if v == 1]
            inactive_genes = [g for g, v in state_dict.items() if v == 0]
            
            print(f"    State {j+1}:")
            print(f"      Active genes ({len(active_genes)}): {', '.join(active_genes) if active_genes else 'none'}")
            print(f"      Inactive genes ({len(inactive_genes)}): {', '.join(inactive_genes) if inactive_genes else 'none'}")
            print(f"      Binary: {state_to_string(state)}")
            if j < len(attractor) - 1:  # Add separator between states except for the last one
                print()
    
    print(f"\n{'='*60}")
    if basins:
        print("BASIN SUMMARY:")
        for i, basin in enumerate(basins, 1):
            percentage = (basin / sum(basins)) * 100 if sum(basins) > 0 else 0
            print(f"  Attractor {i}: {basin:6d} states ({percentage:6.2f}%)")
        print(f"  TOTAL: {sum(basins):6d} states")
    print(f"{'='*60}\n")


def print_rules_summary(gene_rules, verbose=False):
    """
    Print the rules being used, highlighting constant rules.
    
    Args:
        gene_rules (dict): Dictionary of gene rules
        verbose (bool): Enable verbose output
    """
    if verbose:
        print("Rules used:")
        for gene, rule in gene_rules.items():
            is_const, const_value = is_constant_rule(rule)
            if is_const:
                value_str = "ON (1)" if const_value else "OFF (0)"
                print(f"  {gene}: {rule} [CONSTANT - {value_str}]")
            else:
                print(f"  {gene}: {rule}")
        print()


def save_selected_rules_tsv(selected_rules, genes, output_path, filename, verbose=False, mutations=None):
    """
    Save selected rules to TSV format for reuse in other scripts.
    Handles constant rules and mutations properly.
    
    Args:
        selected_rules (dict): Dictionary of selected rules per gene
        genes (list): List of gene names
        output_path (str): Output directory path
        filename (str): Output filename
        verbose (bool): Enable verbose output
        mutations (dict): Dictionary of gene mutations
        
    Returns:
        pd.DataFrame: DataFrame with selected rules
    """
    rules_data = []
    mutations = mutations or {}

    for gene in genes:
        # Check if gene is mutated - override with constant value
        if gene in mutations:
            rule = str(mutations[gene])  # Convert to string: "1" or "0"
        elif gene in selected_rules:
            rule = selected_rules[gene]
        else:
            rule = gene  # Self-regulation fallback
        
        rules_data.append({
            'Gene': gene,
            'Rule': rule
        })
    
    # Create DataFrame and save
    df_rules = pd.DataFrame(rules_data)
    full_path = os.path.join(output_path, filename)
    df_rules.to_csv(full_path, sep='\t', index=False)
    
    log_message(f"Selected rules saved to: {filename}", verbose)
    return df_rules


def analyze_attractors(args):
    """
    Main function to analyze Boolean network attractors.
    
    Args:
        args: Parsed command line arguments
    """
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory structure
    if args.output_dir:
        base_output_dir = args.output_dir
    else:
        base_output_dir = os.path.dirname(args.input_file)
    
    # Create attractors subdirectory
    output_dir = os.path.join(base_output_dir, "attractors")
    os.makedirs(output_dir, exist_ok=True)
    log_message(f"Output directory created: {output_dir}", args.verbose)
    
    try:
        # Read rules table
        df = read_rules_table(args.input_file, args.verbose)
        
        # Select rules based on criteria
        gene_rules = select_rules_by_criteria(df, args.criteria, args.positions, args.verbose)
        
        if not gene_rules:
            print("ERROR: No rules selected", file=sys.stderr)
            sys.exit(1)
        
        log_message(f"Rules selected for {len(gene_rules)} genes", args.verbose)
        
        # Extract all genes (excluding boolean constants)
        genes = extract_genes_from_rules(gene_rules)
        log_message(f"Genes in network: {genes}", args.verbose)
        log_message(f"Total genes: {len(genes)}", args.verbose)

        # Parse mutations
        mutations = parse_mutations(args.mutations, genes)
        if mutations:
            log_message(f"Mutations applied: {mutations}", args.verbose)
            print(f"Genes with fixed values: {', '.join([f'{gene}={value}' for gene, value in mutations.items()])}")
        
        # Print rules summary
        print_rules_summary(gene_rules, args.verbose)
        
        # Find attractors with multiprocessing
        if mutations:
            print(f"Calculating attractor states and basins with {len(mutations)} gene mutation(s)...")
        else:
            print("Calculating attractor states and basins with multiprocessing...")
        attractors, basins = find_attractors_parallel(
            gene_rules, genes, 
            max_iterations=args.max_iterations, 
            n_procs=args.processes,
            verbose=args.verbose,
            mutations=mutations
        )
        
        # Print results
        print_attractors_summary(attractors, genes, basins, args.verbose, mutations)

        # Save results
        print("=== SAVING RESULTS ===")
        mutation_suffix = generate_mutation_suffix(mutations)

        # Generate filenames with mutation suffix
        if args.attractors_output:
            attractors_file = args.attractors_output
        else:
            base_name = "attractors"
            attractors_file = f"{base_name}{mutation_suffix}.tsv"

        if args.rules_output:
            rules_file = args.rules_output
        else:
            base_name = "selected_rules"
            rules_file = f"{base_name}{mutation_suffix}.tsv"

        save_attractors_tsv(attractors, genes, basins, output_dir, attractors_file, args.verbose)
        save_selected_rules_tsv(gene_rules, genes, output_dir, rules_file, args.verbose, mutations)

        # Print final summary
        print(f"\n{'='*80}")
        print("PROCESS COMPLETED SUCCESSFULLY")
        print('='*80)
        print(f"Files generated in: {output_dir}")
        print(f"- {attractors_file} (attractor states and cycles)")
        print(f"- {rules_file} (selected rules for reuse)")
        if mutations:
            print(f"- Mutations applied: {', '.join([f'{gene}={value}' for gene, value in mutations.items()])}")
        print('='*80)

    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Analyze Boolean network attractors and basins of attraction from logical rules',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 BNI3_Attractors.py -i rules_by_gene.tsv
  python3 BNI3_Attractors.py -i selected_rules.tsv -o output_dir/ -v
  python3 BNI3_Attractors.py -i rules.tsv -c custom_positions -p "1,1,2,3,1" -v
  python3 BNI3_Attractors.py -i rules.tsv --max-iter 2000 -n 8
  python3 BNI3_Attractors.py -i rules.tsv -m "GENE1:1,GENE2:0" -v
  python3 BNI3_Attractors.py -i rules.tsv -m "ABF3:1" --processes 4

Input File Format:
  The input file should be a TSV with at minimum:
  - Gene: Gene name
  - Rule: Boolean rule (can be logical expression or constant 0/1/TRUE/FALSE)
  - Position: (optional) Rule ranking - if missing, all rules get Position=1

Boolean Constants in Rules:
  Rules can be constants to fix gene states:
  - "1" or "TRUE" = Gene always ON
  - "0" or "FALSE" = Gene always OFF
  Example TSV:
    Gene    Rule
    GENE1   GENE2 & GENE3
    GENE2   FALSE          <- Gene2 always OFF
    GENE3   TRUE           <- Gene3 always ON

Selection Criteria:
  top_position     : Select rule with Position=1 for each gene (default)
  custom_positions : Select rules using custom position list (requires --positions)

Custom Positions Format:
  Use comma-separated integers corresponding to each gene in alphabetical order
  Example: "1,1,2,3,1" means:
  - 1st gene: use Position 1
  - 2nd gene: use Position 1  
  - 3rd gene: use Position 2
  - 4th gene: use Position 3
  - 5th gene: use Position 1

Mutations Format:
  Use comma-separated gene mutations in format "GENE1:VALUE,GENE2:VALUE"
  Where VALUE is:
  - 1 = Gene overexpression (always ON)
  - 0 = Gene knockout (always OFF)
  Example: "ABF3:1,MLP7:0" means:
  - ABF3 is forced to be always active (overexpressed)
  - MLP7 is forced to be always inactive (knockout)
  
  Mutations reduce the state space from 2^n to 2^(n-k) where k is number of mutated genes
  All analysis (attractors, basins) will reflect the mutated network behavior

Notes:
  - Uses synchronous Boolean network dynamics
  - Automatically detects number of CPU cores for parallel processing
  - Large networks (>15 genes) may require significant computation time
  - If specified position not found for a gene, uses minimum available position
  - Creates attractors/ subdirectory for organized output
  - Mutations fix specific genes to constant values, reducing computational complexity
  - Mutated genes maintain their fixed values throughout all simulations
  - Constant rules (0/1/TRUE/FALSE) are functionally equivalent to mutations
        """
    )
    
    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input_file', type=str, required=True,
                         help='Input TSV file containing Boolean rules')
    
    # Optional parameters
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-o', '--output_dir', type=str, default=None,
                         help='Output directory (default: same as input file)')
    optional.add_argument('-c', '--criteria', type=str, default='top_position',
                         choices=['top_position', 'custom_positions'],
                         help='Rule selection criteria (default: top_position)',
                         metavar='{top_position,custom_positions}')
    optional.add_argument('-p', '--positions', type=str, default=None,
                         help='Comma-separated list of positions for custom_positions criteria (e.g., "1,1,2,3,1")',
                         metavar='POSITIONS')
    optional.add_argument('-ao', '--attractors_output', type=str, default=None,
                         help='Attractors output filename (default: attractors.tsv)')
    optional.add_argument('--max-iter', '--max_iterations', type=int, default=1000,
                         dest='max_iterations',
                         help='Maximum iterations for trajectory simulation (default: 1000)')
    optional.add_argument('-n', '--processes', type=int, default=None,
                         help='Number of processes for parallel computation (default: auto-detect)')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')
    optional.add_argument('-ro', '--rules_output', type=str, default=None,
                         help='Selected rules output filename (default: selected_rules.tsv)')
    optional.add_argument('-m', '--mutations', type=str, default=None,
                         help='Gene mutations in format "GENE1:1,GENE2:0" (1=overexpression, 0=knockout)',
                         metavar='MUTATIONS')
    
    parser.add_argument('--version', action='version', version='Boolean Network Attractor Analyzer v1.1')
    
    args = parser.parse_args()
    
    # Validate processes argument
    if args.processes is not None and args.processes < 1:
        print("ERROR: Number of processes must be >= 1", file=sys.stderr)
        sys.exit(1)
    
    # Run analysis
    try:
        analyze_attractors(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()