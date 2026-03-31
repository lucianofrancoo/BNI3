#!/usr/bin/env python3
"""
GEP Results Analysis Tool
Analyzes results from GEP Boolean network inference and creates summary of best rules per gene.
"""

import argparse
import os
import glob
import re
import pandas as pd
import sys
from pathlib import Path


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def calculate_regulator_score(n_regulators):
    """
    Calculate regulator score based on number of regulators
    Similar to penalty but inverted (higher is better)
    
    Args:
        n_regulators (int): Number of regulators
        
    Returns:
        float: Score from 0 to 1 (1 = optimal, 0 = worst)
    """
    if pd.isna(n_regulators):
        return 0.0
    
    # Inverted penalty: optimal at 2-3 regulators
    score = 1 - ((n_regulators - 2.5) / 2) ** 2
    return max(0, min(1, score))


def detect_mlp_usage(results_dir):
    """
    Detect if MLP was used by checking if any result file contains mlp_loss column
    
    Args:
        results_dir (str): Results directory path
        
    Returns:
        bool: True if MLP was used
    """
    # Find any TSV file in results
    pattern = os.path.join(results_dir, "best_individuals_*.tsv")
    files = glob.glob(pattern)
    
    if not files:
        return False
    
    try:
        # Check first file for mlp_loss column
        df = pd.read_csv(files[0], sep='\t', nrows=1)
        return 'mlp_loss' in df.columns
    except:
        return False


def find_result_files(results_dir, target_genes, verbose):
    """
    Find all result files for target genes
    
    Args:
        results_dir (str): Results directory path
        target_genes (list): List of target gene names
        verbose (bool): Enable verbose output
        
    Returns:
        dict: Dictionary mapping gene names to list of file paths
    """
    log_message("Finding result files...", verbose)
    
    files_by_gene = {}
    
    for gene in target_genes:
        # Find all files for this gene
        pattern = os.path.join(results_dir, f"best_individuals_{gene}_*.tsv")
        files = glob.glob(pattern)
        
        if files:
            files_by_gene[gene] = files
            log_message(f"Found {len(files)} files for gene {gene}", verbose)
        else:
            log_message(f"No files found for gene {gene}", verbose)
    
    return files_by_gene


def collect_all_rules(files_by_gene, use_mlp, verbose):
    """
    Collect all rules from result files
    
    Args:
        files_by_gene (dict): Dictionary mapping genes to file paths
        use_mlp (bool): Whether MLP was used
        verbose (bool): Enable verbose output
        
    Returns:
        list: List of rule dictionaries
    """
    log_message("Collecting rules from all files...", verbose)
    
    all_rules = []
    
    for gene, files in files_by_gene.items():
        log_message(f"Processing {len(files)} files for gene {gene}...", verbose)
        
        for file_path in files:
            try:
                df = pd.read_csv(file_path, sep='\t')
                
                for _, row in df.iterrows():
                    rule_info = {
                        'gene': gene,
                        'expression': row['expression'],
                        'n_correct': int(row['n_correct']),
                        'n_regulators': int(row['n_regulators']),
                        'regulator_penalty': float(row['regulator_penalty']),
                        'mlp_loss': float(row['mlp_loss']) if use_mlp and 'mlp_loss' in row else None
                    }
                    all_rules.append(rule_info)
                    
            except Exception as e:
                log_message(f"Error reading file {file_path}: {str(e)}", verbose)
                continue
    
    log_message(f"Collected {len(all_rules)} total rules", verbose)
    return all_rules


def analyze_rules_by_gene(all_rules, target_genes, use_mlp, n_correct_levels, verbose):
    """
    Analyze rules by gene and create summary
    
    Args:
        all_rules (list): List of all rule dictionaries
        target_genes (list): List of target gene names
        use_mlp (bool): Whether MLP was used
        n_correct_levels (int): Number of n_correct levels to include (from max down)
        verbose (bool): Enable verbose output
        
    Returns:
        list: List of summary dictionaries
    """
    log_message("Analyzing rules by gene...", verbose)
    
    summary_data = []
    
    # First pass: find global maximum n_correct across all genes
    global_max_correct = max(r['n_correct'] for r in all_rules) if all_rules else 0
    log_message(f"Global maximum n_correct found: {global_max_correct}", verbose)
    log_message(f"Including top {n_correct_levels} n_correct levels", verbose)
    
    for gene in target_genes:
        # Get rules for this gene
        gene_rules = [r for r in all_rules if r['gene'] == gene]
        
        if not gene_rules:
            log_message(f"No rules found for gene {gene}", verbose)
            continue
            
        # Find maximum n_correct for this gene
        gene_max_correct = max(r['n_correct'] for r in gene_rules)
        
        # Calculate which n_correct values to include
        target_n_correct_values = []
        for i in range(n_correct_levels):
            n_correct_value = gene_max_correct - i
            if n_correct_value >= 0:  # Don't go below 0
                target_n_correct_values.append(n_correct_value)
        
        # Debug: show what values we're including
        log_message(f"Gene {gene}: max_correct={gene_max_correct}, including n_correct values: {target_n_correct_values}", verbose)
        
        best_rules = [r for r in gene_rules if r['n_correct'] in target_n_correct_values]
        
        # Remove duplicates (same expression)
        unique_rules = {}
        for rule in best_rules:
            expr = rule['expression']
            if expr not in unique_rules:
                unique_rules[expr] = rule
        
        best_rules = list(unique_rules.values())
        
        log_message(f"Gene {gene}: {len(best_rules)} unique rules with n_correct in {target_n_correct_values}", verbose)
        
        # Calculate scores for ranking
        for rule in best_rules:
            # Calculate proportional correct score based on global maximum
            # Formula: 1 - (global_max - current) / global_max
            if global_max_correct > 0:
                correct_score = 1 - (global_max_correct - rule['n_correct']) / global_max_correct
            else:
                correct_score = 1.0
            
            regulator_score = calculate_regulator_score(rule['n_regulators'])
            
            if use_mlp and rule['mlp_loss'] is not None:
                # Convert MSE to score (lower MSE = higher score)
                # Simple normalization: assume MSE range 0-1, could be improved
                mse_score = max(0, 1 - min(1, rule['mlp_loss']))
                total_score = correct_score + regulator_score + mse_score
            else:
                total_score = correct_score + regulator_score
            
            rule['score'] = total_score
            rule['correct_score'] = correct_score  # Store for debugging if needed
        
        # Sort by score (descending)
        best_rules.sort(key=lambda x: x['score'], reverse=True)
        
        # Add to summary
        for i, rule in enumerate(best_rules, 1):
            summary_data.append({
                'Gene': rule['gene'],
                'Position': i,
                'Rule': rule['expression'],
                'Correct': rule['n_correct'],
                'N_Regulators': rule['n_regulators'],
                'MSE': f"{rule['mlp_loss']:.6f}" if rule['mlp_loss'] is not None else "N/A",
                'Score': f"{rule['score']:.4f}"
            })
    
    return summary_data


def save_summary(summary_data, output_file, args, verbose):
    """
    Save summary data to TSV file
    
    Args:
        summary_data (list): List of summary dictionaries
        output_file (str): Output file path
        verbose (bool): Enable verbose output
    """
    if not summary_data:
        log_message("No summary data to save", verbose)
        return
    
    # Create DataFrame
    summary_df = pd.DataFrame(summary_data)
    
    # Save to TSV
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    log_message(f"Summary saved to: {output_file}", verbose)
    
    # Print statistics
    total_rules = len(summary_data)
    unique_genes = len(set(r['Gene'] for r in summary_data))
    
    print(f"\n{'='*60}")
    print("RULES ANALYSIS SUMMARY")
    print('='*60)
    print(f"Genes analyzed:       {unique_genes}")
    print(f"Total best rules:     {total_rules}")
    print(f"Rules per gene:       {total_rules/unique_genes:.1f} average" if unique_genes > 0 else "Rules per gene:       N/A")
    print(f"Output file:          {os.path.basename(output_file)}")
    print('='*60)

    # Construct the command used for this execution
    script_name = "python3 BNI3_Analyze_results.py"
    command_parts = [f"{script_name}"]
    command_parts.append(f"-i {args.input_dir}")
    current_command = " ".join(command_parts)
    
    print(f"\nResults obtained with the script:")
    print(f"{current_command}")
    print(f"\nIf you want to obtain a different result please check {script_name} -h")
    print('='*60)


def analyze_results(args):
    """
    Main function to analyze GEP results
    
    Args:
        args: Parsed command line arguments
    """
    # Validate input directory
    results_dir = os.path.join(args.input_dir, "results")
    if not os.path.exists(results_dir):
        print(f"ERROR: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Get target genes
    if args.target_genes:
        target_genes = [gene.strip() for gene in args.target_genes.split(',')]
    else:
        # Auto-detect genes from file names
        pattern = os.path.join(results_dir, "best_individuals_*.tsv")
        files = glob.glob(pattern)
        gene_set = set()
        
        for file_path in files:
            filename = os.path.basename(file_path)
            # Extract gene name from filename: best_individuals_GENENAME_...
            match = re.match(r'best_individuals_([^_]+)_', filename)
            if match:
                gene_set.add(match.group(1))
        
        target_genes = sorted(list(gene_set))
        
        if not target_genes:
            print("ERROR: No target genes found in results directory", file=sys.stderr)
            sys.exit(1)
    
    log_message(f"Target genes: {', '.join(target_genes)}", args.verbose)
    
    # Detect if MLP was used
    use_mlp = detect_mlp_usage(results_dir)
    log_message(f"MLP usage detected: {use_mlp}", args.verbose)
    
    # Find result files
    files_by_gene = find_result_files(results_dir, target_genes, args.verbose)
    
    if not files_by_gene:
        print("ERROR: No result files found", file=sys.stderr)
        sys.exit(1)
    
    # Collect all rules
    all_rules = collect_all_rules(files_by_gene, use_mlp, args.verbose)
    
    if not all_rules:
        print("ERROR: No rules collected from files", file=sys.stderr)
        sys.exit(1)
    
    # Analyze rules
    summary_data = analyze_rules_by_gene(all_rules, target_genes, use_mlp, args.n_correct_levels, args.verbose)
    
    # Determine output file
    if args.output:
        output_file = args.output
    else:
        output_file = os.path.join(args.input_dir, "rules_by_gene.tsv")
    
    # Save summary
    save_summary(summary_data, output_file, args, args.verbose)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Analyze GEP Boolean network inference results and create summary of best rules per gene',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 analyze_results.py -i results_directory/
  python3 analyze_results.py -i output/ -targets GENE1,GENE2,GENE3 -v
  python3 analyze_results.py -i results/ -o custom_summary.tsv -n 3
  python3 analyze_results.py -i results/ -n 1 -v

Notes:
  - Automatically detects if MLP was used in the experiments
  - Default shows rules with top 2 n_correct levels (max and max-1)
  - Use -n 1 to show only maximum n_correct rules
  - Use -n 3 to show max, max-1, and max-2 n_correct rules
  - If target genes not specified, auto-detects from result file names
  - Creates rules_by_gene.tsv in input directory by default
        """
    )
    
    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input_dir', type=str, required=True,
                         help='Input directory containing GEP results (should contain results/ subdirectory)')
    
    # Optional parameters
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-targets', '--target_genes', type=str, default=None,
                         help='Target genes to analyze (comma-separated). Default: auto-detect from files')
    optional.add_argument('-o', '--output', type=str, default=None,
                         help='Output file path. Default: reglas_por_gen.tsv in input directory')
    optional.add_argument('-n', '--n_correct_levels', type=int, default=2,
                         help='Number of n_correct levels to include from maximum downward (default: 2, meaning max and max-1)')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')
    
    parser.add_argument('--version', action='version', version='GEP Results Analyzer v1.0')
    
    args = parser.parse_args()
    
    # Validate input directory exists
    if not os.path.exists(args.input_dir):
        print(f"ERROR: Input directory '{args.input_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Run analysis
    try:
        analyze_results(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()