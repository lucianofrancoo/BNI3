#!/usr/bin/env python3
"""
Generate test datasets for binarization algorithms
Creates TSV files with simulated gene expression data varying in:
- Number of genes (5, 10, 20, 50, 100)
- Number of timepoints (3, 5, 10, 20, 50)
"""

import numpy as np
import pandas as pd
from pathlib import Path

def generate_expression_pattern(n_timepoints, pattern_type='random'):
    """
    Generate different expression patterns
    
    Pattern types:
    - 'random': Random values
    - 'increasing': Gradual increase
    - 'decreasing': Gradual decrease
    - 'peak': Peak in the middle
    - 'valley': Valley in the middle
    - 'constant_low': Constantly low expression
    - 'constant_high': Constantly high expression
    - 'oscillating': Oscillating pattern
    """
    
    if pattern_type == 'random':
        return np.random.uniform(0.01, 0.95, n_timepoints)
    
    elif pattern_type == 'increasing':
        base = np.linspace(0.1, 0.9, n_timepoints)
        noise = np.random.normal(0, 0.05, n_timepoints)
        return np.clip(base + noise, 0.01, 0.99)
    
    elif pattern_type == 'decreasing':
        base = np.linspace(0.9, 0.1, n_timepoints)
        noise = np.random.normal(0, 0.05, n_timepoints)
        return np.clip(base + noise, 0.01, 0.99)
    
    elif pattern_type == 'peak':
        middle = n_timepoints // 2
        x = np.arange(n_timepoints)
        base = 0.9 * np.exp(-((x - middle) ** 2) / (n_timepoints / 2))
        noise = np.random.normal(0, 0.03, n_timepoints)
        return np.clip(base + noise + 0.1, 0.01, 0.99)
    
    elif pattern_type == 'valley':
        middle = n_timepoints // 2
        x = np.arange(n_timepoints)
        base = 0.8 - 0.6 * np.exp(-((x - middle) ** 2) / (n_timepoints / 2))
        noise = np.random.normal(0, 0.03, n_timepoints)
        return np.clip(base + noise + 0.1, 0.01, 0.99)
    
    elif pattern_type == 'constant_low':
        base = np.random.uniform(0.01, 0.15, n_timepoints)
        noise = np.random.normal(0, 0.01, n_timepoints)
        return np.clip(base + noise, 0.001, 0.2)
    
    elif pattern_type == 'constant_high':
        base = np.random.uniform(0.75, 0.95, n_timepoints)
        noise = np.random.normal(0, 0.03, n_timepoints)
        return np.clip(base + noise, 0.7, 0.99)
    
    elif pattern_type == 'oscillating':
        x = np.arange(n_timepoints)
        base = 0.5 + 0.3 * np.sin(2 * np.pi * x / (n_timepoints / 3))
        noise = np.random.normal(0, 0.05, n_timepoints)
        return np.clip(base + noise, 0.01, 0.99)
    
    else:
        return np.random.uniform(0.01, 0.95, n_timepoints)


def generate_dataset(n_genes, n_timepoints, output_file):
    """
    Generate a complete dataset with diverse expression patterns
    """
    
    # Pattern distribution for realism
    pattern_types = ['random', 'increasing', 'decreasing', 'peak', 'valley', 
                     'constant_low', 'constant_high', 'oscillating']
    
    # Create data structure
    data = {}
    
    for i in range(1, n_genes + 1):
        gene_name = f"Gene{i}"
        
        # Randomly select a pattern type for each gene
        pattern = np.random.choice(pattern_types)
        
        # Generate expression values
        expression = generate_expression_pattern(n_timepoints, pattern)
        
        data[gene_name] = expression
    
    # Create DataFrame
    df = pd.DataFrame(data).T
    df.columns = [f"t{i}" for i in range(1, n_timepoints + 1)]
    df.index.name = 'ID'
    
    # Save to file
    df.to_csv(output_file, sep='\t', float_format='%.7f')
    
    print(f"Generated: {output_file.name} ({n_genes} genes × {n_timepoints} timepoints)")


def main():
    """
    Generate multiple test datasets with different dimensions
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate test datasets for binarization algorithms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 generate_test_datasets.py -g 10,20,50 -t 3,5,10 -o ./datasets
  python3 generate_test_datasets.py -g 5,10 -t 4,5 -o /path/to/output
  python3 generate_test_datasets.py --n_genes 10,20,50,100 --n_timepoints 3,4,5,10,20
        """
    )
    
    parser.add_argument('-g', '--n_genes',
                        required=True,
                        help='Comma-separated list of gene numbers (e.g., 10,20,50,100)')
    
    parser.add_argument('-t', '--n_timepoints',
                        required=True,
                        help='Comma-separated list of timepoint numbers (e.g., 3,5,10,20)')
    
    parser.add_argument('-o', '--output',
                        default='./test_datasets',
                        help='Output directory for generated datasets (default: ./test_datasets)')
    
    parser.add_argument('-s', '--seed',
                        type=int,
                        default=42,
                        help='Random seed for reproducibility (default: 42)')
    
    parser.add_argument('--version',
                        action='version',
                        version='Test Dataset Generator v1.0')
    
    args = parser.parse_args()
    
    # Parse gene numbers
    try:
        gene_numbers = [int(x.strip()) for x in args.n_genes.split(',')]
    except ValueError:
        print("ERROR: Invalid format for n_genes. Use comma-separated integers (e.g., 10,20,50)")
        return 1
    
    # Parse timepoint numbers
    try:
        timepoint_numbers = [int(x.strip()) for x in args.n_timepoints.split(',')]
    except ValueError:
        print("ERROR: Invalid format for n_timepoints. Use comma-separated integers (e.g., 3,5,10)")
        return 1
    
    # Validate inputs
    if any(g <= 0 for g in gene_numbers):
        print("ERROR: All gene numbers must be positive integers")
        return 1
    
    if any(t <= 0 for t in timepoint_numbers):
        print("ERROR: All timepoint numbers must be positive integers")
        return 1
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set random seed for reproducibility
    np.random.seed(args.seed)
    
    print("="*60)
    print("GENERATING TEST DATASETS")
    print("="*60)
    print(f"Gene numbers:      {gene_numbers}")
    print(f"Timepoint numbers: {timepoint_numbers}")
    print(f"Output directory:  {output_dir}")
    print(f"Random seed:       {args.seed}")
    print("="*60 + "\n")
    
    # Generate all combinations
    datasets_generated = 0
    for n_genes in gene_numbers:
        for n_timepoints in timepoint_numbers:
            output_file = output_dir / f"test_data_g{n_genes}_t{n_timepoints}.tsv"
            generate_dataset(n_genes, n_timepoints, output_file)
            datasets_generated += 1
    
    print("\n" + "="*60)
    print(f"Total datasets generated: {datasets_generated}")
    print(f"Output directory: {output_dir.absolute()}")
    print("="*60)
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())