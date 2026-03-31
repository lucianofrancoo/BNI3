#!/usr/bin/env python3
"""
GEP-based Boolean Network Inference Tool
Implementation of Gene Expression Programming for inferring Boolean regulatory networks.
"""

import argparse
import re
import os
import operator
import numpy as np
import pandas as pd
import subprocess
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing
from sklearn.metrics import mean_squared_error
import geppy as gep
from deap import creator, base, tools
from multiprocessing import Pool
import matplotlib.pyplot as plt
import sys
import time
from pathlib import Path
import psutil


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def auto_mlp_config(g, t):
    """
    Automatic MLP configuration based on number of regulators (g)
    and number of transitions (t). Compact architectures suitable
    for small biological datasets.
    
    Args:
        g (int): Number of unique regulators
        t (int): Number of transitions (timepoints - 1)
        
    Returns:
        dict: MLP configuration with keys:
            - hidden_layer_sizes (tuple): Architecture
            - solver (str): 'lbfgs' or 'adam'
            - alpha (float): Regularization parameter
            - max_iter (int): Maximum iterations
            - tol (float): Tolerance for optimization
            - learning_rate_init (float): Initial learning rate (for adam)
    
    Notes:
        For small datasets (t ≤ 40) with many regulators (g > 15), LBFGS may 
        produce convergence warnings due to ill-conditioned optimization landscapes.
        This is expected behavior and does not affect fitness evaluation quality.
    """
    # Solver selection
    if g <= 30:
        solver = "lbfgs"
    else:
        solver = "adam"
    
    # Architecture selection based on number of transitions
    if t <= 15:
        layers = (3,)
        tol = 1e-3
        max_iter = 500 if solver == "lbfgs" else 1000
    elif t <= 40:
        layers = (5,)
        tol = 1e-4
        # Increased max_iter for better convergence in small datasets with many regulators
        max_iter = 1200 if solver == "lbfgs" else 2500
    elif t <= 80:
        layers = (10,)
        tol = 1e-4
        max_iter = 2000 if solver == "lbfgs" else 3000
    elif t <= 150:
        layers = (10, 5)
        tol = 5e-5
        max_iter = 3000 if solver == "lbfgs" else 4000
    else:
        layers = (20, 10)
        tol = 5e-5
        max_iter = 4000 if solver == "lbfgs" else 5000
    
    # Regularization based on number of regulators
    if g <= 30:
        alpha = 0.01
    elif g <= 200:
        alpha = 0.05
    else:
        alpha = 0.1
    
    return {
        "hidden_layer_sizes": layers,
        "solver": solver,
        "alpha": alpha,
        "max_iter": max_iter,
        "tol": tol,
        "learning_rate_init": 0.0005 if solver == "adam" else None
    }

def load_and_clean_data(raw_path, binary_path, verbose):
    """
    Load and clean raw and binary data matrices
    
    Args:
        raw_path (str): Path to raw count matrix
        binary_path (str): Path to binary matrix
        verbose (bool): Enable verbose output
        
    Returns:
        tuple: (raw_data, binary_data, use_mlp)
    """
    log_message("Loading and cleaning data...", verbose)
    
    try:
        # Load raw data
        log_message(f"Reading raw data from: {raw_path}", verbose)
        raw_data = pd.read_csv(raw_path, sep="\t", decimal=".")
        raw_data = raw_data.T
        raw_data.columns = raw_data.iloc[0]
        raw_data = raw_data[1:]
        raw_data = raw_data.reset_index(drop=True)
        raw_data = raw_data.apply(pd.to_numeric) 
        raw_data = raw_data.dropna()
        # Clean column names
        raw_data.columns = [re.sub(r'\W|^(?=\d)', '_', col) for col in raw_data.columns]

        # Load binary data
        log_message(f"Reading binary data from: {binary_path}", verbose)
        binary_data = pd.read_csv(binary_path, sep="\t", decimal=".")
        # Clean column names
        binary_data.columns = [re.sub(r'\W|^(?=\d)', '_', col) for col in binary_data.columns]
        
        rows, cols = binary_data.shape
        use_mlp = rows >= 6
        
        log_message(f"Data loaded successfully: {rows} timepoints, {cols} genes", verbose)
        log_message(f"MLP validation enabled: {use_mlp}", verbose)
        
        return raw_data, binary_data, use_mlp
        
    except Exception as e:
        print(f"ERROR: Failed to load data: {str(e)}", file=sys.stderr)
        sys.exit(1)


def get_target_genes(args, gene_list, verbose):
    """
    Get target genes list from arguments or use all genes
    
    Args:
        args: Parsed arguments
        gene_list: List of all available genes
        verbose (bool): Enable verbose output
        
    Returns:
        list: Target genes to process
    """
    if args.target_genes is None:
        target_genes = gene_list
        log_message(f"Using all {len(target_genes)} genes as targets", verbose)
    else:
        target_genes = [gene.strip() for gene in args.target_genes.split(',')]
        # Validate that target genes exist
        missing_genes = [gene for gene in target_genes if gene not in gene_list]
        if missing_genes:
            print(f"ERROR: Target genes not found in data: {missing_genes}", file=sys.stderr)
            print(f"Available genes: {', '.join(gene_list[:10])}{'...' if len(gene_list) > 10 else ''}", file=sys.stderr)
            sys.exit(1)
        log_message(f"Using {len(target_genes)} specified target genes", verbose)
    
    # Show target genes list for debugging
    if verbose:
        log_message("Target genes to process:", verbose)
        for i, gene in enumerate(target_genes, 1):
            log_message(f"  {i:3d}. {gene}", verbose)
    else:
        # Always show target genes summary even without verbose
        if len(target_genes) <= 20:
            print(f"Target genes ({len(target_genes)}): {', '.join(target_genes)}")
        else:
            print(f"Target genes ({len(target_genes)}): {', '.join(target_genes[:20])}... (and {len(target_genes)-20} more)")
    
    return target_genes


def prepare_data_splits(raw_data, binary_data, verbose):
    """
    Prepare data splits for input and output
    
    Args:
        raw_data (pd.DataFrame): Raw count data
        binary_data (pd.DataFrame): Binary data
        verbose (bool): Enable verbose output
        
    Returns:
        tuple: (raw_data_in, raw_data_out, boolean_data_in, boolean_data_out, all_genes)
    """
    log_message("Preparing data splits...", verbose)
    
    rows = binary_data.shape[0]
    
    # Binary splits
    binary_data_in = binary_data.loc[0:rows-2]
    binary_data_out = binary_data.loc[1:rows-1]
    
    # Convert to boolean
    boolean_data_in = binary_data_in == 1
    boolean_data_out = binary_data_out == 1
    
    # Raw data splits
    raw_data_in = raw_data.loc[0:rows-2]
    raw_data_out = raw_data.loc[1:rows-1]
    
    # Gene list
    all_genes = list(binary_data.columns)
    
    log_message(f"Data splits created: {len(all_genes)} genes, {rows-1} transitions", verbose)
    
    return raw_data_in, raw_data_out, boolean_data_in, boolean_data_out, all_genes


def setup_gep_primitives(all_genes, verbose):
    """
    Configure GEP primitives
    
    Args:
        all_genes (list): List of all gene names
        verbose (bool): Enable verbose output
        
    Returns:
        gep.PrimitiveSet: Configured primitive set
    """
    log_message("Setting up GEP primitives...", verbose)
    
    pset = gep.PrimitiveSet('Main', input_names=all_genes)
    pset.add_function(operator.and_, 2)
    pset.add_function(operator.or_, 2)
    pset.add_function(operator.not_, 1)
    
    log_message(f"Primitives configured: {len(all_genes)} inputs + AND, OR, NOT functions", verbose)
    
    return pset


def setup_fitness_and_creator(weights, use_mlp, use_reg_penalty, verbose):
    """
    Configure fitness function based on MLP usage and regulator penalty
    
    Args:
        weights (tuple): Fitness weights
        use_mlp (bool): Whether to use MLP validation
        use_reg_penalty (bool): Whether to use regulator penalty
        verbose (bool): Enable verbose output
    """
    log_message("Setting up fitness function...", verbose)
    
    # Check if classes already exist and delete them if they do
    if hasattr(creator, 'FitnessMulti'):
        del creator.FitnessMulti
    if hasattr(creator, 'Individual'):
        del creator.Individual
    
    # Determine number of objectives based on configuration
    if use_mlp and use_reg_penalty:
        # All 3 objectives
        creator.create("FitnessMulti", base.Fitness, weights=weights)
        log_message(f"Three-objective fitness (correct + reg_penalty + MLP): weights = {weights}", verbose)
    elif use_mlp and not use_reg_penalty:
        # Only correct predictions and MLP
        weights_no_reg = (weights[0], weights[2])
        creator.create("FitnessMulti", base.Fitness, weights=weights_no_reg)
        log_message(f"Two-objective fitness (correct + MLP): weights = {weights_no_reg}", verbose)
    elif not use_mlp and use_reg_penalty:
        # Only correct predictions and regulator penalty
        weights_no_mlp = weights[:2]
        creator.create("FitnessMulti", base.Fitness, weights=weights_no_mlp)
        log_message(f"Two-objective fitness (correct + reg_penalty): weights = {weights_no_mlp}", verbose)
    else:
        # Only correct predictions
        weights_single = (weights[0],)
        creator.create("FitnessMulti", base.Fitness, weights=weights_single)
        log_message(f"Single-objective fitness (correct only): weights = {weights_single}", verbose)
    
    # Create Individual class
    creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMulti)


def mlp_evaluation(raw_data_in, target_raw_data_out, regulators, manual_config, verbose):
    """
    Evaluate regulator set using MLP with automatic configuration
    
    Args:
        raw_data_in (pd.DataFrame): Input raw data
        target_raw_data_out (np.array): Target output data
        regulators (list): List of regulator gene names
        manual_config (dict or None): Manual MLP configuration to override auto-config
        verbose (bool): Enable verbose output
        
    Returns:
        float: Mean squared error
    """
    try:
        # Get unique regulators and calculate dimensions
        unique_regulators = list(set(regulators))
        g = len(unique_regulators)  # Number of regulators
        t = len(raw_data_in)  # Number of transitions
        
        # Get MLP configuration (auto or manual)
        if manual_config is not None:
            config = manual_config
            if verbose:
                log_message(f"Using MANUAL MLP config: {config}", verbose)
        else:
            config = auto_mlp_config(g, t)
            if verbose:
                log_message(f"Auto MLP config (g={g}, t={t}): {config}", verbose)
        
        # Prepare data
        X = np.array(raw_data_in[unique_regulators])
        
        # Separate scalers for X and y
        scaler_x = preprocessing.MinMaxScaler()
        X_normalized = scaler_x.fit_transform(X)
        
        scaler_y = preprocessing.MinMaxScaler()
        y_normalized = scaler_y.fit_transform(target_raw_data_out).ravel()
        
        # Create MLP model with configuration
        mlp_params = {
            'hidden_layer_sizes': config['hidden_layer_sizes'],
            'activation': 'relu',
            'solver': config['solver'],
            'alpha': config['alpha'],
            'max_iter': config['max_iter'],
            'tol': config['tol'],
            'random_state': 42,
            'verbose': False
        }
        
        # Add learning_rate_init only for adam solver
        if config['solver'] == 'adam':
            mlp_params['learning_rate_init'] = config['learning_rate_init']
            mlp_params['early_stopping'] = True
            mlp_params['validation_fraction'] = 0.2
            mlp_params['n_iter_no_change'] = 50
        
        model = MLPRegressor(**mlp_params)
        
        # Train and predict
        model.fit(X_normalized, y_normalized)
        predictions = model.predict(X_normalized)
        mse = mean_squared_error(predictions, y_normalized)
        
        return mse
        
    except Exception as e:
        log_message(f"Warning: MLP evaluation failed: {str(e)}", verbose)
        return 1.0  # Return high error if MLP fails


def calculate_regulator_penalty(n_regulators, optimal_k, verbose):
    """
    Calculate penalty factor based on deviation from optimal number of regulators.
    Uses symmetric linear penalty: penalty = 0.3 × |n_regulators - optimal_k|
    
    Args:
        n_regulators (int): Number of unique regulators
        optimal_k (int): Optimal number of regulators (user-defined)
        verbose (bool): Enable verbose output
        
    Returns:
        float: Penalty factor (0 = optimal, higher = more deviation)
        
    Penalty scheme (symmetric):
        - n = optimal_k: penalty = 0.0 (optimal)
        - Each regulator away from optimal_k adds 0.3 to penalty
        
    Examples:
        optimal_k=2, n=2 → penalty = 0.0
        optimal_k=2, n=1 → penalty = 0.3
        optimal_k=2, n=3 → penalty = 0.3
        optimal_k=2, n=5 → penalty = 0.9
        optimal_k=3, n=1 → penalty = 0.6
        optimal_k=3, n=3 → penalty = 0.0
        optimal_k=3, n=6 → penalty = 0.9
    """
    deviation = abs(n_regulators - optimal_k)
    penalty = 0.3 * deviation
    
    log_message(f"Regulator penalty: {n_regulators} regulators (optimal={optimal_k}) → "
               f"deviation={deviation} → penalty={penalty:.3f}", verbose)
    
    return penalty


# Global variables for evaluation function (needed for multiprocessing)
eval_data = {}


def evaluate_individual(individual):
    """
    Global evaluation function for individuals (needed for multiprocessing)
    Uses global eval_data dictionary to access necessary data
    """
    # Extract regulators
    regulators = []
    for gene in individual:
        input_variables = gene.kexpression
        gene_regulators = [item.name for item in input_variables 
                         if item.name not in ['and_', 'or_', 'not_']]
        regulators.extend(gene_regulators)

    n_regulators = len(set(regulators))
    
    # Calculate regulator penalty if enabled
    if eval_data['use_reg_penalty']:
        regulator_penalty = calculate_regulator_penalty(
            n_regulators, 
            eval_data['optimal_k'], 
            False  # Disable verbose for individual evaluations
        )
    else:
        regulator_penalty = 0.0
    
    # Compile individual
    func = eval_data['toolbox'].compile(individual)

    # Calculate correct predictions (first objective)
    n_correct = 0
    for input_vals, expected_output in zip(eval_data['vector_boolean_data_in'], eval_data['target_boolean_out_data']):
        try:
            prediction = func(*input_vals)
            if prediction == expected_output:
                n_correct += 1
        except:
            # If evaluation fails, penalize heavily
            n_correct = 0
            break

    # Determine return values based on configuration
    if eval_data['use_mlp'] and eval_data['use_reg_penalty']:
        # All 3 objectives
        if len(set(regulators)) > 0:
            mse = mlp_evaluation(eval_data['raw_data_in'], eval_data['target_raw_data_out'], 
                               regulators, eval_data['manual_mlp_config'], eval_data['verbose'])
            scaling_factor = len(eval_data['vector_boolean_data_in'])
            mse_scaled = mse * scaling_factor
        else:
            mse_scaled = 1.0 * len(eval_data['vector_boolean_data_in'])
        return (n_correct, regulator_penalty, mse_scaled)
    
    elif eval_data['use_mlp'] and not eval_data['use_reg_penalty']:
        # Only correct predictions and MLP
        if len(set(regulators)) > 0:
            mse = mlp_evaluation(eval_data['raw_data_in'], eval_data['target_raw_data_out'], 
                               regulators, eval_data['manual_mlp_config'], eval_data['verbose'])
            scaling_factor = len(eval_data['vector_boolean_data_in'])
            mse_scaled = mse * scaling_factor
        else:
            mse_scaled = 1.0 * len(eval_data['vector_boolean_data_in'])
        return (n_correct, mse_scaled)
    
    elif not eval_data['use_mlp'] and eval_data['use_reg_penalty']:
        # Only correct predictions and regulator penalty
        return (n_correct, regulator_penalty)
    
    else:
        # Only correct predictions
        return (n_correct,)


def create_evaluation_function(vector_boolean_data_in, target_boolean_out_data, 
                              raw_data_in, target_raw_data_out, use_mlp, use_reg_penalty,
                              optimal_k, manual_mlp_config, verbose):
    """
    Create evaluation function specific for this experiment
    
    Returns:
        function: Evaluation function for individuals
    """
    # Store data in global dictionary for multiprocessing
    global eval_data
    eval_data = {
        'vector_boolean_data_in': vector_boolean_data_in,
        'target_boolean_out_data': target_boolean_out_data,
        'raw_data_in': raw_data_in,
        'target_raw_data_out': target_raw_data_out,
        'use_mlp': use_mlp,
        'use_reg_penalty': use_reg_penalty,
        'optimal_k': optimal_k,
        'manual_mlp_config': manual_mlp_config,
        'verbose': verbose,
        'toolbox': None  # Will be set later
    }
    
    return evaluate_individual


def setup_toolbox(pset, head_length, n_genes, verbose):
    """
    Configure GEP toolbox
    
    Args:
        pset: Primitive set
        head_length (int): Head length of chromosomes
        n_genes (int): Number of genes per chromosome
        verbose (bool): Enable verbose output
        
    Returns:
        gep.Toolbox: Configured toolbox
    """
    log_message("Setting up GEP toolbox...", verbose)
    
    toolbox = gep.Toolbox()
    toolbox.register('gene_gen', gep.Gene, pset=pset, head_length=head_length)
    toolbox.register('individual', creator.Individual, gene_gen=toolbox.gene_gen, 
                    n_genes=n_genes, linker=None)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    # Genetic operators
    toolbox.register('select', tools.selTournament, tournsize=3)
    toolbox.register('mut_uniform', gep.mutate_uniform, pset=pset, 
                    ind_pb=2 / (2 * head_length + 1))
    toolbox.pbs['mut_uniform'] = 0.5
    toolbox.register('mut_invert', gep.invert, pb=0.5)
    toolbox.register('mut_is_ts', gep.is_transpose, pb=0.5)
    toolbox.register('mut_ris_ts', gep.ris_transpose, pb=0.5)
    toolbox.register('mut_gene_ts', gep.gene_transpose, pb=0.5)
    toolbox.register('cx_1p', gep.crossover_one_point, pb=0.5)
    toolbox.register('cx_2p', gep.crossover_two_point, pb=0.5)
    toolbox.register('cx_gene', gep.crossover_gene, pb=0.5)
    
    toolbox.register('compile', gep.compile_, pset=pset)
    
    log_message(f"Toolbox configured with head length: {head_length}, genes per individual: {n_genes}", verbose)
    
    return toolbox


def save_evolution_metrics(log, output_dir, target_gene, args, use_mlp, use_reg_penalty, rep_num, verbose):
    """
    Save evolution metrics to file
    
    Args:
        log: Evolution log from GEP
        output_dir (str): Output directory
        target_gene (str): Target gene name
        args: Parsed arguments
        use_mlp (bool): Whether MLP was used
        use_reg_penalty (bool): Whether regulator penalty was used
        rep_num (int): Repetition number
        verbose (bool): Enable verbose output
    """
    log_message(f"Saving evolution metrics for {target_gene} (rep {rep_num})...", verbose)
    
    evolution_metrics = []
    
    for gen_data in log:
        if 'avg' in gen_data:
            metric = {
                'generation': gen_data['gen'],
                'avg_correct': round(gen_data['avg'][0], 2),
                'std_correct': round(gen_data['std'][0], 2),
                'min_correct': int(gen_data['min'][0]),
                'max_correct': int(gen_data['max'][0])
            }
            
            # Add metrics based on configuration
            idx = 1
            if use_reg_penalty:
                metric['avg_regulator_penalty'] = round(gen_data['avg'][idx], 3)
                metric['min_regulator_penalty'] = round(gen_data['min'][idx], 1)
                metric['max_regulator_penalty'] = round(gen_data['max'][idx], 1)
                idx += 1
            
            if use_mlp:
                metric['avg_mlp_loss'] = round(gen_data['avg'][idx], 4)
                metric['min_mlp_loss'] = round(gen_data['min'][idx], 4)
                metric['max_mlp_loss'] = round(gen_data['max'][idx], 4)
            
            evolution_metrics.append(metric)

    evolution_df = pd.DataFrame(evolution_metrics)
    
    # Format columns for consistent decimal places
    evolution_df['avg_correct'] = evolution_df['avg_correct'].map('{:.2f}'.format)
    evolution_df['std_correct'] = evolution_df['std_correct'].map('{:.2f}'.format)
    
    if use_reg_penalty:
        evolution_df['avg_regulator_penalty'] = evolution_df['avg_regulator_penalty'].map('{:.3f}'.format)
        evolution_df['min_regulator_penalty'] = evolution_df['min_regulator_penalty'].map('{:.1f}'.format)
        evolution_df['max_regulator_penalty'] = evolution_df['max_regulator_penalty'].map('{:.1f}'.format)
    
    if use_mlp:
        evolution_df['avg_mlp_loss'] = evolution_df['avg_mlp_loss'].map('{:.4f}'.format)
        evolution_df['min_mlp_loss'] = evolution_df['min_mlp_loss'].map('{:.4f}'.format)
        evolution_df['max_mlp_loss'] = evolution_df['max_mlp_loss'].map('{:.4f}'.format)
    
    # Create metrics directory
    metrics_dir = os.path.join(output_dir, "metrics")
    os.makedirs(metrics_dir, exist_ok=True)
    
    # Save metrics
    metrics_file = os.path.join(metrics_dir, 
                               f"evolution_metrics_{target_gene}_"
                               f"pop{args.population}_gen{args.generations}_rep{rep_num:03d}.tsv")
    
    evolution_df.to_csv(metrics_file, sep="\t", index=False)
    log_message(f"Evolution metrics saved: {metrics_file}", verbose)


def save_best_individuals(population, output_dir, target_gene, args, use_mlp, use_reg_penalty, 
                         optimal_k, rep_num, verbose):
    """
    Save best individuals from final population
    
    Args:
        population: Final population from evolution
        output_dir (str): Output directory
        target_gene (str): Target gene name
        args: Parsed arguments
        use_mlp (bool): Whether MLP was used
        use_reg_penalty (bool): Whether regulator penalty was used
        optimal_k (int): Optimal number of regulators
        rep_num (int): Repetition number
        verbose (bool): Enable verbose output
    """
    log_message(f"Saving best individuals for {target_gene} (rep {rep_num})...", verbose)
    
    # Sort population by fitness based on configuration
    if use_mlp and use_reg_penalty:
        # All 3 objectives
        sorted_pop = sorted(population, key=lambda ind: 
                          (-ind.fitness.values[0], ind.fitness.values[1], ind.fitness.values[2]))
    elif use_mlp and not use_reg_penalty:
        # Correct predictions and MLP only
        sorted_pop = sorted(population, key=lambda ind: 
                          (-ind.fitness.values[0], ind.fitness.values[1]))
    elif not use_mlp and use_reg_penalty:
        # Correct predictions and regulator penalty only
        sorted_pop = sorted(population, key=lambda ind: 
                          (-ind.fitness.values[0], ind.fitness.values[1]))
    else:
        # Only correct predictions
        sorted_pop = sorted(population, key=lambda ind: -ind.fitness.values[0])
    
    best_individuals_data = []
    
    for i, individual in enumerate(sorted_pop):
        try:
            # Simplify expression
            simplified = gep.simplify(individual)
            expression = str(simplified)
            
            # Extract regulators from simplified expression
            regulators_found = re.findall(r'\b[a-zA-Z_]\w*\b', expression)
            unique_regulators = sorted(set(regulators_found) - {'and_', 'or_', 'not_'})
            
            # Get actual number of regulators from the expression
            actual_n_regulators = len(unique_regulators)
            
            # Calculate regulator penalty for display (even if not used in fitness)
            regulator_penalty = calculate_regulator_penalty(actual_n_regulators, optimal_k, False)

            individual_info = {
                'rank': i + 1,
                'expression': expression,
                'n_correct': int(individual.fitness.values[0]),
                'n_regulators': actual_n_regulators,
                'regulator_penalty': regulator_penalty,
                'regulators': ', '.join(unique_regulators),
                'n_unique_regulators': len(unique_regulators),
                'repetition': rep_num
            }
            
            # Add MLP loss if used
            if use_mlp and use_reg_penalty:
                individual_info['mlp_loss'] = float(individual.fitness.values[2])
            elif use_mlp and not use_reg_penalty:
                individual_info['mlp_loss'] = float(individual.fitness.values[1])
            
            best_individuals_data.append(individual_info)
        except Exception as e:
            log_message(f"Warning: Could not process individual {i+1}: {str(e)}", verbose)
            continue

    best_df = pd.DataFrame(best_individuals_data)
    
    # Create results directory
    results_dir = os.path.join(output_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Save results
    results_file = os.path.join(results_dir,
                               f"best_individuals_{target_gene}_"
                               f"pop{args.population}_gen{args.generations}_rep{rep_num:03d}.tsv")
    
    best_df.to_csv(results_file, sep="\t", index=False)
    log_message(f"Best individuals saved: {results_file}", verbose)
    
    # Print top 3 to console if verbose
    if verbose and len(best_df) > 0:
        log_message(f"TOP 3 INDIVIDUALS for {target_gene} (rep {rep_num}):", verbose)
        for _, row in best_df.head(3).iterrows():
            log_message(f"  Rank {row['rank']}: {row['expression']}", verbose)
            fitness_str = f"    Fitness: {row['n_correct']} correct"
            if use_reg_penalty:
                fitness_str += f", {row['n_regulators']} regulators (penalty: {row['regulator_penalty']})"
            if use_mlp:
                fitness_str += f", MLP loss: {row['mlp_loss']:.4f}"
            log_message(fitness_str, verbose)


def run_automatic_analysis(output_dir, target_genes, verbose):
    """
    Run automatic analysis after inference experiments complete
    
    Args:
        output_dir (str): Output directory containing results
        target_genes (list): List of target genes (None for auto-detect)
        verbose (bool): Enable verbose output
    """
    log_message("Starting automatic results analysis...", verbose)
    
    # Look for the analysis script in the same directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    analysis_script = os.path.join(script_dir, "BNI3_Analyze_results.py")
    
    if not os.path.exists(analysis_script):
        print(f"WARNING: Analysis script not found at {analysis_script}")
        print("Please run analysis manually with: python3 BNI3_Analyze_results.py -i <output_directory>")
        return False
    
    # Prepare command
    cmd = [
        "python3", analysis_script,
        "-i", output_dir
    ]
    
    # Add target genes if specific ones were used
    if target_genes:
        cmd.extend(["-targets", ",".join(target_genes)])
    
    # Add verbose flag if enabled
    if verbose:
        cmd.append("-v")

    try:
        print(f"\nRunning automatic analysis...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 minute timeout
        
        if result.returncode == 0:
            print("Automatic analysis completed successfully!")
            # Print the analysis output
            if result.stdout:
                print(result.stdout)
            return True
        else:
            print(f"Analysis failed with return code {result.returncode}")
            if result.stderr:
                print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("Analysis timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"Error running analysis: {str(e)}")
        return False


def run_single_evolution(target_gene, raw_data_in, raw_data_out, boolean_data_in, 
                        boolean_data_out, all_genes, args, use_mlp, rep_num, verbose):
    """
    Run a single evolution for one target gene
    
    Args:
        target_gene (str): Target gene to model
        raw_data_in, raw_data_out: Raw data splits
        boolean_data_in, boolean_data_out: Boolean data splits
        all_genes (list): List of all genes
        args: Parsed arguments
        use_mlp (bool): Whether to use MLP validation
        rep_num (int): Repetition number
        verbose (bool): Enable verbose output
        
    Returns:
        tuple: (final_population, evolution_log)
    """
    log_message(f"Starting evolution for {target_gene} (repetition {rep_num})...", verbose)
    
    # Prepare target-specific data
    vector_boolean_data_in = boolean_data_in[all_genes].values.tolist()
    target_boolean_out_data = boolean_data_out[target_gene].tolist()
    target_raw_data_out = np.array(raw_data_out[target_gene].tolist()).reshape(-1, 1)
    
    # Parse weights
    weights = tuple(map(float, args.fitness_weights.split(',')))
    
    # Prepare manual MLP config if provided
    manual_mlp_config = None
    if args.mlp_layers is not None or args.mlp_alpha is not None:
        manual_mlp_config = {}
        
        if args.mlp_layers is not None:
            manual_mlp_config['hidden_layer_sizes'] = tuple(map(int, args.mlp_layers.split(',')))
        
        if args.mlp_alpha is not None:
            manual_mlp_config['alpha'] = args.mlp_alpha
        
        # Set defaults for other parameters when manual config is used
        manual_mlp_config.setdefault('solver', 'lbfgs')
        manual_mlp_config.setdefault('max_iter', 5000)
        manual_mlp_config.setdefault('tol', 1e-4)
        manual_mlp_config.setdefault('learning_rate_init', 0.001)
        
        log_message(f"Using manual MLP configuration: {manual_mlp_config}", verbose)
    else:
        log_message("Using automatic MLP configuration", verbose)
    
    log_message(f"Evolution parameters: pop={args.population}, gen={args.generations}, "
               f"weights={weights}, reg_penalty={args.reg_penalty}, optimal_k={args.reg_optimal}", verbose)
    
    # Setup GEP components
    pset = setup_gep_primitives(all_genes, verbose)
    setup_fitness_and_creator(weights, use_mlp, args.reg_penalty, verbose)
    
    global toolbox  # Make toolbox globally accessible
    toolbox = setup_toolbox(pset, args.individual_head, args.individual_genes, verbose)
    
    # Create evaluation function
    evaluate_func = create_evaluation_function(
        vector_boolean_data_in, target_boolean_out_data,
        raw_data_in, target_raw_data_out, use_mlp, args.reg_penalty,
        args.reg_optimal, manual_mlp_config, verbose
    )
    
    # Set toolbox reference in global eval_data
    eval_data['toolbox'] = toolbox
    toolbox.register('evaluate', evaluate_func)
    
    # Setup statistics
    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register('avg', lambda x: np.mean(x, axis=0))
    stats.register('std', lambda x: np.std(x, axis=0))
    stats.register('min', lambda x: np.min(x, axis=0))
    stats.register('max', lambda x: np.max(x, axis=0))
    
    # Initialize population
    pop = toolbox.population(n=args.population)
    hof = tools.HallOfFame(args.elites)
    
    # Setup parallel processing with CPU affinity
    try:
        process = psutil.Process(os.getpid())
        # Limit to first n_processors cores
        available_cores = list(range(args.processes))
        process.cpu_affinity(available_cores)
        log_message(f"CPU affinity set to cores: {available_cores}", verbose)
    except Exception as e:
        log_message(f"Warning: Could not set CPU affinity: {e}", verbose)
    
    pool = Pool(processes=args.processes)
    toolbox.register("map", pool.map)
    
    log_message(f"Running evolution with {args.processes} processes on cores {list(range(args.processes))}...", verbose)
    
    try:
        # Run evolution
        final_pop, log = gep.gep_simple(
            pop, toolbox,
            n_generations=args.generations,
            n_elites=args.elites,
            stats=stats,
            hall_of_fame=hof,
            verbose=False
        )
        
        pool.close()
        pool.join()
        
        log_message(f"Evolution completed for {target_gene} (rep {rep_num})", verbose)
        
        return final_pop, log
        
    except Exception as e:
        pool.close()
        pool.join()
        print(f"ERROR: Evolution failed for {target_gene} (rep {rep_num}): {str(e)}", file=sys.stderr)
        return None, None


def run_gep_experiment(args):
    """
    Main function to run complete GEP experiment
    
    Args:
        args: Parsed command line arguments
    """
    start_time = time.time()
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    log_message(f"Output directory created: {args.output}", args.verbose)
    
    # Load and prepare data
    raw_data, binary_data, use_mlp = load_and_clean_data(args.input, args.input_binary, args.verbose)
    raw_data_in, raw_data_out, boolean_data_in, boolean_data_out, all_genes = prepare_data_splits(
        raw_data, binary_data, args.verbose)
    
    # Get target genes
    target_genes = get_target_genes(args, all_genes, args.verbose)

    # Calculate example MLP config for display (if MLP is enabled)
    example_mlp_config = None
    if use_mlp and args.mlp_layers is None and args.mlp_alpha is None:
        # Use median number of genes as estimate for g
        example_g = len(all_genes)
        example_t = len(raw_data_in) + 1
        example_mlp_config = auto_mlp_config(example_g, example_t)
    
    total_experiments = len(target_genes) * args.n_repetitions
    
    print(f"\n{'='*60}")
    print("GEP BOOLEAN NETWORK INFERENCE")
    print('='*60)
    print(f"Input files:       {os.path.basename(args.input)} + {os.path.basename(args.input_binary)}")
    print(f"Output directory:  {args.output}")
    print(f"Target genes:      {len(target_genes)}")
    print(f"Repetitions:       {args.n_repetitions}")
    print(f"Total experiments: {total_experiments}")
    print(f"Population size:   {args.population}")
    print(f"Generations:       {args.generations}")
    print(f"Processes:         {args.processes}")
    print(f"Regulator penalty: {'Enabled' if args.reg_penalty else 'Disabled'}")
    if args.reg_penalty:
        print(f"  Optimal K:       {args.reg_optimal}")
        print(f"  Penalty formula: 0.3 × |n_regulators - {args.reg_optimal}|")
    print(f"MLP validation:    {'Enabled' if use_mlp else 'Disabled'}")
    if use_mlp:
        if args.mlp_layers is None and args.mlp_alpha is None:
            print(f"MLP configuration: AUTO (adaptive)")
            if example_mlp_config:
                print(f"  - architecture:  {example_mlp_config['hidden_layer_sizes']} (Automatic for {example_g} genes and {example_t} time points)")
                print(f"  - solver:        {example_mlp_config['solver']}")
                print(f"  - alpha:         {example_mlp_config['alpha']}")
                print(f"  - max_iter:      {example_mlp_config['max_iter']}")
                print(f"  - tol:           {example_mlp_config['tol']}")
        else:
            print(f"MLP configuration: MANUAL")
            if args.mlp_layers:
                print(f"  - layers:        {args.mlp_layers}")
            if args.mlp_alpha:
                print(f"  - alpha:         {args.mlp_alpha}")
            print(f"  - solver:        lbfgs (default)")
            print(f"  - max_iter:      5000 (default)")
            print(f"  - tol:           1e-4 (default)")
    # Configuration parameters
    print(f"Individual head:   {args.individual_head}")
    print(f"Individual genes:  {args.individual_genes}")
    print(f"Fitness weights:   {args.fitness_weights}")
    print(f"Elite individuals: {args.elites}")
    print('='*60)
    
    # Run experiments for each target gene
    experiment_count = 0
    failed_experiments = 0
    
    for target_gene in target_genes:
        log_message(f"\nProcessing target gene: {target_gene}", args.verbose)
        
        for rep in range(args.n_repetitions):
            experiment_count += 1
            
            print(f"\nExperiment {experiment_count}/{total_experiments}: "
                  f"{target_gene} | Rep {rep+1}/{args.n_repetitions}")
            
            # Run evolution
            final_pop, log = run_single_evolution(
                target_gene, raw_data_in, raw_data_out, 
                boolean_data_in, boolean_data_out, all_genes, 
                args, use_mlp, rep+1, args.verbose
            )
            
            if final_pop is not None and log is not None:
                # Save results
                save_evolution_metrics(log, args.output, target_gene, args, use_mlp, 
                                     args.reg_penalty, rep+1, args.verbose)
                save_best_individuals(final_pop, args.output, target_gene, args, use_mlp, 
                                    args.reg_penalty, args.reg_optimal, rep+1, args.verbose)
            else:
                failed_experiments += 1
                print(f"FAILED: Experiment {experiment_count}")
    
    # Final summary
    total_time = time.time() - start_time
    successful_experiments = total_experiments - failed_experiments
    
    print(f"\n{'='*60}")
    print("EXPERIMENT SUMMARY")
    print('='*60)
    print(f"Total experiments:    {total_experiments}")
    print(f"Successful:           {successful_experiments}")
    print(f"Failed:               {failed_experiments}")
    print(f"Success rate:         {100*successful_experiments/total_experiments:.1f}%")
    print(f"Total time:           {total_time:.2f} seconds")
    print(f"Average time/exp:     {total_time/total_experiments:.2f} seconds")
    print(f"Results saved in:     {args.output}")
    print('='*60)
    
    if failed_experiments == 0:
        print("STATUS: ALL EXPERIMENTS COMPLETED SUCCESSFULLY")
        
        # Run automatic analysis unless disabled
        if not args.no_analysis:
            analysis_success = run_automatic_analysis(args.output, target_genes, args.verbose)
            if not analysis_success:
                print("\nAutomatic analysis failed. You can run it manually with:")
                print(f"python3 BNI3_Analyze_results.py -i {args.output}")
        else:
            print("\nAutomatic analysis skipped. You can run it manually with:")
            print(f"python3 BNI3_Analyze_results.py -i {args.output}")
            
    else:
        print(f"STATUS: {failed_experiments} EXPERIMENTS FAILED")
        print("\nSkipping analysis due to failed experiments")
        print("Fix issues and run analysis manually if needed:")
        print(f"python3 BNI3_Analyze_results.py -i {args.output}")

def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='GEP-based Boolean Network Inference for Gene Regulatory Networks',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with automatic MLP configuration and default K=2
  python3 gep_algorithm.py -i counts.tsv -i_binary binary.tsv -o results/
  
  # Disable regulator penalty (free exploration)
  python3 gep_algorithm.py -i data.tsv -i_binary data_bin.tsv -o out/ --no_reg_penalty
  
  # Specific target genes with more repetitions
  python3 gep_algorithm.py -i data.tsv -i_binary data_bin.tsv -o out/ -targets GENE1,GENE2 -n_rep 50 -v
  
  # Custom population and generations with parallel processing
  python3 gep_algorithm.py -i data.tsv -i_binary data_bin.tsv -o out/ -pop 200 -gen 150 -p 8
  
  # Manual MLP configuration (overrides automatic configuration)
  python3 gep_algorithm.py -i data.tsv -i_binary data_bin.tsv -o out/ -mlp_layers 20,10 -mlp_alpha 0.5
  
  # Skip automatic analysis
  python3 gep_algorithm.py -i data.tsv -i_binary data_bin.tsv -o out/ --no_analysis

Notes:
  - For small datasets (3-4 time points), consider using -n_rep 50 or higher
  - MLP validation is automatically enabled for datasets with 6+ time points
  - MLP configuration is AUTOMATIC by default (adapts to dataset size and regulator count)
  - Use -mlp_layers and/or -mlp_alpha to manually override automatic MLP configuration
  - Regulator penalty is ENABLED by default with optimal K=2 (Kauffman principle)
  - Use --no_reg_penalty to disable regulator penalty for free exploration
  - Use -reg_optimal to set different optimal number of regulators (e.g., K=3 for stress, K=5 for cell cycle)
  - Target genes are taken from column names of the input matrices
  - Use -v flag for detailed verbose output during processing
  - Analysis runs automatically after successful completion unless --no_analysis is used
  - Analysis creates rules_by_gene.tsv in the output directory
        """
    )
    
    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input', type=str, required=True,
                         help='Input count matrix file (tab-separated)')
    required.add_argument('-i_binary', '--input_binary', type=str, required=True,
                         help='Input binarized matrix file (tab-separated)')
    required.add_argument('-o', '--output', type=str, required=True,
                         help='Output directory (will be created if it doesn\'t exist)')
    
    # Flexible parameters with defaults
    flexible = parser.add_argument_group('Flexible parameters')
    flexible.add_argument('-targets', '--target_genes', type=str, default=None,
                         help='Target genes to model (comma-separated). Default: all genes from input matrix')
    flexible.add_argument('-p', '--processes', type=int, default=4,
                         help='Number of processes for parallel computation (default: 4)')
    flexible.add_argument('-n_rep', '--n_repetitions', type=int, default=10,
                         help='Number of evolution repetitions. For small datasets (3-4 timepoints) consider 50+ (default: 10)')
    flexible.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')
    flexible.add_argument('--no_analysis', action='store_true',
                         help='Skip automatic analysis after inference completion')
    
    # Regulator penalty parameters
    reg_params = parser.add_argument_group('Regulator penalty parameters')
    reg_params.add_argument('--no_reg_penalty', dest='reg_penalty', action='store_false',
                           help='Disable regulator penalty (allow free exploration of regulator numbers)')
    reg_params.add_argument('-reg_optimal', '--reg_optimal', type=int, default=2,
                           help='Optimal number of regulators (K parameter). Default: 2 (Kauffman principle)')
    reg_params.set_defaults(reg_penalty=True)
    
    # GEP technical parameters
    gep_params = parser.add_argument_group('GEP algorithm parameters')
    gep_params.add_argument('-pop', '--population', type=int, default=100,
                           help='Population size (default: 100)')
    gep_params.add_argument('-gen', '--generations', type=int, default=100,
                           help='Number of generations (default: 100)')
    gep_params.add_argument('-elites', '--elites', type=int, default=10,
                           help='Number of elite individuals (default: 10)')
    gep_params.add_argument('-ind_head', '--individual_head', type=int, default=5,
                           help='Head length of individual chromosomes (default: 5)')
    gep_params.add_argument('-ind_genes', '--individual_genes', type=int, default=1,
                           help='Number of genes per individual chromosome (default: 1)')
    gep_params.add_argument('-weights', '--fitness_weights', type=str, default="1,-1,-1",
                           help='Fitness weights (comma-separated). Default: "1,-1,-1" for (correct, reg_penalty, mlp_loss)')
    
    # MLP parameters (optional manual override)
    mlp_params = parser.add_argument_group('MLP validation parameters (optional manual override)')
    mlp_params.add_argument('-mlp_layers', '--mlp_layers', type=str, default=None,
                           help='MANUAL OVERRIDE: MLP hidden layer sizes (comma-separated, e.g., "20,10"). '
                                'If not specified, MLP architecture is configured AUTOMATICALLY based on '
                                'dataset size and number of regulators.')
    mlp_params.add_argument('-mlp_alpha', '--mlp_alpha', type=float, default=None,
                           help='MANUAL OVERRIDE: MLP regularization parameter alpha (e.g., 0.1). '
                                'If not specified, alpha is configured AUTOMATICALLY based on '
                                'number of regulators (ranges from 0.01 to 1.0).')
    
    parser.add_argument('--version', action='version', version='GEP Boolean Network Inference v2.1 (Configurable Regulator Penalty)')
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not os.path.exists(args.input):
        print(f"ERROR: Input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(args.input_binary):
        print(f"ERROR: Binary input file '{args.input_binary}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Validate reg_optimal
    if args.reg_penalty and args.reg_optimal < 1:
        print(f"ERROR: reg_optimal must be >= 1 (got {args.reg_optimal})", file=sys.stderr)
        sys.exit(1)
    
    # Process experiments
    try:
        run_gep_experiment(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()