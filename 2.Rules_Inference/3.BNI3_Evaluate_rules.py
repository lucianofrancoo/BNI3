#!/usr/bin/env python3
"""
Boolean Rules Evaluator with Attractor Metrics
Evaluates different combinations of Boolean rules using attractor analysis
without generating intermediate files.

Author: Luciano
"""

import argparse
import time
import pandas as pd
import numpy as np
import itertools
import math
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import sys
import os
import multiprocessing as mp
from functools import partial

try:
    from numba import njit as _njit
    _NUMBA_AVAILABLE = True
except ImportError:
    def _njit(func):          # transparent no-op when numba is absent
        return func
    _NUMBA_AVAILABLE = False

_RULE_CACHE: Dict = None      # set per-process by _worker_init; None falls back to eval()
_PRIORITY_REGULATORS: set = None   # set per-process by _worker_init; None disables the metric


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
    
    # Edges landing on a prioritized regulator (e.g. known transcription factors).
    # Counted as a preference, never as a requirement: rules without them stay in the
    # running, they just lose this tie-break. See CASCADE_CRITERIA for where it ranks.
    n_priority_regulators = 0
    if _PRIORITY_REGULATORS:
        for rule in gene_rules.values():
            if rule in ['True', 'False', 'true', 'false']:
                continue
            n_priority_regulators += len(
                extract_regulators(rule) & _PRIORITY_REGULATORS)

    return {
        'total_literals': total_literals,
        'total_nots': total_nots,
        'avg_k': avg_k,
        'std_k': std_k,
        'k_distance_from_2': k_distance_from_2,
        'n_priority_regulators': n_priority_regulators
    }


def _precompute_rule_cache(rules_by_gene: Dict[str, List], genes: List[str]) -> Dict:
    """
    Evaluate every unique (gene, rule) pair once across all 2^N states.

    Returns {(gene, rule_str): np.ndarray(2^N, bool)} so _build_transition_table
    can skip eval() entirely for the whole run.
    Called once in the main process; forwarded to each worker via pool initializer.
    """
    N = len(genes)
    total = 1 << N
    states = np.arange(total, dtype=np.int64)
    gene_vals = {
        gene: ((states >> i) & 1).astype(np.bool_)
        for i, gene in enumerate(genes)
    }

    cache: Dict[tuple, np.ndarray] = {}
    for gene in genes:
        for _pos, rule in rules_by_gene[gene]:
            rule_stripped = rule.strip()
            key = (gene, rule_stripped)
            if key in cache:
                continue
            if rule_stripped in ('True', 'true'):
                cache[key] = np.ones(total, dtype=np.bool_)
            elif rule_stripped in ('False', 'false'):
                cache[key] = np.zeros(total, dtype=np.bool_)
            else:
                result = eval(rule_stripped, {'__builtins__': {}}, gene_vals)
                result = np.asarray(result, dtype=np.bool_)
                if result.shape == ():
                    result = np.full(total, bool(result), dtype=np.bool_)
                cache[key] = result
    return cache


def _build_transition_table(gene_rules: Dict[str, str], genes: List[str],
                            rule_cache: Dict = None) -> np.ndarray:
    """
    Build full Boolean network transition table for all 2^N states at once.

    When rule_cache is provided (normal case after precomputation), each gene's
    bool array is a direct dict lookup — no eval() call needed.
    gene_vals is built lazily only on a cache miss (fallback path).

    Returns array T of shape (2^N,) where T[s] = next state integer.
    """
    N = len(genes)
    total = 1 << N

    gene_vals = None  # lazily initialised on first cache miss

    next_bits = []
    for gene in genes:
        rule = gene_rules.get(gene, gene)
        rule_stripped = rule.strip()

        # Fast path: use precomputed bool array
        if rule_cache is not None:
            key = (gene, rule_stripped)
            if key in rule_cache:
                next_bits.append(rule_cache[key])
                continue

        # Fallback: evaluate on the fly (cache miss or no cache)
        if rule_stripped in ('True', 'true'):
            bit = np.ones(total, dtype=np.bool_)
        elif rule_stripped in ('False', 'false'):
            bit = np.zeros(total, dtype=np.bool_)
        else:
            if gene_vals is None:
                states = np.arange(total, dtype=np.int64)
                gene_vals = {
                    g: ((states >> i) & 1).astype(np.bool_)
                    for i, g in enumerate(genes)
                }
            result = eval(rule_stripped, {'__builtins__': {}}, gene_vals)
            result = np.asarray(result, dtype=np.bool_)
            if result.shape == ():
                result = np.full(total, bool(result), dtype=np.bool_)
            bit = result

        next_bits.append(bit)

    # Assemble next-state integers: next_state = Σ next_bits[i] * 2^i
    bit_matrix = np.stack(next_bits, axis=1)
    powers = np.int64(1) << np.arange(N, dtype=np.int64)
    T = bit_matrix.astype(np.int64) @ powers
    return T.astype(np.int32)


@_njit
def _assign_attractors_jit(T):
    """
    Assign every state to an attractor using fixed-size array buffers instead of
    Python dicts — allows Numba to compile this to native code.

    Returns:
        state_attractor : int32 array, state_attractor[s] = attractor index
        in_cycle        : bool array, True for states that ARE the attractor cycle
        n_attractors    : number of distinct attractors found
    """
    N = len(T)
    state_attractor = np.full(N, -1, dtype=np.int32)
    in_cycle        = np.zeros(N, dtype=np.bool_)
    path            = np.empty(N, dtype=np.int32)
    visited_at      = np.full(N, -1, dtype=np.int32)  # step index in current path
    n_attractors    = 0

    for start in range(N):
        if state_attractor[start] >= 0:
            continue

        path_len = 0
        s = start

        while state_attractor[s] < 0 and visited_at[s] < 0:
            visited_at[s] = path_len
            path[path_len] = s
            path_len += 1
            s = T[s]

        if state_attractor[s] >= 0:
            aid = state_attractor[s]
        else:
            # s was already visited in this traversal: new cycle
            aid = n_attractors
            n_attractors += 1
            cycle_begin = visited_at[s]
            for i in range(cycle_begin, path_len):
                in_cycle[path[i]] = True

        for i in range(path_len):
            state_attractor[path[i]] = aid
            visited_at[path[i]] = -1   # reset for future traversals

    return state_attractor, in_cycle, n_attractors


def _find_attractors_from_table(T: np.ndarray, N_genes: int) -> Tuple[List, Dict]:
    """
    Find all attractors given a precomputed transition table.

    Delegates the core O(2^N) traversal to _assign_attractors_jit (Numba-compiled
    when available), then extracts cycle states and basin sizes in Python.
    """
    state_attractor, in_cycle, n_attractors = _assign_attractors_jit(T)

    n_attr = int(n_attractors)
    attractors: List[List[List[bool]]] = [[] for _ in range(n_attr)]
    for s in range(len(T)):
        if in_cycle[s]:
            aid = int(state_attractor[s])
            attractors[aid].append([(s >> i) & 1 == 1 for i in range(N_genes)])

    unique, counts = np.unique(state_attractor, return_counts=True)
    basins = {int(u): int(c) for u, c in zip(unique, counts)}
    return attractors, basins


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
    
    # Metric 4: Attractor concordance with final timepoint (if matrix provided)
    # Every attractor participates, not only fixed points: for a cyclic attractor the
    # best-matching state of the cycle is used. Restricting this to fixed points scored
    # every cyclic network at 0 by construction rather than by disagreeing with the data.
    metrics['final_state_concordance'] = 0
    metrics['concordance_basin_frac'] = 0.0
    if binarized_matrix is not None:
        last_timepoint = binarized_matrix.iloc[-1]

        for idx, attractor in enumerate(attractors):
            best_in_attractor = 0.0
            for state in attractor:
                matches = sum(1 for i, gene in enumerate(genes)
                              if gene in last_timepoint.index and
                              state[i] == bool(last_timepoint[gene]))
                best_in_attractor = max(best_in_attractor, matches / len(genes))

            if best_in_attractor > metrics['final_state_concordance']:
                metrics['final_state_concordance'] = best_in_attractor
                # Share of state space draining into the attractor that matched.
                # Informative only — a match in a tiny basin is weaker evidence than
                # the same match in the dominant one.
                metrics['concordance_basin_frac'] = (
                    basins[idx] / total_states if total_states > 0 else 0.0
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


# ── Lexicographic tie-break cascade ─────────────────────────────────────────────────
# Criteria are compared IN ORDER: the first one that differs decides, and the remaining
# ones are never consulted. This replaces the weighted sum of normalized metrics, where
# a small parsimony gain could outweigh a real loss of agreement with the data.
#
# Each entry is (metric name, sign); sign = +1 when lower is better, -1 when higher is
# better. Metrics are used raw — no normalization — so keys are comparable across runs.
CASCADE_CRITERIA = [
    ('final_state_concordance', -1),   # 1. agreement with the observed final state
    ('n_basins',                +1),   # 2. fewer basins
    ('avg_cycle_length',        +1),   # 3. shorter attractors
    ('fixed_point_ratio',       -1),   # 4. more fixed points
    ('basin_entropy',           +1),   # 5. one dominant basin
    ('k_distance_from_2',       +1),   # 6. Kauffman K ≈ 2
    ('total_literals',          +1),   # 7. simpler rules
    ('total_nots',              +1),   # 8. less repression
]

# Decimals kept when building a key. Guards against float noise (1e-16 differences)
# splitting combinations that are in fact tied.
_CASCADE_ROUND = 9

# Where the optional priority-regulator preference is inserted when
# --priority-regulators is given: after the attractor block, before parsimony.
# It therefore never overrides how well a network fits the data or how clean its
# dynamics are — it only chooses among networks already tied on those, preferring
# the one built on known regulators before preferring the simplest one. Move the
# index to 8 to make it a last-resort tie-break instead.
_PRIORITY_CASCADE_POSITION = 5
_PRIORITY_CRITERION = ('n_priority_regulators', -1)

# The criteria actually in force this run; set by configure_cascade().
_ACTIVE_CASCADE = list(CASCADE_CRITERIA)


def configure_cascade(priority_regulators: set) -> List[tuple]:
    """Activate the priority-regulator criterion when a regulator list was given."""
    global _ACTIVE_CASCADE
    _ACTIVE_CASCADE = list(CASCADE_CRITERIA)
    if priority_regulators:
        _ACTIVE_CASCADE.insert(_PRIORITY_CASCADE_POSITION, _PRIORITY_CRITERION)
    return _ACTIVE_CASCADE


def load_priority_regulators(value: str) -> set:
    """
    Parse the prioritized regulators given on the command line.

    The normal form is a comma-separated list of gene names, so the regulators can
    travel in the pipeline command itself:

        --priority-regulators HB6,MYB44,ABF3

    A path to a file with one name per line is also accepted, for lists too long to
    type comfortably; blank lines and lines starting with # are ignored there.

    Names are sanitized like everywhere else, so they match the identifiers used
    inside rules whether they are written "SnRK2.8" or "SnRK2_8".
    """
    if not value:
        return None

    if ',' not in value and os.path.isfile(value):
        with open(value, encoding='utf-8') as handle:
            names = [line.strip() for line in handle]
        names = [n for n in names if n and not n.startswith('#')]
        source = f"file {value}"
    else:
        names = [n.strip() for n in value.split(',') if n.strip()]
        source = "--priority-regulators"

    regulators = {sanitize_gene_name(n) for n in names}

    if not regulators:
        print(f"ERROR: no gene names found in {source}", file=sys.stderr)
        sys.exit(1)

    return regulators


def cascade_key(metrics: Dict) -> tuple:
    """
    Build the lexicographic ranking key for one combination.

    Smaller tuple = better combination, so the key can be used directly with min(),
    sorted() or pandas ascending sorts.
    """
    return tuple(round(sign * float(metrics[name]), _CASCADE_ROUND)
                 for name, sign in _ACTIVE_CASCADE)


def load_rules_table(rules_file: str) -> pd.DataFrame:
    """Load rules table from TSV file"""
    df = pd.read_csv(rules_file, sep='\t', encoding='utf-8')
    required_cols = ['Gene', 'Position', 'Rule']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    return df


# Same expression the inference script uses (load_and_clean_data): any character
# that is not a letter, digit or underscore becomes "_", and a leading digit gets an
# "_" prefix. Gene names must be valid Python identifiers because rules are evaluated
# with eval(). "SnRK2.8" -> "SnRK2_8".
_GENE_NAME_RE = r'\W|^(?=\d)'


def sanitize_gene_name(name) -> str:
    """Turn one gene name into a valid Python identifier (idempotent)."""
    return re.sub(_GENE_NAME_RE, '_', str(name))


def sanitize_gene_names(names, source: str) -> List[str]:
    """
    Sanitize a list of gene names, refusing silently ambiguous results.

    Sanitization is many-to-one: "SnRK2.8", "SnRK2 8" and "SnRK2-8" all collapse to
    "SnRK2_8", and so do "GEN.1" and "GEN_1". Two distinct genes mapping to the same
    identifier would make every later lookup ambiguous, so that is an error rather
    than a warning.
    """
    clean = [sanitize_gene_name(n) for n in names]

    groups: Dict[str, set] = {}
    for original, cleaned in zip(names, clean):
        groups.setdefault(cleaned, set()).add(str(original))
    collisions = {k: v for k, v in groups.items() if len(v) > 1}

    if collisions:
        print(f"\nERROR: gene names in {source} become ambiguous once sanitized to "
              f"Python identifiers.", file=sys.stderr)
        for cleaned, originals in list(collisions.items())[:10]:
            print(f"  {sorted(originals)} all become '{cleaned}'", file=sys.stderr)
        print(f"  Rename them upstream so they stay distinct.", file=sys.stderr)
        sys.exit(1)

    return clean


def load_binarized_matrix(matrix_file: str) -> pd.DataFrame:
    """
    Load a binarized expression matrix: rows are timepoints, columns are genes.

    An optional first column holding timepoint or sample labels is moved to the
    index. Whether it is such a label column is decided from its VALUES, not from
    its name: a real gene column contains only 0/1. Guessing from the name (the
    previous heuristic accepted only names starting with G/AT/LOC/ENSG) silently
    swallowed the first gene of any matrix using plain symbols like "RCAR1".
    """
    df = pd.read_csv(matrix_file, sep='\t')

    # Match the inference script's sanitization, otherwise a matrix column
    # "AT1G01360 (RCAR1)" never matches the rule gene "AT1G01360__RCAR1_" and
    # final_state_concordance silently scores 0 for every combination.
    df.columns = sanitize_gene_names(df.columns, os.path.basename(matrix_file))

    def is_binary(series) -> bool:
        values = pd.to_numeric(series, errors='coerce')
        if values.isna().any():
            return False
        return set(values.unique()) <= {0, 1}

    if len(df.columns) > 1 and not is_binary(df.iloc[:, 0]):
        df = df.set_index(df.columns[0])

    non_binary = [c for c in df.columns if not is_binary(df[c])]
    if non_binary:
        print(f"\nERROR: {len(non_binary)} column(s) of {os.path.basename(matrix_file)} "
              f"are not binary (0/1): {non_binary[:5]}"
              f"{' ...' if len(non_binary) > 5 else ''}", file=sys.stderr)
        print(f"  This file should be the binarized matrix, not the expression matrix.",
              file=sys.stderr)
        sys.exit(1)

    return df


def _rule_value_at_state(rule: str, gene_state: Dict[str, bool]) -> bool:
    """Evaluate one rule at a single network state. Raises on unknown regulators."""
    stripped = rule.strip()
    if stripped in ('True', 'true'):
        return True
    if stripped in ('False', 'false'):
        return False
    return bool(eval(stripped, {'__builtins__': {}}, gene_state))


def prune_rules_by_final_state(rules_by_gene: Dict[str, List],
                               genes: List[str],
                               binarized_matrix: pd.DataFrame,
                               verbose: bool = True) -> Tuple[Dict[str, List], Dict]:
    """
    Drop candidate rules that cannot reproduce the observed final state.

    This is an exact prune of the first cascade criterion, not a heuristic.
    final_state_concordance reaches its maximum of 1.0 exactly when some attractor
    state equals the observed final state s*. The simplest way that happens is for
    s* to be a fixed point, and THAT condition factorizes over genes: s* is a fixed
    point iff rule_g(s*) == s*_g for every gene g independently.

    So the subset of combinations achieving concordance 1.0 through a fixed point can
    be found with sum(n_g) rule evaluations instead of enumerating prod(n_g)
    combinations. Since concordance ranks first in CASCADE_CRITERIA and 1.0 is its
    ceiling, every combination outside this subset is dominated and can never win.

    The prune is skipped when any gene keeps no rule: s* is then not a fixed point of
    any reachable combination, concordance < 1.0 everywhere, and the discarded rules
    would still be in contention.

    Returns:
        (rules_by_gene, info) — info carries the before/after sizes for reporting.
    """
    last_timepoint = binarized_matrix.iloc[-1]
    missing = [g for g in genes if g not in last_timepoint.index]
    if missing:
        return rules_by_gene, {'applied': False,
                               'reason': f'{len(missing)} gene(s) absent from the matrix'}

    target = {g: bool(last_timepoint[g]) for g in genes}
    gene_state = dict(target)

    kept: Dict[str, List] = {}
    for gene in genes:
        survivors = []
        for position, rule in rules_by_gene[gene]:
            try:
                if _rule_value_at_state(rule, gene_state) == target[gene]:
                    survivors.append((position, rule))
            except Exception:
                # Unknown regulator or malformed rule: keep it rather than prune
                # something the full evaluation might still be able to score.
                survivors.append((position, rule))
        kept[gene] = survivors

    empty = [g for g in genes if not kept[g]]
    if empty:
        return rules_by_gene, {
            'applied': False,
            'reason': (f'{len(empty)} gene(s) have no rule reproducing the final state, '
                       f'so it is not a fixed point of any combination')
        }

    before = math.prod(len(rules_by_gene[g]) for g in genes)
    after = math.prod(len(kept[g]) for g in genes)
    info = {'applied': True, 'before': before, 'after': after,
            'evaluations': sum(len(rules_by_gene[g]) for g in genes),
            'per_gene': {g: (len(rules_by_gene[g]), len(kept[g])) for g in genes}}

    if verbose:
        print(f"   Pruning to rules consistent with the observed final state...")
        for gene in genes:
            n_before, n_after = info['per_gene'][gene]
            print(f"   - {gene}: {n_before} -> {n_after}")
        # A fixed 2-decimal percentage reads "0.00%" for any serious reduction, which
        # says nothing. Report the reduction factor, which stays readable at any scale.
        kept_total = sum(len(kept[g]) for g in genes)
        print(f"   Candidate rules: {info['evaluations']:,} -> {kept_total:,} "
              f"(one evaluation each, at the observed final state)")
        print(f"   Search space: {before:.4e} -> {after:.4e}, "
              f"{before / after:,.0f}x smaller")
        print(f"   Every surviving combination scores the maximum 1.00 on "
              f"final_state_concordance")

    return kept, info


def calculate_intelligent_defaults(rules_by_gene: Dict[str, List],
                                  target_time_minutes: int = 10,
                                  n_processes: int = None) -> Dict[str, int]:
    """
    Calculate intelligent default values for top_n and max_combinations
    based on the problem size and desired execution time.

    Args:
        rules_by_gene: Dictionary with rules per gene
        target_time_minutes: Target execution time in minutes (default: 10)
        n_processes: Number of processes (default: None = auto-detect all available cores)

    Returns:
        Dictionary with recommended 'top_n' and 'max_combinations'
    """
    # Get number of rules per gene
    rules_counts = [len(rules) for rules in rules_by_gene.values()]
    n_genes = len(rules_by_gene)

    # ===== LÍMITES =====
    MAX_COMBINATIONS_ABSOLUTE = 50_000_000  # Safety cap (50M)
    MAX_RULES_PER_GENE_ABSOLUTE = 30        # Máximo de reglas por gen
    TARGET_COMBINATIONS_FOR_SAMPLING = 5_000_000  # Target cuando se necesita muestreo

    # Calculate total possible combinations (with overflow protection)
    try:
        total_combinations = math.prod(rules_counts)
    except Exception:
        total_combinations = float('inf')

    # Estimate throughput (combinations per second per core).
    # The flat 300/core this used to assume ignored network size, so the printed time
    # estimate was off by ~4x on a 15-gene network (predicted 4.3 min, took 17). Cost
    # per combination is dominated by building the transition table over all 2^N states
    # for N genes, so throughput scales as 1/(N * 2^N). The constant is calibrated from
    # a measured run (15 genes, 64 processes, 4,865 combos/s) and predicts an
    # independent 13-gene run to within 1%.
    _THROUGHPUT_CONSTANT = 37_000_000
    throughput_per_core = max(1.0, _THROUGHPUT_CONSTANT / (n_genes * (2 ** n_genes))) \
        if n_genes < 40 else 1.0
    if n_processes is None:
        n_processes = mp.cpu_count()
    total_throughput = throughput_per_core * n_processes
    
    # Calculate how many combinations we can evaluate in target time
    max_feasible = int(total_throughput * target_time_minutes * 60)
    
    # Apply absolute limit (no more than 50K regardless of time)
    max_feasible = min(max_feasible, MAX_COMBINATIONS_ABSOLUTE)
    
    # Strategy 1: total fits within time budget → evaluate ALL
    if total_combinations <= max_feasible:
        return {
            'top_n': None,
            'max_combinations': None,
            'reason': 'small_problem',
            'estimated_time_minutes': total_combinations / total_throughput / 60 if total_combinations != float('inf') else target_time_minutes
        }

    # Strategy 2: too many to evaluate all → sample, but keep ALL max-score rules
    # top_n is never reduced: all rules with maximum score per gene are always kept.
    recommended_combos = min(TARGET_COMBINATIONS_FOR_SAMPLING, int(total_combinations) if total_combinations != float('inf') else TARGET_COMBINATIONS_FOR_SAMPLING)
    return {
        'top_n': None,
        'max_combinations': int(recommended_combos),
        'reason': 'sample_combinations',
        'estimated_time_minutes': recommended_combos / total_throughput / 60
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


def _worker_init(rule_cache: Dict, priority_regulators: set = None) -> None:
    """Set per-worker globals: rule cache, priority regulators, Numba JIT warm-up."""
    global _RULE_CACHE, _PRIORITY_REGULATORS
    _RULE_CACHE = rule_cache
    _PRIORITY_REGULATORS = priority_regulators
    if _NUMBA_AVAILABLE:
        _assign_attractors_jit(np.array([0, 0], dtype=np.int32))




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

    # Build transition table using the precomputed rule cache (zero eval() calls
    # when the cache covers all rules in this combination).
    T = _build_transition_table(gene_rules, genes, rule_cache=_RULE_CACHE)
    attractors, basins = _find_attractors_from_table(T, len(genes))

    # Calculate attractor-based metrics
    metrics = calculate_attractor_metrics(attractors, basins, genes, binarized_matrix)

    # Calculate parsimony-based metrics
    parsimony_metrics = calculate_parsimony_metrics(gene_rules)
    metrics.update(parsimony_metrics)

    # Add position information to metrics
    metrics['positions'] = [rules_by_gene[gene][combo_indices[i]][0]
                            for i, gene in enumerate(genes)]

    return combo_indices, metrics


def write_evaluation_log(path: str, header: Dict, trace: List[Dict]) -> None:
    """
    Write a compact trace of how the search progressed.

    One row per restart (hill climbing) or per batch (sampling), plus provenance as
    '#' comment lines so the file stays a single self-describing artifact that
    pandas reads with comment='#'. The point is to be able to reconstruct afterwards
    how the answer was reached — when the best key improved, how ambiguity evolved,
    how the coverage bound tightened — without keeping millions of rows.
    """
    if not trace:
        return

    columns = list(trace[0].keys())

    with open(path, 'w', encoding='utf-8') as handle:
        for key, value in header.items():
            handle.write(f"# {key}: {value}\n")
        handle.write('\t'.join(columns) + '\n')
        for row in trace:
            handle.write('\t'.join(str(row.get(c, '')) for c in columns) + '\n')


def _hill_climb_search(genes: List[str],
                       rules_by_gene: Dict[str, List],
                       eval_func,
                       rule_cache: Dict,
                       priority_regulators: set,
                       n_processes: int,
                       budget: int,
                       restarts: int = None,
                       target_restarts: int = None,
                       coverage_conf: float = 0.99,
                       seed: int = 42,
                       verbose: bool = True) -> Tuple[list, Dict]:
    """
    Random-restart steepest-ascent hill climbing over the combination space.

    Uniform sampling scales terribly here: drawing 15,000,000 points out of 10^27
    independently explores nothing, and typically ends with k=1 — a single isolated
    best, not a reproducible optimum. Hill climbing spends the same budget walking
    toward optima instead.

    One sweep evaluates every single-gene change of the current combination — about
    sum(n_g) candidates, so ~1,100 for a 15-gene network rather than prod(n_g) — and
    moves to the best one. When no single-gene change improves the cascade key the
    combination is a local optimum, and the search restarts from a fresh random point.

    Evaluating all genes' alternatives as one batch (steepest ascent) rather than gene
    by gene also keeps the worker pool busy with a reasonably sized chunk of work.

    Returns:
        (results, info) — results is the list of (combo, metrics) actually evaluated,
        deduplicated; info carries restart statistics.
    """
    # Seed the starting points. Without this, hill climbing drew a different set of
    # restarts on every invocation, so the same inputs gave different restart counts,
    # different convergence rates and potentially a different winning network — the
    # sampling path was seeded but this one never was.
    if seed is not None:
        np.random.seed(seed)

    n_options = [len(rules_by_gene[g]) for g in genes]
    n_genes = len(genes)
    seen: Dict[tuple, Dict] = {}

    cap = budget if budget else float('inf')

    def evaluate(combos, pool):
        """Evaluate the combinations not seen yet, respecting the budget."""
        fresh = []
        for combo in combos:
            if combo not in seen and len(seen) + len(fresh) < cap:
                fresh.append(combo)
        if not fresh:
            return
        if pool is None:
            for combo in fresh:
                seen[combo] = eval_func(combo)[1]
        else:
            chunksize = max(1, len(fresh) // (n_processes * 4))
            for combo, metrics in pool.imap_unordered(eval_func, fresh,
                                                      chunksize=chunksize):
                seen[tuple(combo)] = metrics

    def neighbours(combo):
        """Every combination one rule swap away from this one."""
        out = []
        for i in range(n_genes):
            for j in range(n_options[i]):
                if j != combo[i]:
                    candidate = list(combo)
                    candidate[i] = j
                    out.append(tuple(candidate))
        return out

    pool = None
    best_key = None
    restart_keys: list = []

    # Combinations observed at the current best key. Their count is the ambiguity left
    # in the answer, and per gene the number of distinct rules among them says which
    # genes the cascade has settled. Both reset when a strictly better optimum appears
    # and grow as ties accumulate — a sawtooth, not a monotone shrink. The search space
    # itself never shrinks during the climb; the prune already did that once, up front.
    #
    # The count is of combinations actually seen. Multiplying the per-gene rule counts
    # instead would give the Cartesian closure of those ties, which includes
    # combinations never evaluated and overstates the ambiguity by many orders of
    # magnitude.
    best_combos: set = set()
    trace: List[Dict] = []

    def note(combo, key):
        """Track the combinations tied at the best key found so far."""
        nonlocal best_key
        if best_key is None or key < best_key:
            best_key = key
            best_combos.clear()
            best_combos.add(combo)
        elif key == best_key:
            best_combos.add(combo)

    try:
        if n_processes > 1:
            pool = mp.Pool(processes=n_processes, initializer=_worker_init,
                           initargs=(rule_cache, priority_regulators))

        # Restarts, not evaluations, are the statistically meaningful unit here: the
        # coverage guarantee is (1-x)^R, independent of how many combinations each
        # restart happened to touch. So target_restarts is the stopping criterion and
        # budget is only a safety cap.
        restart_num = 0
        def done():
            if restarts is not None:
                return restart_num >= restarts
            if target_restarts is not None:
                return restart_num >= target_restarts
            return False

        while len(seen) < cap and not done():
            restart_num += 1
            key_before_restart = best_key
            current = tuple(np.random.randint(0, n) for n in n_options)
            evaluate([current], pool)
            if current not in seen:
                break                      # budget exhausted mid-restart
            current_key = cascade_key(seen[current])
            note(current, current_key)

            while len(seen) < cap:
                candidates = neighbours(current)
                evaluate(candidates, pool)

                scored = [(cascade_key(seen[c]), c) for c in candidates if c in seen]
                if not scored:
                    break
                for key, combo in scored:
                    note(combo, key)
                step_key, step_combo = min(scored)
                if step_key >= current_key:
                    break                  # local optimum: no single swap improves it
                current, current_key = step_combo, step_key

            restart_keys.append(current_key)

            fixed = sum(1 for i in range(n_genes)
                        if len({c[i] for c in best_combos}) == 1)
            # Restarts are uniform draws of a starting point, so after R of them any
            # optimum still unseen must sit in a basin smaller than this. The
            # confidence is the same one the stopping criterion uses, so the live
            # number converges on the guarantee the run promises instead of
            # describing a different one.
            basin = 1 - (1 - coverage_conf) ** (1 / restart_num)

            row = {'restart': restart_num,
                   'evaluations_total': len(seen),
                   'improved': int(key_before_restart is None or
                                   best_key < key_before_restart),
                   'genes_fixed': fixed,
                   'combos_tied_at_best': len(best_combos),
                   'missed_basin': round(basin, 6),
                   'basin_confidence': coverage_conf}
            # The best key so far, as readable metric values (key holds sign*value)
            for (name, sign), component in zip(_ACTIVE_CASCADE, best_key):
                row[f'best_{name}'] = round(sign * component, 6)
            trace.append(row)

            if verbose:
                progress = (f"Restart {restart_num}/{target_restarts}"
                            if target_restarts else f"Restart {restart_num}")
                evals = (f"{len(seen):,}/{budget:,} evals" if budget
                         else f"{len(seen):,} evals")
                print(f"   {progress} | {evals} "
                      f"| {fixed}/{n_genes} genes fixed | {len(best_combos):,} tied "
                      f"| missed basin <{basin:.3%} @{coverage_conf:.1%}   ", end='\r')
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    at_best = sum(1 for k in restart_keys if k == best_key)
    reached_target = (target_restarts is not None
                      and len(restart_keys) >= target_restarts)
    info = {'restarts': len(restart_keys), 'restarts_at_best': at_best,
            'evaluated': len(seen), 'trace': trace,
            'target_restarts': target_restarts,
            'stopped_on': ('restart target' if reached_target
                           else 'evaluation budget')}

    if verbose:
        print(f"\n   Hill climbing finished — {len(seen):,} combinations evaluated "
              f"in {len(restart_keys)} restarts ({info['stopped_on']} reached)")
        print(f"   {at_best} of {len(restart_keys)} restarts reached the best key "
              f"({at_best / max(1, len(restart_keys)):.1%})")

    return [(list(combo), metrics) for combo, metrics in seen.items()], info


def evaluate_rule_combinations(rules_file: str,
                               output_file: str = None,
                               binarized_matrix_file: str = None,
                               top_n: int = None,
                               max_combinations: int = None,
                               max_iterations: int = 1000,
                               n_processes: int = None,
                               top_results: int = 10000,
                               score_column: str = 'Score',
                               confidence_target: float = 0.999,
                               prune_concordance: bool = True,
                               search: str = 'auto',
                               restarts: int = None,
                               target_basin: float = 0.01,
                               seed: int = 42,
                               priority_regulators: set = None,
                               write_log: bool = True,
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
        n_processes: Number of processes for parallelization (default: None = auto-detect all available cores)
        top_results: Number of top-ranked combinations to write to evaluation_results.tsv.
                     If None, writes all evaluated combinations (can be very large).
        score_column: Column to use for rule ranking
        confidence_target: Target P(found global optimum). Adaptive batching adds batches
                           until k ≥ ceil(-ln(1-p)) combinations share the best cascade key
                           are observed (confidence ≈ 1−e^{−k}). Default: 0.99 (k≥5).
                           Set to 0 to disable adaptive batching (single batch only).
        verbose: Print progress information
    """
    # Resolve process count before any downstream use
    if n_processes is None:
        n_processes = mp.cpu_count()

    t_start = time.perf_counter()

    if verbose:
        print("="*80)
        print("BOOLEAN RULES EVALUATION WITH INTEGRATED SCORE")
        print("="*80)

    # Load data
    if verbose:
        print(f"\n1. Loading rules from: {rules_file}")
    df = load_rules_table(rules_file)

    # Activate the optional priority-regulator criterion before anything builds a key
    configure_cascade(priority_regulators)
    if priority_regulators and verbose:
        print(f"   Prioritizing {len(priority_regulators)} regulators as cascade "
              f"criterion #{_PRIORITY_CASCADE_POSITION + 1} (a preference, not a filter)")
    
    # Generate output file path if not provided
    if output_file is None:
        rules_dir = os.path.dirname(rules_file)
        output_file = os.path.join(rules_dir, "evaluation_results.tsv")
        if verbose:
            print(f"   Output will be saved to: {output_file}")
    
    # Whether the evaluation cap was asked for, as opposed to filled in by the
    # intelligent defaults. In hill mode the restart target is the criterion, so an
    # automatic cap must not be allowed to cut the run short of its guarantee.
    max_combos_explicit = max_combinations is not None

    binarized_matrix = None
    if binarized_matrix_file:
        if verbose:
            print(f"2. Loading binarized matrix from: {binarized_matrix_file}")
        binarized_matrix = load_binarized_matrix(binarized_matrix_file)

        # The concordance criterion compares attractor states against this matrix
        # by gene name. A name mismatch makes every comparison miss, scoring 0 for
        # every combination — a silent failure that invalidates the top criterion
        # of the ranking cascade without any visible error. Check it up front.
        rule_genes = set(df['Gene'].unique())
        matrix_genes = set(binarized_matrix.columns)
        overlap = rule_genes & matrix_genes

        if not overlap:
            print(f"\nERROR: none of the {len(rule_genes)} genes in the rules table "
                  f"match a column of the binarized matrix.", file=sys.stderr)
            print(f"  rules:  {sorted(rule_genes)[:3]}", file=sys.stderr)
            print(f"  matrix: {sorted(matrix_genes)[:3]}", file=sys.stderr)
            print(f"  final_state_concordance would be 0.00 for every combination, "
                  f"silently disabling the first cascade criterion.", file=sys.stderr)
            print(f"  Check that the matrix is the one used for the inference run.",
                  file=sys.stderr)
            sys.exit(1)

        if len(overlap) < len(rule_genes):
            missing = sorted(rule_genes - matrix_genes)
            print(f"\nWARNING: {len(missing)} of {len(rule_genes)} rule genes have no "
                  f"column in the binarized matrix: {missing[:5]}"
                  f"{' ...' if len(missing) > 5 else ''}")
            print(f"  They cannot contribute to final_state_concordance, capping it "
                  f"at {len(overlap)/len(rule_genes):.2f} instead of 1.00.")
        elif verbose:
            print(f"   All {len(overlap)} rule genes matched to matrix columns")

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
            strategy_labels = {
                'small_problem':      'Evaluate ALL combinations (fits within time budget)',
                'sample_combinations':'Sample combinations — all max-score rules kept, too many to evaluate all',
            }
            label = strategy_labels.get(defaults['reason'], defaults['reason'])
            print(f"   Problem size analysis:")
            print(f"   - Strategy: {label}")
            if top_n is not None:
                print(f"   - Recommended top_n: {top_n} rules per gene")
            else:
                print(f"   - Using ALL rules with maximum score per gene")
            if max_combinations is not None:
                print(f"   - Recommended max_combinations: {max_combinations:,}")
            else:
                print(f"   - Will evaluate ALL possible combinations")
            # In hill mode the budget is only a cap, so this is an upper bound the run
            # usually stops well short of — the restart target normally fires first.
            bound = (" at the evaluation cap" if search in ('auto', 'hill')
                     and defaults['max_combinations'] else "")
            print(f"   - Estimated time: ~{defaults['estimated_time_minutes']:.1f} "
                  f"minutes{bound}")
    
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

    # Exact prune on the first cascade criterion — see prune_rules_by_final_state
    prune_info = {'applied': False, 'reason': 'no binarized matrix given'}
    if binarized_matrix is not None and prune_concordance:
        rules_by_gene, prune_info = prune_rules_by_final_state(
            rules_by_gene, genes, binarized_matrix, verbose)
        if verbose and not prune_info['applied']:
            print(f"   Prune skipped: {prune_info['reason']}")

    # Calculate total possible combinations
    total_combinations = math.prod([len(rules_by_gene[gene]) for gene in genes])
    
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
                how = ('hill climbing' if search in ('auto', 'hill') else 'random sampling')
                print(f"   Too many to enumerate — will explore by {how}")
        else:
            if verbose:
                print(f"   Evaluating all {total_combinations:,} combinations (less than max_combinations)")
    
    # Precompute all unique (gene, rule) → bool array mappings once for the whole run.
    if verbose:
        print(f"\n6. Precomputing rule evaluation cache...")
    rule_cache = _precompute_rule_cache(rules_by_gene, genes)
    if verbose:
        print(f"   {len(rule_cache)} unique rules cached")

    # Set cache in main process for single-process fallback
    global _RULE_CACHE
    _RULE_CACHE = rule_cache

    # Warm up Numba JIT on single-process runs (avoids first-call latency mid-evaluation)
    if _NUMBA_AVAILABLE and n_processes == 1:
        _assign_attractors_jit(np.array([0, 0], dtype=np.int32))

    # Create partial function
    eval_func = partial(evaluate_single_combination,
                        genes=genes,
                        rules_by_gene=rules_by_gene,
                        binarized_matrix=binarized_matrix,
                        max_iterations=max_iterations)

    # ── Adaptive evaluation loop ────────────────────────────────────────────────────
    # Adds batches until k ≥ target_k, where confidence ≈ 1−e^{−k} ≥ confidence_target.
    # k = number of evaluated combos sharing the best cascade key (raw metrics, stable
    # across batches). The n and M cancel in the coverage formula, leaving only k.
    is_exhaustive       = (n_combos_to_eval >= total_combinations)
    adaptive_cap        = (min(n_combos_to_eval * 3, 50_000_000)
                           if not is_exhaustive else total_combinations)
    target_k            = math.ceil(-math.log(1.0 - min(confidence_target, 0.9999)))
    all_results:        list  = []
    sample_trace:       list  = []
    batch_num:          int   = 0
    k_converged:        int   = 0
    confidence_achieved: float = 0.0

    # Hill climbing replaces uniform sampling when the space cannot be enumerated
    # 'auto': enumerate everything when it fits, otherwise climb. Uniform sampling
    # is never the automatic choice — on spaces large enough to need sampling, hill
    # climbing found strictly better networks with a fraction of the budget. Ask for
    # 'sample' explicitly when an unbiased sample is what you want, for instance
    # because you need the 1-e^-k coverage estimate, which assumes uniform draws.
    use_hill  = (not is_exhaustive) and search in ('auto', 'hill')
    hill_info = None

    if verbose and not is_exhaustive:
        chosen = 'hill climbing' if use_hill else 'uniform sampling'
        why = 'automatic' if search == 'auto' else 'requested'
        print(f"\n   Search strategy: {chosen} ({why})")
    elif verbose and search == 'hill':
        print(f"\n   Search strategy: exhaustive — the whole space fits, "
              f"so it beats hill climbing")

    if use_hill:
        # Restarts needed so that, at confidence_target, any optimum still unseen
        # must sit in a basin smaller than target_basin: R = ln(1-conf)/ln(1-x)
        target_restarts = None
        if restarts is None and target_basin and 0 < target_basin < 1:
            conf = min(confidence_target, 0.9999) if confidence_target else 0.99
            target_restarts = math.ceil(math.log(1 - conf) / math.log(1 - target_basin))

        # The restart target defines when the run is done; a cap the defaults filled in
        # would only stop it short of the guarantee it just promised, so it is dropped
        # unless the user asked for one with --max-combos.
        hill_budget = n_combos_to_eval
        if target_restarts and not max_combos_explicit:
            hill_budget = None

        if verbose:
            print(f"\n7. Starting hill-climbing search with {n_processes} CPU "
                  f"processes...")
            if target_restarts:
                print(f"   Stopping criterion: {target_restarts:,} restarts — enough "
                      f"to guarantee at {conf:.1%} confidence that any optimum not "
                      f"found")
                print(f"   sits in a basin under {target_basin:.3%} of starting points.")
            if hill_budget:
                print(f"   Evaluation cap: {hill_budget:,} (from --max-combos). The "
                      f"run stops at whichever comes first.")
            elif target_restarts:
                print(f"   No evaluation cap — the restart target ends the run. "
                      f"Pass --max-combos to impose one.")
            else:
                print(f"   Evaluation budget: {n_combos_to_eval:,}")

        all_results, hill_info = _hill_climb_search(
            genes, rules_by_gene, eval_func, rule_cache, priority_regulators,
            n_processes, hill_budget, restarts, target_restarts,
            min(confidence_target, 0.9999) if confidence_target else 0.99,
            seed, verbose)
        batch_num = 1
        best_key = None
        for _, m in all_results:
            key = cascade_key(m)
            if best_key is None or key < best_key:
                best_key, k_converged = key, 1
            elif key == best_key:
                k_converged += 1
        # The coverage formula assumes uniform random sampling; hill climbing is
        # deliberately biased, so 1-e^-k would not mean anything here. Restart
        # convergence is reported instead.
        confidence_achieved = 0.0
    elif verbose:
        if is_exhaustive:
            print(f"\n7. Starting evaluation with {n_processes} CPU processes...")
        else:
            print(f"\n7. Starting adaptive evaluation — "
                  f"target confidence: {confidence_target:.1%} (need k ≥ {target_k})...")

    while not use_hill:
        batch_num += 1

        # ── Generate this batch's combinations ───────────────────────────────────
        if batch_num == 1:
            if is_exhaustive:
                if verbose:
                    print(f"   Generating all {total_combinations:,} combinations...")
                batch_combos = list(itertools.product(
                    *[range(len(rules_by_gene[gene])) for gene in genes]
                ))
            else:
                np.random.seed(42)
                seen: set = set()
                for _ in range(int(n_combos_to_eval * 1.05) + 100):
                    seen.add(tuple(np.random.randint(0, len(rules_by_gene[gene]))
                                   for gene in genes))
                    if len(seen) >= n_combos_to_eval:
                        break
                batch_combos = list(seen)
                if verbose:
                    print(f"   Batch 1: {len(batch_combos):,} unique combinations")
        else:
            n_total_needed = math.ceil(target_k * len(all_results) / k_converged)
            n_additional   = min(n_total_needed - len(all_results),
                                 adaptive_cap - len(all_results))
            if n_additional <= 0:
                break
            if verbose:
                print(f"\n   Batch {batch_num}: {n_additional:,} more combinations "
                      f"(k={k_converged}, confidence so far: {confidence_achieved:.1%})...")
            seen = set()
            for _ in range(int(n_additional * 1.05) + 100):
                seen.add(tuple(np.random.randint(0, len(rules_by_gene[gene]))
                               for gene in genes))
                if len(seen) >= n_additional:
                    break
            batch_combos = list(seen)

        # ── Evaluate this batch ──────────────────────────────────────────────────
        batch_results:  list = []
        offset:         int  = len(all_results)
        print_interval: int  = max(100, len(batch_combos) // 50)
        chunksize = max(1, len(batch_combos) // (n_processes * 5))

        if n_processes > 1:
            with mp.Pool(processes=n_processes,
                         initializer=_worker_init, initargs=(rule_cache, priority_regulators)) as pool:
                for i, result in enumerate(
                    pool.imap_unordered(eval_func, batch_combos, chunksize=chunksize), 1
                ):
                    batch_results.append(result)
                    if verbose and i % print_interval == 0:
                        print(f"   Progress: {offset + i:,} combinations evaluated",
                              end='\r')
        else:
            for i, combo in enumerate(batch_combos, 1):
                batch_results.append(eval_func(combo))
                if verbose and i % print_interval == 0:
                    print(f"   Progress: {offset + i:,} combinations evaluated",
                          end='\r')

        all_results.extend(batch_results)

        # ── Compute k and confidence ─────────────────────────────────────────────
        # k counts the combinations sharing the best cascade key — the very key that
        # picks the winner, so the coverage estimate applies to the ranking actually
        # used. (It used to be counted on composite_score, a different metric: the
        # reported confidence then described a combination that need not be the one
        # the script returned.) Raw metrics → keys are stable across batches.
        best_key    = None
        k_converged = 0
        for _, m in all_results:
            key = cascade_key(m)
            if best_key is None or key < best_key:
                best_key, k_converged = key, 1
            elif key == best_key:
                k_converged += 1
        confidence_achieved = 1.0 if is_exhaustive else 1.0 - math.exp(-k_converged)

        sample_trace.append({
            'batch': batch_num,
            'evaluations_total': len(all_results),
            'k_at_best': k_converged,
            'confidence': round(confidence_achieved, 6),
            **{f'best_{name}': round(sign * component, 6)
               for (name, sign), component in zip(_ACTIVE_CASCADE, best_key)}
        })

        if verbose:
            print(f"\n   Batch {batch_num} complete — "
                  f"{len(all_results):,} total | k={k_converged} | "
                  f"confidence={confidence_achieved:.1%}")

        # ── Stopping conditions ──────────────────────────────────────────────────
        if (is_exhaustive or confidence_achieved >= confidence_target
                or len(all_results) >= adaptive_cap):
            break

    results = all_results
    if verbose:
        batches_label = f" in {batch_num} batch(es)" if batch_num > 1 else ""
        print(f"\n   Total: {len(results):,} combinations evaluated{batches_label}")
    
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
            'concordance_basin_frac': metrics.get('concordance_basin_frac', 0.0),
            'fixed_point_ratio': metrics['fixed_point_ratio'],
            # Parsimony metrics
            'total_literals': metrics['total_literals'],
            'total_nots': metrics['total_nots'],
            'avg_k': metrics['avg_k'],
            'std_k': metrics['std_k'],
            'k_distance_from_2': metrics['k_distance_from_2'],
            'n_priority_regulators': metrics.get('n_priority_regulators', 0)
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
    
    # FINAL INTEGRATED SCORE — informative only since v1.2.0
    # This weighted sum no longer drives the ranking (CASCADE_CRITERIA does). It is still
    # written to evaluation_results.tsv so runs can be compared against earlier results.
    # Note it is normalized within the evaluated sample, so its values are NOT comparable
    # between runs — only the ordering within one run is meaningful.

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

    # ── Rank by the lexicographic cascade ───────────────────────────────────────
    # Criteria in CASCADE_CRITERIA are compared in order; the first that differs decides.
    # Signs are folded into helper columns so a plain ascending sort implements the
    # cascade, and rounding keeps float noise from splitting genuinely tied rows.
    rank_cols = []
    for name, sign in _ACTIVE_CASCADE:
        col = f'_rank_{name}'
        results_df[col] = (sign * results_df[name]).round(_CASCADE_ROUND)
        rank_cols.append(col)

    results_df = results_df.sort_values(rank_cols, ascending=True)

    # Winners: every row sharing the best key tuple. Compared on the rounded helper
    # columns, so combinations differing only by float noise still count as tied.
    best_key_values   = results_df.iloc[0][rank_cols]
    is_best           = (results_df[rank_cols] == best_key_values).all(axis=1)
    best_combinations = results_df[is_best].copy()

    # Helper columns are internal — drop them before writing anything out
    results_df = results_df.drop(columns=rank_cols)

    # Reorder columns to put final_score first
    cols = list(results_df.columns)
    cols.remove('final_score')
    cols.insert(1, 'final_score')  # After combination_id
    results_df = results_df[cols]

    # Save results — truncate to top_results rows to keep file manageable
    if verbose:
        print(f"\n9. Saving results to: {output_file}")
    if top_results is not None and len(results_df) > top_results:
        if verbose:
            print(f"   Saving top {top_results:,} of {len(results_df):,} combinations evaluated")
        results_df.head(top_results).to_csv(output_file, sep='\t', index=False)
    else:
        if verbose:
            print(f"   Saving all {len(results_df):,} combinations evaluated")
        results_df.to_csv(output_file, sep='\t', index=False)
    
    # Generate rules_by_gene_evaluated.tsv with winning rules
    if verbose:
        print(f"\n10. Generating rules_by_gene_evaluated.tsv with winning rules...")
    
    # best_combinations was determined by the cascade above

    # Create a new dataframe with the winning rules
    winning_rules_data = []

    for gene in genes:
        # Count how many top combinations used each (pos, rule) for this gene
        rule_counts: Dict[tuple, int] = {}
        for _, row in best_combinations.iterrows():
            key = (row[f'{gene}_position'], row[f'{gene}_rule'])
            rule_counts[key] = rule_counts.get(key, 0) + 1

        # Look up original GEP scores for each winning rule
        for (pos, rule), count in rule_counts.items():
            original_row = df[(df['Gene'] == gene) &
                              (df['Position'] == pos) &
                              (df['Rule'] == rule)]

            if not original_row.empty:
                orig = original_row.iloc[0]

                entry = {
                    'Gene': gene,
                    'Position': int(pos),
                    'Rule': rule,
                }

                # Add GEP columns in the original order
                for col in df.columns:
                    if col not in ['Gene', 'Position', 'Rule']:
                        entry[col] = orig[col]

                # How many of the tied best combinations use this rule for this gene.
                # Per gene these counts sum to n_tied_combos, so a single row carrying
                # the full count means the cascade settled that gene outright, while
                # several rows mean those alternatives are indistinguishable under
                # every criterion — genuinely ambiguous, not arbitrarily chosen.
                entry['n_top_combos'] = count
                entry['n_tied_combos'] = len(best_combinations)
                entry['consensus'] = round(count / len(best_combinations), 4)

                winning_rules_data.append(entry)

    # Create dataframe and sort by gene, then descending consensus (most used rule first)
    winning_rules_df = pd.DataFrame(winning_rules_data)
    winning_rules_df = winning_rules_df.sort_values(
        ['Gene', 'n_top_combos', 'Position'], ascending=[True, False, True]
    )

    # Reorder columns: GEP columns first (original order), then the consensus block
    original_columns = ['Gene', 'Position', 'Rule']

    for col in df.columns:
        if col not in original_columns and col in winning_rules_df.columns:
            original_columns.append(col)

    original_columns.extend(['n_top_combos', 'n_tied_combos', 'consensus'])

    winning_rules_df = winning_rules_df[original_columns]
    
    # Save to file in same directory as output_file
    output_dir = os.path.dirname(output_file)
    winning_rules_file = os.path.join(output_dir, "rules_by_gene_evaluated.tsv")
    winning_rules_df.to_csv(winning_rules_file, sep='\t', index=False)

    # Search trace — how the answer was reached, not just what it was
    if write_log:
        log_file = os.path.join(output_dir, "evaluation_log.tsv")
        trace = (hill_info['trace'] if hill_info is not None else sample_trace)
        log_header = {
            'generated':   time.strftime('%Y-%m-%d %H:%M:%S'),
            'rules_file':  rules_file,
            'matrix_file': binarized_matrix_file or 'none',
            'genes':       len(genes),
            'strategy':    ('exhaustive' if is_exhaustive
                            else ('hill' if hill_info is not None else 'sample')),
            'search_space_before_prune': f"{prune_info.get('before', total_combinations):.6e}",
            'search_space_after_prune':  f"{total_combinations:.6e}",
            'prune_applied':   prune_info['applied'],
            'prune_evaluations': prune_info.get('evaluations', 0),
            'combinations_evaluated': len(results),
            'processes':   n_processes,
            'cascade':     ' > '.join(f"{n}({'min' if sg > 0 else 'max'})"
                                      for n, sg in _ACTIVE_CASCADE),
            'priority_regulators': (','.join(sorted(priority_regulators))
                                    if priority_regulators else 'none'),
        }
        write_evaluation_log(log_file, log_header, trace)
        if verbose and trace:
            print(f"   Search trace ({len(trace)} rows) saved to: {log_file}")
    
    if verbose:
        print(f"   Saved {len(winning_rules_df)} winning rules to: {winning_rules_file}")
        n_best_combos = len(best_combinations)
        if n_best_combos > 1:
            per_gene = winning_rules_df.groupby('Gene').size()
            settled = int((per_gene == 1).sum())
            ambiguous = per_gene[per_gene > 1]

            print(f"   {n_best_combos} combinations tied on ALL cascade criteria, so "
                  f"the file lists every rule they use.")
            print(f"   Per gene, n_top_combos sums to {n_best_combos} "
                  f"(consensus = n_top_combos / {n_best_combos}):")
            print(f"   - {settled} gene(s) settled outright: one rule used by all "
                  f"{n_best_combos} tied combinations")
            if len(ambiguous):
                print(f"   - {len(ambiguous)} gene(s) still ambiguous, no criterion "
                      f"separates their alternatives:")
                for gene, n_rules in ambiguous.items():
                    top = winning_rules_df[winning_rules_df['Gene'] == gene].iloc[0]
                    print(f"       {gene}: {n_rules} rules, best one used by "
                          f"{top['n_top_combos']}/{n_best_combos} "
                          f"({top['consensus']:.0%})")
    
    # Print summary
    if verbose:
        print("\n" + "="*80)
        print("EVALUATION SUMMARY")
        print("="*80)
        print(f"Best combination (lexicographic tie-break cascade):")
        best = results_df.iloc[0]
        print(f"  Criteria in priority order — first one that differs decides:")
        for rank, (name, sign) in enumerate(_ACTIVE_CASCADE, 1):
            direction = "lower" if sign > 0 else "higher"
            print(f"    {rank}. {name:24s} {best[name]:>9.4f}  ({direction} is better)")
        print(f"\n  Attractor metrics:")
        print(f"    Number of basins: {int(best['n_basins'])}")
        print(f"    Avg cycle length: {best['avg_cycle_length']:.2f}")
        print(f"    Basin entropy: {best['basin_entropy']:.4f}")
        print(f"    Fixed point ratio: {best['fixed_point_ratio']:.2f}")
        if binarized_matrix is not None:
            print(f"    Final state concordance: {best['final_state_concordance']:.4f}")
            print(f"    Basin share of matching attractor: "
                  f"{best['concordance_basin_frac']:.4f}")
        print(f"\n  Parsimony metrics:")
        print(f"    Total literals: {int(best['total_literals'])} (fewer = simpler)")
        print(f"    Total NOT operators: {int(best['total_nots'])} (fewer = less repression)")
        print(f"    Avg K (connectivity): {best['avg_k']:.2f} (optimal ≈ 2.0)")
        print(f"    K distance from 2.0: {best['k_distance_from_2']:.3f} (closer to 0 = better)")
        
        # Legacy scores — informative only, they no longer decide the ranking
        print(f"\n  Reference (not used for ranking):")
        print(f"    final_score (weighted sum): {best['final_score']:.4f}")
        if 'composite_score_old' in best:
            print(f"    old composite_score: {best['composite_score_old']:.4f}")
        if results_df['final_score'].idxmin() != results_df.index[0]:
            print(f"    Note: the old weighted sum would have picked a "
                  f"different combination")

        print(f"\nAll {len(results_df)} combinations saved to output file")
        print(f"Ranked by lexicographic cascade over {len(_ACTIVE_CASCADE)} criteria")

        elapsed = time.perf_counter() - t_start
        n_evaluated = len(results)
        combos_per_sec = n_evaluated / elapsed if elapsed > 0 else 0
        mins, secs = divmod(elapsed, 60)
        print(f"\n{'─'*40}")
        print(f"  Execution time : {int(mins):02d}m {secs:05.2f}s  ({elapsed:.1f}s total)")
        batches_str = f" in {batch_num} batches" if batch_num > 1 else ""
        print(f"  Combinations   : {n_evaluated:,} evaluated{batches_str}")
        print(f"  Throughput     : {combos_per_sec:,.0f} combinations/sec")
        if combos_per_sec > 0:
            print(f"  Per combination: {1000/combos_per_sec:.2f} ms")
        if is_exhaustive:
            print(f"  Coverage       : 100% (exhaustive search)")
        elif hill_info is not None:
            n_restarts = hill_info['restarts']
            at_best    = hill_info['restarts_at_best']
            rate       = at_best / max(1, n_restarts)
            print(f"  ── Hill-climbing coverage ────────")
            print(f"  Restarts       : {n_restarts}")
            print(f"  Reached best   : {at_best} ({rate:.1%} of restarts)")
            print(f"  Combos at best : {k_converged}")
            # The 1-e^-k formula needs uniform sampling and does not apply here. But
            # restarts ARE independent uniform draws of a STARTING point, so they bound
            # something better: basin size. An optimum reachable from a fraction x of
            # starting points survives R restarts unseen with probability (1-x)^R.
            # Unlike 1-e^-k this is not circular — it describes the landscape rather
            # than the best value already observed.
            target = hill_info.get('target_restarts')
            if target:
                if n_restarts >= target:
                    print(f"  Stopped on   : restart target reached "
                          f"({target:,} needed)")
                else:
                    print(f"  Stopped on   : evaluation budget — only {n_restarts:,} "
                          f"of the {target:,} restarts needed")
                    print(f"                 for the requested coverage. Raise "
                          f"--max-combos to reach it.")
            if n_restarts > 0:
                print(f"  Any optimum missed must have a basin smaller than:")
                for label, conf in [("95%  ", 0.95), ("99%  ", 0.99), ("99.9%", 0.999)]:
                    x = 1 - (1 - conf) ** (1 / n_restarts)
                    print(f"    at {label} conf: {x:>8.3%} of starting points")
                if at_best:
                    print(f"  Missing one as reachable as this best: "
                          f"p = {(1 - rate) ** n_restarts:.1e}")
                if rate < 0.10:
                    print(f"  Note: a {rate:.1%} convergence rate means a rugged "
                          f"landscape —")
                    print(f"        most restarts settle on worse local optima.")
        else:
            print(f"  ── Statistical coverage ──────────")
            print(f"  k (top rank)   : {k_converged} combos sharing the best cascade key")
            print(f"  Confidence     : {confidence_achieved:.2%}")
            if k_converged > 0:
                for label, p in [("99%  ", 0.99), ("99.5%", 0.995), ("99.99%", 0.9999)]:
                    n_need = math.ceil(-math.log(1.0 - p) * n_evaluated / k_converged)
                    extra  = max(0, n_need - n_evaluated)
                    status = "(done)" if extra == 0 else f"({extra:,} more needed)"
                    print(f"  For {label}    : {n_need:,} total {status}")
        print(f"{'─'*40}")
        print("="*80)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Evaluate Boolean rule combinations using attractor metrics.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: Evaluate ALL rules with max score, ALL combinations
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv -m binarized_matrix.tsv -v

  # Limit to top 5 rules per gene (instead of all tied for best)
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --top-n 5 -v

  # Force an unbiased uniform sample (the only mode with a confidence estimate)
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --search sample -v

  # Limit to 500 combinations max (instead of all)
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --max-combos 500 -v

  # Both limits combined
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --top-n 3 --max-combos 200 -v

  # Custom output directory
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv -o /path/to/results.tsv -v

  # Use custom score column for ranking
  python3 3.BNI3_Evaluate_rules.py -i rules_by_gene.tsv --score-col MSE -v

Ranking — lexicographic tie-break cascade:
  Combinations are compared criterion by criterion, IN ORDER. The first criterion
  that differs decides the winner; the remaining ones are never consulted. Metrics
  are compared raw (no normalization), so rankings are reproducible across runs.

    1. final_state_concordance  (higher) — agreement with the observed final state
    2. n_basins                 (lower)  — fewer basins
    3. avg_cycle_length         (lower)  — shorter attractors
    4. fixed_point_ratio        (higher) — more fixed points
    5. basin_entropy            (lower)  — one dominant basin
    6. k_distance_from_2        (lower)  — Kauffman K ≈ 2
    7. total_literals           (lower)  — simpler rules
    8. total_nots               (lower)  — less repression

  With --priority-regulators, a 9th criterion (n_priority_regulators, higher is
  better) is inserted at position 6 — after data fit and attractor structure,
  before parsimony. It prefers rules built on the listed regulators without ever
  requiring them.

  The order is defined once in CASCADE_CRITERIA near the top of this file; edit it
  there to re-prioritize. Data agreement comes first because it is the only
  criterion anchored in observation — the rest encode priors about network shape.

  Individual metrics (for reference):
  - n_basins: Number of attractor basins
  - avg_cycle_length: Average attractor cycle length
  - basin_entropy: Basin size distribution uniformity
  - final_state_concordance: Best match with last timepoint across ALL attractors
      (for a cyclic attractor, its best-matching state is used)
  - concordance_basin_frac: Basin share of the attractor that produced that match
      (informative: a match in a tiny basin is weaker evidence than in the dominant one)
  - fixed_point_ratio: Proportion of fixed point attractors
  - total_literals: Total gene mentions (parsimony)
  - total_nots: Total NOT operators (repression level)
  - avg_k: Average regulators per gene (connectivity)
  - k_distance_from_2: Distance from optimal K=2
  - final_score: Legacy weighted sum of normalized metrics. Informative ONLY — it no
      longer ranks anything. Normalized within the evaluated sample, so its values are
      not comparable between runs.
  - composite_score_old: Original composite score (reference)

Output format:
  Two TSV files are generated:

  1. evaluation_results.tsv:
     - Sorted by the lexicographic cascade described above
     - First row is the overall best combination
     - Truncated to --top-results rows (default 10,000) to keep file manageable

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
                              'Default: None (AUTOMATIC - calculates intelligent default, up to 1,000,000)')
    optional.add_argument('--max-iter', type=int, default=1000,
                         help='Maximum iterations for attractor search (default: 1000)')
    optional.add_argument('-n', '--processes', type=int, default=None,
                         help='Number of parallel worker processes (default: None = all available cores)')
    optional.add_argument('--top-results', type=int, default=10000,
                         help='Number of top-ranked combinations to write to evaluation_results.tsv '
                              '(default: 10000). Use 0 to save all.')
    optional.add_argument('--score-col', type=str, default='Score',
                         help='Column name for rule ranking (default: Score)')
    optional.add_argument('--confidence', type=float, default=0.999,
                         help='Target statistical confidence for finding the global optimum '
                              '(default: 0.99). Adaptive batching adds samples until this '
                              'confidence is reached. Set to 0 to disable. Range: 0.0–1.0.')
    # Search strategy and pruning
    strategy = parser.add_argument_group('Search strategy')
    strategy.add_argument('--search', choices=['auto', 'sample', 'hill'],
                         default='auto',
                         help="How to explore the combination space. 'auto' (default) "
                              "enumerates everything when it fits and otherwise runs "
                              "hill climbing. 'sample' forces uniform random sampling, "
                              "which is unbiased and is the only mode the 1-e^-k "
                              "confidence estimate applies to. 'hill' forces "
                              "random-restart hill climbing")
    strategy.add_argument('--restarts', type=int, default=None,
                         help='Number of random restarts for --search hill '
                              '(default: as many as the evaluation budget allows)')
    strategy.add_argument('--no-prune-concordance', dest='prune_concordance',
                         action='store_false',
                         help='Disable the exact prune that keeps only rules able to '
                              'reproduce the observed final state. The prune is safe '
                              '(pruned combinations can never win) and costs sum(n) '
                              'rule evaluations; disable it only to reproduce older runs')
    strategy.add_argument('--target-basin', type=float, default=0.01,
                         help='Stopping criterion for --search hill: keep restarting '
                              'until any optimum not found must sit in a basin smaller '
                              'than this fraction of starting points, at --confidence. '
                              'Default 0.01 (1%%). Restarts, not evaluations, are what '
                              'the guarantee depends on; --max-combos is only a safety '
                              'cap. Set 0 to fall back to the evaluation budget')
    strategy.add_argument('--seed', type=int, default=42,
                         help='Random seed for the hill-climbing starting points '
                              '(default: 42). Runs are reproducible; change it to '
                              'check how stable the answer is across independent searches')
    strategy.add_argument('--no-log', dest='write_log', action='store_false',
                         help='Do not write evaluation_log.tsv, the per-restart (or '
                              'per-batch) trace of how the search progressed')
    strategy.add_argument('--priority-regulators', type=str, default=None,
                         help='Comma-separated gene names to prefer as regulators, '
                              'e.g. transcription factors: HB6,MYB44,ABF3. Adds a '
                              'cascade criterion at position 6 — after data fit and '
                              'attractor structure, before parsimony — so it is a '
                              'preference, never a requirement. A file with one name '
                              'per line is also accepted')

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
            top_results=args.top_results if args.top_results != 0 else None,
            score_column=args.score_col,
            confidence_target=args.confidence,
            prune_concordance=args.prune_concordance,
            search=args.search,
            restarts=args.restarts,
            target_basin=args.target_basin,
            seed=args.seed,
            priority_regulators=load_priority_regulators(args.priority_regulators),
            write_log=args.write_log,
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