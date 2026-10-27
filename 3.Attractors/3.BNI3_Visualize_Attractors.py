#!/usr/bin/env python3
"""
Boolean Network Attractor Visualizer
Creates visualizations of Boolean network attractors from attractor analysis results.
Part of the Boolean Network Inference (BNI) pipeline.
"""

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import seaborn as sns
import networkx as nx
from pathlib import Path


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)

def extract_mutation_suffix(filename):
    """
    Extract mutation suffix from filename if present.
    
    Args:
        filename (str): Input filename
        
    Returns:
        str: Mutation suffix (e.g., "_ABF3_1") or empty string
    """
    # Get base name without extension
    base_name = Path(filename).stem
    
    # Check if filename follows pattern: attractors_GENE_VALUE_GENE_VALUE...
    if base_name.startswith('attractors_'):
        # Extract everything after 'attractors'
        suffix = base_name[len('attractors'):]
        if suffix:  # If there's content after 'attractors'
            return suffix
    
    return ""

# Same expression the inference and evaluation scripts use: any character that is not
# a letter, digit or underscore becomes "_", and a leading digit gets an "_" prefix.
# The attractors file carries sanitized names, a binarized matrix does not, so without
# this "SnRK2.8" never matches the network gene "SnRK2_8" and reads as 0.
_GENE_NAME_RE = r'\W|^(?=\d)'


def sanitize_gene_name(name):
    """Turn one gene name into a valid Python identifier (idempotent)."""
    return re.sub(_GENE_NAME_RE, '_', str(name))


def sanitize_gene_names(names, source):
    """Sanitize gene names, refusing silently ambiguous results."""
    clean = [sanitize_gene_name(n) for n in names]

    groups = {}
    for original, cleaned in zip(names, clean):
        groups.setdefault(cleaned, set()).add(str(original))
    collisions = {k: v for k, v in groups.items() if len(v) > 1}

    if collisions:
        print(f"\nERROR: gene names in {source} become ambiguous once sanitized to "
              f"Python identifiers.", file=sys.stderr)
        for cleaned, originals in list(collisions.items())[:10]:
            print(f"  {sorted(originals)} all become '{cleaned}'", file=sys.stderr)
        print("  Rename them upstream so they stay distinct.", file=sys.stderr)
        sys.exit(1)

    return clean


def read_final_matrix_state(matrix_file, gene_cols, verbose):
    """
    Read binarized matrix and return the last row as a boolean state vector,
    aligned to the gene order defined by gene_cols.

    Args:
        matrix_file (str): Path to binarized matrix TSV
        gene_cols (list): Ordered gene names from attractors file
        verbose (bool): Enable verbose output

    Returns:
        tuple: (boolean state vector, sample label string)
    """
    df = pd.read_csv(matrix_file, sep='\t')
    first_col = df.columns[0]
    if not pd.api.types.is_numeric_dtype(df[first_col]):
        sample_name = str(df[first_col].iloc[-1])
        df = df.drop(columns=[first_col])
    else:
        sample_name = f"Sample_{len(df)}"

    # Match the sanitization used everywhere else; otherwise a column named
    # "SnRK2.8" never matches the network gene "SnRK2_8" and silently reads as 0.
    df.columns = sanitize_gene_names(df.columns, os.path.basename(matrix_file))

    matrix_genes = list(df.columns)
    missing = sorted(set(gene_cols) - set(matrix_genes))
    if missing:
        # Defaulting to 0 corrupts the state vector, so this is never silent.
        print(f"\nWARNING: {len(missing)} network gene(s) absent from "
              f"{os.path.basename(matrix_file)} and defaulted to 0: "
              f"{', '.join(missing)}", file=sys.stderr)

    last_row = df.iloc[-1]
    state = [bool(int(last_row[g])) if g in matrix_genes else False for g in gene_cols]
    log_message(f"Final matrix state: {sample_name}", verbose)
    return state, sample_name


def build_attractors_dict(df, gene_cols):
    """Build a simple attractors lookup dict from the attractors dataframe."""
    attractors_dict = {}
    for att_id in df['attractor_id'].unique():
        att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
        states = [[bool(row[g]) for g in gene_cols] for _, row in att_data.iterrows()]
        attractors_dict[att_id] = {'states': states}
    return attractors_dict


def find_target_attractor_id(final_state, attractors_dict, gene_cols,
                              gene_rules=None, max_steps=1000, verbose=False):
    """
    Determine which attractor the final matrix state belongs to or converges to.

    Args:
        final_state (list): Boolean state vector of the last matrix row
        attractors_dict (dict): Attractor states keyed by attractor_id
        gene_cols (list): Ordered gene names
        gene_rules (dict|None): Boolean rules per gene (required for transient states)
        max_steps (int): Maximum simulation steps
        verbose (bool): Enable verbose output

    Returns:
        int|None: attractor_id if found, else None
    """
    def to_str(s):
        return ''.join('1' if x else '0' for x in s)

    target_str = to_str(final_state)

    # Direct membership check
    for att_id, info in attractors_dict.items():
        if any(to_str(s) == target_str for s in info['states']):
            log_message(f"Final matrix state IS attractor {att_id}", verbose)
            return att_id

    if gene_rules is None:
        log_message("Final state not in any attractor; no rules provided — cannot simulate.", verbose)
        return None

    # Simulate forward
    current = final_state[:]
    visited = set()
    for step in range(max_steps):
        s = to_str(current)
        for att_id, info in attractors_dict.items():
            if any(to_str(a) == s for a in info['states']):
                log_message(f"Final matrix state → attractor {att_id} in {step} steps", verbose)
                return att_id
        if s in visited:
            log_message("Final matrix state enters an unknown cycle.", verbose)
            return None
        visited.add(s)

        next_state = {}
        for gene in gene_cols:
            if gene in gene_rules:
                state_dict = {gene_cols[i]: current[i] for i in range(len(gene_cols))}
                rule = (gene_rules[gene]
                        .replace('&', ' and ')
                        .replace('|', ' or ')
                        .replace('~', ' not '))
                for g, v in state_dict.items():
                    rule = rule.replace(g, str(v))
                try:
                    next_state[gene] = bool(eval(rule))
                except Exception:
                    next_state[gene] = False
            else:
                next_state[gene] = current[gene_cols.index(gene)]
        current = [next_state[g] for g in gene_cols]

    log_message(f"Final matrix state did not converge within {max_steps} steps.", verbose)
    return None


def read_rules_simple(rules_file, verbose):
    """
    Minimal rules reader that returns a gene→rule dict.
    Supports both selected_rules.tsv (Gene/Rule columns) and full rules table (with Position).
    """
    df = pd.read_csv(rules_file, sep='\t', encoding='utf-8')
    if 'Position' in df.columns:
        selected = {}
        for gene in df['Gene'].unique():
            sub = df[df['Gene'] == gene]
            top = sub[sub['Position'] == 1]
            row = top.iloc[0] if not top.empty else sub.loc[sub['Position'].idxmin()]
            selected[str(row['Gene'])] = str(row['Rule'])
        log_message(f"Rules loaded from full table ({len(selected)} genes)", verbose)
        return selected
    else:
        rules = {str(r['Gene']): str(r['Rule']) for _, r in df.iterrows()}
        log_message(f"Rules loaded from selected rules file ({len(rules)} genes)", verbose)
        return rules


# ---------------------------------------------------------------------------
# Predecessor states around the attractors
#
# BoolNet-style figures draw the whole state transition graph. That is 2^N nodes
# and is unreadable past ~10 genes, so instead a few representative predecessors
# are drawn per attractor node: states that map INTO it in one update.
#
# The fan-in is the problem, not the cost. On a 15-gene network the Control fixed
# point has 2,559 direct predecessors and 24,576 within three steps — essentially
# the whole space. So the predecessors are SELECTED, not enumerated:
#
#   1. states observed in the binarized matrix are always kept (they are the only
#      predecessors that are actual measurements rather than possibilities);
#   2. the remaining slots go to the states closest to their successor in Hamming
#      distance, which read as "flip these few genes and the system still returns".
#
# Ties break on the state code, so the figure is deterministic.
# ---------------------------------------------------------------------------

_MAX_PREDECESSOR_GENES = 24     # 2^24 states; beyond this the table stops being cheap


def _popcount(values):
    """Vectorized bit count, for Hamming distances between packed states."""
    v = np.asarray(values).astype(np.uint32)
    v = v - ((v >> 1) & np.uint32(0x55555555))
    v = (v & np.uint32(0x33333333)) + ((v >> 2) & np.uint32(0x33333333))
    v = (v + (v >> 4)) & np.uint32(0x0F0F0F0F)
    return ((v * np.uint32(0x01010101)) >> 24).astype(np.int64)


def state_to_code(state, n_genes):
    """Pack a boolean state vector into an integer, first gene in the high bit."""
    code = 0
    for i, value in enumerate(state):
        if value:
            code |= 1 << (n_genes - 1 - i)
    return code


def code_to_state(code, n_genes):
    """Unpack an integer state code back into a boolean vector."""
    return [bool((code >> (n_genes - 1 - i)) & 1) for i in range(n_genes)]


def _rule_to_numpy_expr(rule):
    """Rewrite a rule so it evaluates element-wise over numpy boolean arrays."""
    expr = str(rule)
    expr = re.sub(r'\bnot\b', ' ~ ', expr)
    expr = re.sub(r'\band\b', ' & ', expr)
    expr = re.sub(r'\bor\b', ' | ', expr)
    return expr


def build_transition_table(genes, gene_rules, verbose=False):
    """
    Evaluate the synchronous update of every one of the 2^N states.

    Each rule is evaluated once against numpy boolean arrays spanning the whole
    state space rather than once per state, so the cost is N vectorized ops
    instead of N * 2^N interpreted ones.

    Returns:
        np.ndarray: nxt[s] is the state code reached from state code s.
    """
    n_genes = len(genes)
    if n_genes > _MAX_PREDECESSOR_GENES:
        raise ValueError(
            f"Predecessor states need the full 2^{n_genes} transition table "
            f"({2 ** n_genes:,} states), which is beyond the {_MAX_PREDECESSOR_GENES}-gene "
            f"limit. Drop --predecessors for this network."
        )

    missing = [g for g in genes if g not in gene_rules]
    if missing:
        raise ValueError(f"No rule for {len(missing)} gene(s): {', '.join(missing)}")

    total = 1 << n_genes
    log_message(f"Building transition table over {total:,} states...", verbose)

    index = np.arange(total, dtype=np.int64)
    env = {g: ((index >> (n_genes - 1 - i)) & 1).astype(bool)
           for i, g in enumerate(genes)}

    nxt = np.zeros(total, dtype=np.int64)
    for i, gene in enumerate(genes):
        rule = str(gene_rules[gene]).strip()
        if rule in ('True', 'true', '1'):
            values = np.ones(total, dtype=bool)
        elif rule in ('False', 'false', '0'):
            values = np.zeros(total, dtype=bool)
        else:
            values = eval(_rule_to_numpy_expr(rule), {'__builtins__': {}}, env)
            values = np.broadcast_to(np.asarray(values, dtype=bool), (total,))
        nxt |= values.astype(np.int64) << (n_genes - 1 - i)

    return nxt


def clamped_state_mask(genes, gene_rules, total_states):
    """
    States consistent with every gene that a constant rule pins to a value.

    A mutation is encoded as a constant rule ("MYB44 -> 0"), and
    1.BNI3_Attractors.py enumerates only the states where the clamped gene already
    holds that value: 2^12 rather than 2^13 for one knockout. Counting the other
    half here would make the basin sizes in this figure disagree with the ones in
    the attractors file, which is exactly the kind of mismatch the counts are
    meant to resolve.

    Returns:
        np.ndarray of bool, or None when no gene is clamped.
    """
    n_genes = len(genes)
    clamped = {}
    for i, gene in enumerate(genes):
        rule = str(gene_rules.get(gene, '')).strip()
        if rule in ('1', 'True', 'true', 'TRUE'):
            clamped[i] = 1
        elif rule in ('0', 'False', 'false', 'FALSE'):
            clamped[i] = 0

    if not clamped:
        return None

    index = np.arange(total_states, dtype=np.int64)
    mask = np.ones(total_states, dtype=bool)
    for i, value in clamped.items():
        bit = (index >> (n_genes - 1 - i)) & 1
        mask &= (bit == value)
    return mask


def entry_point_counts(nxt, drawn_codes, valid_mask=None):
    """
    Split every state the figure does not draw by where it enters the figure.

    Dynamics are deterministic, so each state has exactly one forward path, and
    that path must eventually reach the attractor — which is drawn. Every state
    therefore has exactly one *first* drawn state it lands on, and grouping by it
    partitions the undrawn remainder with no overlap:

        sum(counts.values()) + len(drawn_codes) == basin size

    This is what per-node "+N" placeholders could not do. Those counted direct
    predecessors, which are nested inside each other's upstream sets, so adding
    them up overshot the basin.

    Returns:
        {drawn_code: how many undrawn states enter the figure there}
    """
    total = nxt.size
    drawn = np.asarray(sorted(set(int(c) for c in drawn_codes)), dtype=np.int64)

    in_drawn = np.zeros(total, dtype=bool)
    in_drawn[drawn] = True
    owner = np.full(total, -1, dtype=np.int64)
    owner[drawn] = drawn               # a drawn state is its own entry point

    frontier = drawn
    while frontier.size:
        candidates = np.flatnonzero(np.isin(nxt, frontier))
        candidates = candidates[(owner[candidates] < 0) & ~in_drawn[candidates]]
        if valid_mask is not None:
            candidates = candidates[valid_mask[candidates]]
        if candidates.size == 0:
            break
        # The successor is already resolved, so its entry point is this one's too.
        owner[candidates] = owner[nxt[candidates]]
        frontier = candidates

    outside = np.flatnonzero((owner >= 0) & ~in_drawn)
    return {int(code): int(np.count_nonzero(owner[outside] == code))
            for code in drawn}


def basin_layer_sizes(nxt, seed_codes, total_states, valid_mask=None):
    """
    Size of every backward layer of one attractor's basin, walked to exhaustion.

    This is the unconditioned count, and it is NOT what build_predecessor_layers
    reports. That function expands only the handful of states the figure kept, so
    its "available at step 2" means "predecessors of the three states I drew", not
    "states two steps from the attractor". Summing those reconciles with nothing.

    Here every state of a layer is expanded, so the layers partition the basin:
    len(seed_codes) + sum(layers) == basin size.

    Returns:
        (layers, reachable) with layers[d-1] the number of states exactly d steps
        upstream, and reachable the total including the attractor itself.
    """
    seen = np.zeros(total_states, dtype=bool)
    seeds = np.asarray(sorted(set(seed_codes)), dtype=np.int64)
    seen[seeds] = True

    frontier = seeds
    layers = []
    while frontier.size:
        candidates = np.flatnonzero(np.isin(nxt, frontier))
        candidates = candidates[~seen[candidates]]
        if valid_mask is not None:
            candidates = candidates[valid_mask[candidates]]
        if candidates.size == 0:
            break
        layers.append(int(candidates.size))
        seen[candidates] = True
        frontier = candidates

    return layers, int(seen.sum())


def build_predecessor_layers(nxt, seed_codes, depth, breadth,
                             pinned=None, verbose=False, valid_mask=None):
    """
    Walk backwards from the attractor, keeping a few predecessors per state.

    Args:
        nxt: transition table from build_transition_table
        seed_codes: state codes of one attractor (one for a fixed point, k for a cycle)
        depth: how many update steps to walk back
        breadth: how many predecessors to keep per state per step
        pinned: state codes that are kept regardless of Hamming distance
        valid_mask: bool array of states a constant (mutation) rule allows. A
            knockout figure must not show states where the knocked-out gene is ON:
            step 1 never counted them, so drawing them puts states in the figure
            that are outside the state space it reports.

    Returns:
        (layers, stats). layers[d] is {'chosen': [(code, successor_code)],
        'omitted': {successor_code: how many of its predecessors were left out}}
        for states d steps upstream; stats carries the totals before selection.
        The omitted counts are what lets the figure state its own scale instead of
        implying that the few states drawn are all there are.
    """
    pinned = set(pinned or ())
    claimed = set(seed_codes)          # never draw a state twice
    frontier = list(seed_codes)
    layers = []
    stats = {'available': [], 'drawn': [], 'total_available': 0, 'total_drawn': 0}

    for _ in range(depth):
        if not frontier:
            break

        targets = np.asarray(frontier, dtype=np.int64)
        # One pass over the table finds every predecessor of the whole frontier.
        candidates = np.flatnonzero(np.isin(nxt, targets))
        candidates = candidates[~np.isin(candidates, np.asarray(sorted(claimed),
                                                                dtype=np.int64))]
        if valid_mask is not None:
            candidates = candidates[valid_mask[candidates]]
        stats['available'].append(int(candidates.size))
        if candidates.size == 0:
            break

        successors = nxt[candidates]

        chosen = []
        omitted = {}
        for target in frontier:
            mine = candidates[successors == target]
            if mine.size == 0:
                continue
            mine_pinned = [int(c) for c in mine if int(c) in pinned]
            rest = [int(c) for c in mine if int(c) not in pinned]
            rest.sort(key=lambda c: (int(_popcount(np.array([c ^ target]))[0]), c))
            keep = mine_pinned + rest[:max(0, breadth - len(mine_pinned))]
            kept_here = 0
            for code in keep:
                if code not in claimed:
                    claimed.add(code)
                    chosen.append((code, int(target)))
                    kept_here += 1
            left = int(mine.size) - kept_here
            if left > 0:
                omitted[int(target)] = left

        stats['drawn'].append(len(chosen))
        if not chosen:
            # Nothing kept, but the count of what was passed over still matters.
            if omitted:
                layers.append({'chosen': [], 'omitted': omitted})
                stats['total_available'] += sum(omitted.values())
            break
        layers.append({'chosen': chosen, 'omitted': omitted})
        frontier = [code for code, _ in chosen]

    stats['total_available'] = sum(stats['available'])
    stats['total_drawn'] = sum(stats['drawn'])
    return layers, stats


def generate_attractor_colors(n_attractors):
    """
    Generate distinctive colors for each attractor
    
    Args:
        n_attractors (int): Number of attractors
        
    Returns:
        list: List of color codes
    """
    if n_attractors <= 10:
        # Use qualitative color palette for small numbers
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', 
                 '#ffff33', '#a65628', '#f781bf', '#999999', '#1f78b4']
        return colors[:n_attractors]
    else:
        # Generate colors using colormap for larger numbers
        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab20')
        return [cmap(i / n_attractors) for i in range(n_attractors)]


def read_attractors_file(attractors_file, verbose):
    """
    Read the attractors TSV file generated by 1.BNI3_Attractors.py
    
    Args:
        attractors_file (str): Path to the attractors TSV file
        verbose (bool): Enable verbose output
        
    Returns:
        tuple: (DataFrame, list of gene names)
    """
    try:
        log_message(f"Reading attractors file: {attractors_file}", verbose)
        
        df = pd.read_csv(attractors_file, sep='\t')
        log_message(f"Loaded {len(df)} attractor states", verbose)
        
        # Validate required columns
        required_cols = ['attractor_id', 'type', 'cycle_length', 'step_in_cycle', 'binary_state']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Extract gene names (columns that are not metadata)
        metadata_cols = ['attractor_id', 'type', 'cycle_length', 'step_in_cycle', 'binary_state', 'basin_size', 'basin_percentage']
        gene_cols = [col for col in df.columns if col not in metadata_cols]
        
        log_message(f"Found {len(gene_cols)} genes: {', '.join(gene_cols)}", verbose)
        
        return df, gene_cols
        
    except Exception as e:
        raise ValueError(f"Error reading attractors file: {str(e)}")


def get_basin_sizes_from_df(df, verbose):
    """
    Extract basin sizes from the dataframe (now included in the TSV)
    
    Args:
        df (pd.DataFrame): Attractors dataframe
        verbose (bool): Enable verbose output
        
    Returns:
        dict: Mapping of attractor_id to basin size
    """
    basin_sizes = {}
    
    # Check if basin_size column exists
    if 'basin_size' in df.columns:
        for att_id in df['attractor_id'].unique():
            basin_size = df[df['attractor_id'] == att_id]['basin_size'].iloc[0]
            basin_sizes[att_id] = int(basin_size)
        
        log_message(f"Basin sizes loaded from file: {basin_sizes}", verbose)
    else:
        # Fallback to equal distribution if columns don't exist
        log_message("Warning: basin_size column not found. Using equal distribution fallback.", verbose)
        unique_attractors = df['attractor_id'].unique()
        gene_cols = [col for col in df.columns if col not in 
                    ['attractor_id', 'type', 'cycle_length', 'step_in_cycle', 'binary_state', 'basin_size', 'basin_percentage']]
        total_states = 2**len(gene_cols)
        basin_size_per_attractor = total_states // len(unique_attractors)
        basin_sizes = {aid: basin_size_per_attractor for aid in unique_attractors}
    
    return basin_sizes


def create_attractor_heatmap(df, gene_cols, basin_sizes, output_path,
                             verbose=False, svg_output=False):
    """
    Create a heatmap visualization of attractor states with color-coded attractors.

    Args:
        df (pd.DataFrame): Attractors dataframe
        gene_cols (list): List of gene column names
        basin_sizes (dict): Basin sizes for each attractor
        output_path (str): Output file path (without extension)
        verbose (bool): Enable verbose output
        svg_output (bool): Also save SVG format

    The heatmap shows the attractor states and nothing else: the trajectory's
    destination is marked on the network and trajectory figures, where it is the
    point of the plot. Here it only added a blue box and a legend row to a panel
    whose job is to let the states be read off and compared.
    """
    log_message("Creating attractor states heatmap with color-coded attractors...", verbose)
    
    # Prepare data for heatmap
    attractors = df['attractor_id'].unique()
    n_attractors = len(attractors)
    n_genes = len(gene_cols)
    
    # Generate colors for attractors
    attractor_colors = generate_attractor_colors(n_attractors)
    attractor_color_map = {att_id: attractor_colors[i] for i, att_id in enumerate(sorted(attractors))}
    
    # Prepare data
    heatmap_data = []
    row_labels = []
    attractor_assignments = []
    attractor_boundaries = []
    current_row = 0
    
    for att_id in sorted(attractors):
        att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
        att_type = att_data.iloc[0]['type']
        cycle_length = att_data.iloc[0]['cycle_length']
        
        for i, (_, row) in enumerate(att_data.iterrows()):
            # Extract gene states
            gene_states = [row[gene] for gene in gene_cols]
            heatmap_data.append(gene_states)
            attractor_assignments.append(att_id)
            
            # Create row label
            if att_type == 'fixed_point':
                label = f"A{att_id} fixed point"
            else:
                if i == 0:  # First row of cycle
                    label = f"A{att_id} limit cycle ({cycle_length})"
                else:  # Subsequent rows
                    label = ""
            row_labels.append(label)
            current_row += 1
        
        attractor_boundaries.append(current_row - 0.5)
    
    # ── Layout: one fixed cell size for every heatmap ────────────────────────────
    # The grid is placed in absolute inches rather than stretched to fill the figure,
    # so a cell is the same square whether the network has one attractor or five.
    # Panels from different runs can then be set side by side and compared directly.
    n_rows = len(heatmap_data)
    grid_w = n_genes * HEATMAP_CELL_INCHES
    grid_h = n_rows * HEATMAP_CELL_INCHES

    # Margins are sized from the text that has to fit inside them
    longest_gene = max((len(g) for g in gene_cols), default=0)
    longest_row = max((len(lbl) for lbl in row_labels), default=0)
    longest_legend = max([len('Inactive gene')] +
                         [len(f'A{a}: {basin_sizes.get(a, 0):,} states (100.0%)')
                          for a in attractors])
    n_legend_rows = 1 + n_attractors

    margin_left = 0.95 + 0.068 * longest_row        # row labels + y axis title
    margin_bottom = 0.70 + 0.050 * longest_gene     # gene names rotated 45 degrees
    margin_top = 0.55                               # title
    margin_right = max(2.3, 0.75 + 0.075 * longest_legend)   # legend panel

    # With few rows the legend is taller than the grid, so the band that holds both
    # is the taller of the two and the grid is centred inside it.
    legend_h = 0.45 + 0.27 * n_legend_rows
    content_h = max(grid_h, legend_h)

    fig_w = margin_left + grid_w + margin_right
    fig_h = margin_bottom + content_h + margin_top

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([margin_left / fig_w,
                       (margin_bottom + (content_h - grid_h) / 2.0) / fig_h,
                       grid_w / fig_w, grid_h / fig_h])
    ax_legend = fig.add_axes([(margin_left + grid_w + 0.2) / fig_w,
                              margin_bottom / fig_h,
                              (margin_right - 0.35) / fig_w, content_h / fig_h])
    
    # Convert to numpy array
    heatmap_matrix = np.array(heatmap_data)
    
    # Create main heatmap (gene states)
    colored_matrix = np.zeros((len(heatmap_data), n_genes, 3))  # RGB matrix

    for i, (gene_states, att_id) in enumerate(zip(heatmap_data, attractor_assignments)):
        attractor_color = attractor_color_map[att_id]
        # Convert hex color to RGB
        if attractor_color.startswith('#'):
            hex_color = attractor_color[1:]
            rgb = tuple(int(hex_color[i:i+2], 16)/255.0 for i in (0, 2, 4))
        else:
            rgb = attractor_color  # If it's already RGB tuple
        
        for j, gene_state in enumerate(gene_states):
            if gene_state == 1:  # Active gene
                colored_matrix[i, j] = rgb
            else:  # Inactive gene
                colored_matrix[i, j] = [0.94, 0.94, 0.94]  # Light gray (#f0f0f0)

    # The axes box is already exactly grid_w x grid_h inches for an n_rows x n_genes
    # matrix, so 'equal' keeps every cell square without resizing the box.
    im = ax.imshow(colored_matrix, aspect='equal', interpolation='nearest')
    
    # Set ticks and labels for main heatmap
    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(gene_cols, rotation=45, ha='right')
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    
    # Add grid to main heatmap
    ax.set_xticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which='minor', color='white', linestyle='-', linewidth=1)
    
    # Add attractor boundaries to main heatmap
    for boundary in attractor_boundaries[:-1]:
        ax.axhline(y=boundary, color='black', linewidth=1, alpha=0.7)
    
    # Create legend in right panel
    ax_legend.axis('off')

    # Create unified legend elements
    legend_elements = []

    # Add inactive state
    legend_elements.append(mpatches.Patch(color='#f0f0f0', label='Inactive gene'))

    # Add each attractor with its color, basin size and percentage
    total_basin_size = sum(basin_sizes.values())

    for att_id in sorted(attractors):
        color = attractor_color_map[att_id]
        basin_size = basin_sizes.get(att_id, 0)
        percentage = (basin_size / total_basin_size) * 100 if total_basin_size > 0 else 0
        att_type = df[df['attractor_id'] == att_id].iloc[0]['type']
        
        label = f'A{att_id}: {basin_size} states ({percentage:.1f}%)'
        legend_elements.append(mpatches.Patch(color=color, label=label))

    # Plot unified legend
    ax_legend.legend(handles=legend_elements, loc='upper left', 
                    bbox_to_anchor=(0, 1), title='Attractors', title_fontsize=12)
    
    # Set title and labels
    ax.set_title('Boolean Network Attractor States', fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Genes', fontsize=12)
    ax.set_ylabel('Attractor States', fontsize=12)

    # Save plot
    plt.savefig(f"{output_path}_heatmap.png", dpi=300, bbox_inches='tight')
    if svg_output:
        plt.savefig(f"{output_path}_heatmap.svg", bbox_inches='tight')

    plt.close()

    if svg_output:
        log_message(f"Heatmap saved as {output_path}_heatmap.png and .svg", verbose)
    else:
        log_message(f"Heatmap saved as {output_path}_heatmap.png", verbose)


# Inches of figure per data unit of layout, and the blank margin around the
# nominal boxes. Fixing the scale is what makes an attractor occupy the same area
# in every figure, whatever else is in it.
DATA_UNIT_INCHES = 0.55
PLOT_PAD = 0.5

# Side of one heatmap cell, in inches. Fixed so the grid never stretches to fill the
# figure: the canvas grows with the number of genes and attractor states instead.
HEATMAP_CELL_INCHES = 0.34


def create_attractor_network(df, gene_cols, basin_sizes, output_path,
                             verbose=False, svg_output=False,
                             target_attractor_id=None, target_label=None,
                             transition_table=None, pred_depth=0, pred_breadth=3,
                             pinned_codes=None, clamp_mask=None):
    """
    Create a network diagram of attractor transitions with color-coded attractors.

    Args:
        df (pd.DataFrame): Attractors dataframe
        gene_cols (list): List of gene column names
        basin_sizes (dict): Basin sizes for each attractor
        output_path (str): Output file path (without extension)
        verbose (bool): Enable verbose output
        svg_output (bool): Also save SVG format
        target_attractor_id: Attractor the system converges to (reported in the log)
        target_label (str|None): Retained for callers; not drawn on this figure
        transition_table: Full 2^N transition table, required for predecessor states
        pred_depth (int): Update steps to walk back from each attractor node (0 = off)
        pred_breadth (int): Predecessors kept per state per step
        pinned_codes (set|None): State codes always kept (the observed samples)
        clamp_mask: Bool array of states consistent with constant (mutation) rules,
            so the basin counts match the ones in the attractors file
    """
    log_message("Creating attractor transition network with color-coded attractors...", verbose)

    # Every state read from the binarized matrix. A royal-blue ring means exactly
    # this and nothing else, wherever it appears: an attractor state that was
    # measured gets one, an upstream state that was measured gets one, and a state
    # that was only inferred never does. The ring used to double as a marker for
    # the attractor the system converges to, which made it unlabelable — the same
    # colour stood for "observed" on one node and "destination" on another.
    observed_codes = set(pinned_codes or ())

    # Create network graph
    G = nx.DiGraph()
    
    # Process each attractor
    attractors = df['attractor_id'].unique()
    n_attractors = len(attractors)
    node_info = {}
    
    # Generate colors for attractors
    attractor_colors = generate_attractor_colors(n_attractors)
    attractor_color_map = {att_id: attractor_colors[i] for i, att_id in enumerate(sorted(attractors))}
    
    for att_id in attractors:
        att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
        att_type = att_data.iloc[0]['type']
        cycle_length = att_data.iloc[0]['cycle_length']
        
        if att_type == 'fixed_point':
            # Single node for fixed point
            node_id = f"A{att_id}"
            G.add_node(node_id)
            
            # Get active genes
            gene_states = {gene: att_data.iloc[0][gene] for gene in gene_cols}
            active_genes = [gene for gene, state in gene_states.items() if state == 1]
            
            node_info[node_id] = {
                'attractor_id': att_id,
                'type': 'fixed_point',
                'active_genes': active_genes,
                'basin_size': basin_sizes.get(att_id, 0),
                'step': 1,
                'cycle_length': 1,
                # Kept so the node can be checked against the observed states
                'code': state_to_code([bool(att_data.iloc[0][g]) for g in gene_cols],
                                      len(gene_cols))
            }
            
            # Self-loop for fixed point
            G.add_edge(node_id, node_id)
            
        else:
            # Multiple nodes for cycle
            cycle_nodes = []
            for _, row in att_data.iterrows():
                step = row['step_in_cycle']
                node_id = f"A{att_id}_S{step}"
                G.add_node(node_id)
                cycle_nodes.append(node_id)
                
                # Get active genes for this step
                gene_states = {gene: row[gene] for gene in gene_cols}
                active_genes = [gene for gene, state in gene_states.items() if state == 1]
                
                node_info[node_id] = {
                    'attractor_id': att_id,
                    'type': 'cycle',
                    'active_genes': active_genes,
                    'basin_size': basin_sizes.get(att_id, 0),  # Same for all nodes in cycle
                    'step': step,
                    'cycle_length': cycle_length,
                    'code': state_to_code([bool(row[g]) for g in gene_cols],
                                          len(gene_cols))
                }
            
            # Add edges for cycle
            for i in range(len(cycle_nodes)):
                current_node = cycle_nodes[i]
                next_node = cycle_nodes[(i + 1) % len(cycle_nodes)]
                G.add_edge(current_node, next_node)
    
    # Calculate layout with special handling for cycles
    pos = {}
    x_offset = 0

    for att_id in sorted(attractors):
        att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
        att_type = att_data.iloc[0]['type']
        cycle_length = att_data.iloc[0]['cycle_length']
        
        if att_type == 'fixed_point':
            node_id = f"A{att_id}"
            pos[node_id] = (x_offset, 0)
            x_offset += 2  # Minimal horizontal separation
            
        else:
            # Arrange cycle nodes in a circle
            radius = max(0.8, cycle_length * 0.15)  # Smaller radius
            
            for i, (_, row) in enumerate(att_data.iterrows()):
                step = row['step_in_cycle']
                node_id = f"A{att_id}_S{step}"
                
                angle = 2 * np.pi * i / cycle_length - np.pi/2
                x = x_offset + radius * np.cos(angle)
                y = radius * np.sin(angle)
                pos[node_id] = (x, y)
            
            x_offset += (radius * 2.5)  # Move to next position

    # ---- Predecessor states -------------------------------------------------
    # Laid out as a radial tree around each attractor: layer d sits on a ring of
    # radius base + d*RING, and every node's angular wedge is subdivided among its
    # own predecessors, so children stay next to the state they feed.
    # Every attractor is laid out inside a box of its own, and these are its
    # bounds. The axes are set from them rather than from where the nodes happened
    # to land, so an attractor takes up the same space whether it is alone in the
    # figure or one of five — which is what makes two runs comparable side by side.
    pred_bounds = None
    pred_nodes = {}
    # Placeholders standing for the predecessors that exist but were not drawn.
    # Without them the figure would imply that the handful of states shown is all
    # there is, when a single attractor can have thousands.
    ghost_nodes = {}
    pred_stats = {}
    if pred_depth and transition_table is not None:
        RING = 1.5
        pinned_codes = observed_codes
        n_genes = len(gene_cols)

        # Redo the horizontal placement: the attractors now need room for their rings.
        attractor_slots = {}
        boxes = []
        x_offset = 0
        for att_id in sorted(attractors):
            att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
            att_type = att_data.iloc[0]['type']
            cycle_length = int(att_data.iloc[0]['cycle_length'])
            # Predecessor rings stretch the drawing to tens of data units while the
            # figure width is capped, so a node takes up far more data space than it
            # does in the plain diagram and the old cycle radius let the states of a
            # long cycle overlap. Size the ring by arc length per state instead.
            base_radius = (0.0 if att_type == 'fixed_point'
                           else max(0.9, cycle_length * 1.6 / (2 * np.pi)))
            # A receiver with no drawn predecessors puts its hollow node one ring
            # past itself, so the slot has to allow for that extra ring or the node
            # drifts into the neighbouring attractor's space.
            reach = base_radius + (pred_depth + 1.15) * RING
            attractor_slots[att_id] = (x_offset + reach, base_radius)
            boxes.append((x_offset, x_offset + 2 * reach, reach))
            x_offset += 2 * reach + 2.5

        if boxes:
            half = max(b[2] for b in boxes)
            pred_bounds = (min(b[0] for b in boxes), max(b[1] for b in boxes),
                           -half, half)

        for att_id in sorted(attractors):
            centre_x, base_radius = attractor_slots[att_id]
            att_data = df[df['attractor_id'] == att_id].sort_values('step_in_cycle')
            att_type = att_data.iloc[0]['type']
            cycle_length = int(att_data.iloc[0]['cycle_length'])

            # Reposition the attractor nodes around their own centre and give each
            # one the angular wedge its predecessors will grow into.
            seeds = []
            spans = {}
            if att_type == 'fixed_point':
                node_id = f"A{att_id}"
                pos[node_id] = (centre_x, 0.0)
                code = state_to_code([bool(att_data.iloc[0][g]) for g in gene_cols], n_genes)
                seeds.append((node_id, code))
                spans[node_id] = (0.0, np.pi)            # the whole circle
            else:
                for i, (_, row) in enumerate(att_data.iterrows()):
                    node_id = f"A{att_id}_S{row['step_in_cycle']}"
                    angle = 2 * np.pi * i / cycle_length - np.pi / 2
                    pos[node_id] = (centre_x + base_radius * np.cos(angle),
                                    base_radius * np.sin(angle))
                    code = state_to_code([bool(row[g]) for g in gene_cols], n_genes)
                    seeds.append((node_id, code))
                    spans[node_id] = (angle, np.pi / cycle_length)

            layers, stats = build_predecessor_layers(
                transition_table, [c for _, c in seeds],
                pred_depth, pred_breadth, pinned=pinned_codes, verbose=verbose,
                valid_mask=clamp_mask)

            # The true layer sizes, so the figure can say how much of the basin it
            # is leaving out. These partition the basin and therefore reconcile
            # with the basin size in the legend; the counts inside `stats` do not,
            # because they only ever expanded the few states that were kept.
            true_layers, reachable = basin_layer_sizes(
                transition_table, [c for _, c in seeds], transition_table.size,
                valid_mask=clamp_mask)
            stats['basin_layers'] = true_layers
            stats['basin_total'] = reachable
            stats['attractor_states'] = len(seeds)
            pred_stats[att_id] = stats

            # Which drawn states the rest of the basin flows through. Computed
            # before anything is positioned, because a receiver needs a slot of its
            # own in the fan and that changes how the fan is spread.
            #
            # Every undrawn state has exactly one first drawn state on its forward
            # path, so these counts partition the remainder: they sum to it, and
            # remainder + drawn == basin size.
            seed_codes = {c for _, c in seeds}
            drawn_codes = [c for L in layers for c, _ in L['chosen']]
            entries = entry_point_counts(transition_table,
                                         list(seed_codes) + drawn_codes,
                                         valid_mask=clamp_mask)
            receivers = {code: n for code, n in entries.items() if n > 0}
            stats['entry_points'] = receivers

            def add_hollow(code, count, parent_node, x, y, angle):
                """A hollow node standing for the states that enter at parent_node."""
                node_id = f"B{att_id}_{code}"
                G.add_node(node_id)
                G.add_edge(node_id, parent_node)
                pos[node_id] = (x, y)
                ghost_nodes[node_id] = {
                    'count': count, 'attractor_id': att_id, 'angle': angle,
                    'basin_remainder': True, 'centre_x': centre_x,
                }

            code_to_node = {c: n for n, c in seeds}
            placed_hollow = set()
            for depth_index, layer in enumerate(layers, start=1):
                by_parent = {}
                for code, successor in layer['chosen']:
                    by_parent.setdefault(successor, []).append(code)

                radius = base_radius + depth_index * RING
                for successor, children in by_parent.items():
                    parent_node = code_to_node.get(successor)
                    if parent_node is None:
                        continue

                    # A receiver gets a slot beside its own children rather than a
                    # position found by searching for free space. Anything placed
                    # outside the fan has to reach its target across the rings, and
                    # on a cycle every such line crosses the fan it came from.
                    wants_hollow = (successor in receivers
                                    and successor not in placed_hollow)
                    slots = len(children) + (1 if wants_hollow else 0)
                    if not slots:
                        continue

                    angle_centre, half = spans[parent_node]
                    # Subdividing the parent's wedge by angle alone makes a handful of
                    # children fan across it. Capping the spread by arc length instead
                    # keeps each child visibly attached to the state it feeds.
                    spread = min(2 * half, slots * 1.15 / max(radius, 1e-6))
                    step = spread / slots
                    for k, code in enumerate(sorted(children)):
                        child_angle = angle_centre - spread / 2 + step * (k + 0.5)
                        node_id = f"P{att_id}_D{depth_index}_{code}"
                        G.add_node(node_id)
                        G.add_edge(node_id, parent_node)
                        pos[node_id] = (centre_x + radius * np.cos(child_angle),
                                        radius * np.sin(child_angle))
                        spans[node_id] = (child_angle, step / 2)
                        code_to_node[code] = node_id
                        pred_nodes[node_id] = {
                            'code': code, 'depth': depth_index,
                            'attractor_id': att_id, 'pinned': code in pinned_codes,
                        }

                    if wants_hollow:
                        hollow_angle = angle_centre - spread / 2 + step * (slots - 0.5)
                        add_hollow(successor, receivers[successor], parent_node,
                                   centre_x + radius * np.cos(hollow_angle),
                                   radius * np.sin(hollow_angle), hollow_angle)
                        placed_hollow.add(successor)

            # Receivers with no drawn predecessors of their own never came up in the
            # loop above, so they get a slot one ring further out, still inside their
            # own wedge — the line stays short and crosses nothing.
            for code, count in receivers.items():
                if code in placed_hollow:
                    continue
                parent_node = code_to_node.get(code)
                if parent_node is None:
                    continue
                px_, py_ = pos[parent_node]
                node_radius = float(np.hypot(px_ - centre_x, py_))
                angle_centre, _ = spans.get(parent_node,
                                            (np.arctan2(py_, px_ - centre_x), 0.3))
                out = node_radius + RING
                add_hollow(code, count, parent_node,
                           centre_x + out * np.cos(angle_centre),
                           out * np.sin(angle_centre), angle_centre)
                placed_hollow.add(code)


        for att_id, stats in pred_stats.items():
            for d, (avail, drawn) in enumerate(zip(stats['available'],
                                                   stats['drawn']), start=1):
                log_message(f"  A{att_id} step {d} back: {drawn} drawn of "
                            f"{avail:,} predecessors of the states already drawn",
                            verbose)
            # The unconditioned layers, which are the ones that sum to the basin.
            shown = ' · '.join(f"{n:,} at {d}" for d, n
                               in enumerate(stats.get('basin_layers', []), start=1))
            log_message(f"  A{att_id} basin: {stats.get('basin_total', 0):,} states "
                        f"({stats.get('attractor_states', 0)} in the attractor) — "
                        f"steps back: {shown}", verbose)

        # The drawn states are the figure's claim, so they are written out too.
        pred_file = f"{output_path}_network_predecessors.tsv"
        with open(pred_file, 'w', encoding='utf-8') as handle:
            handle.write("node\tattractor_id\tsteps_upstream\tobserved_sample\t"
                         "binary_state\t" + "\t".join(gene_cols) + "\n")
            for node_id, info in sorted(pred_nodes.items(),
                                        key=lambda kv: (kv[1]['attractor_id'],
                                                        kv[1]['depth'], kv[0])):
                bits = code_to_state(info['code'], n_genes)
                handle.write(f"{node_id}\t{info['attractor_id']}\t{info['depth']}\t"
                             f"{'yes' if info['pinned'] else 'no'}\t"
                             + ''.join('1' if b else '0' for b in bits) + "\t"
                             + "\t".join('1' if b else '0' for b in bits) + "\n")
        log_message(f"Predecessor states written to {pred_file}", verbose)

    # Debug: print positions
    if verbose:
        print("Node positions:")
        for node, position in pos.items():
            print(f"  {node}: {position}")
    
    # Calculate plot boundaries based on actual positions
    if pos:
        x_coords = [p[0] for p in pos.values()]
        y_coords = [p[1] for p in pos.values()]
        x_range = max(x_coords) - min(x_coords)
        y_range = max(y_coords) - min(y_coords)
        
        # Set figure size based on content
        if pred_bounds:
            # A fixed number of inches per data unit, applied to the nominal boxes
            # rather than to wherever the nodes ended up. Deriving the width and
            # the height independently and then locking the aspect to equal let
            # the scale drift with the number of attractors and with which way a
            # fan happened to point: one attractor came out tall and narrow with
            # small nodes, two came out wide and short. At a fixed scale an
            # attractor is the same size in every figure, so panels from different
            # runs can sit side by side.
            bx0, bx1, by0, by1 = pred_bounds
            plot_w = (bx1 - bx0) + 2 * PLOT_PAD
            plot_h = (by1 - by0) + 2 * PLOT_PAD
            fig_width = min(40, plot_w * DATA_UNIT_INCHES)
            fig_height = min(40, plot_h * DATA_UNIT_INCHES)
        else:
            cap_w, cap_h = (16, 12)
            fig_width = min(cap_w, max(8, x_range + 2))
            fig_height = min(cap_h, max(6, y_range + 2))
        
        plt.figure(figsize=(fig_width, fig_height))
    else:
        plt.figure(figsize=(10, 8))
    
    node_sizes = []
    node_colors = []
    node_edge_colors = []
    drew_observed = False      # whether any blue ring made it onto the canvas
    node_linewidths = []

    fixed_node_size = 800

    # Node area is in points, which are physical, while the layout is in data
    # units. Deriving the conversion from the figure's own extent made a node a
    # few hundredths of the plot in one figure and half a ring in another. At a
    # fixed scale it follows directly: a node is always the same fraction of the
    # ring spacing, and therefore the same physical size in every figure.
    if pred_bounds:
        radius_points = 0.30 * 1.5 * DATA_UNIT_INCHES * 72.0
        fixed_node_size = float(np.pi * radius_points ** 2)

    node_alphas = []

    for node_id in G.nodes():
        if node_id in ghost_nodes:
            # Hollow and dashed so it reads as "more of these" rather than as a
            # state in its own right. Sized by its count on a log scale, since the
            # counts span three orders of magnitude across one figure.
            info = ghost_nodes[node_id]
            # Outlined in its attractor's colour, and sized by how many states it
            # stands for on a log scale, since one figure mixes counts of tens and
            # counts of tens of thousands.
            node_sizes.append(150 + 105 * np.log10(max(info['count'], 1) + 1))
            node_edge_colors.append(attractor_color_map[info['attractor_id']])
            node_linewidths.append(1.8)
            node_colors.append('white')
            node_alphas.append(1.0)
            continue

        if node_id in pred_nodes:
            # Upstream states: same hue as the attractor they drain into, but smaller
            # and fainter the further back they sit, so the attractor stays dominant.
            info = pred_nodes[node_id]
            att_id = info['attractor_id']
            node_sizes.append(max(90, fixed_node_size * (0.42 ** info['depth'])))
            node_colors.append(attractor_color_map[att_id])
            node_alphas.append(max(0.30, 0.75 - 0.15 * (info['depth'] - 1)))
            if info['pinned']:
                node_edge_colors.append('royalblue')   # an observed sample
                node_linewidths.append(2.5)
                drew_observed = True
            else:
                node_edge_colors.append('dimgray')
                node_linewidths.append(0.8)
            continue

        info = node_info[node_id]
        att_id = info['attractor_id']
        node_sizes.append(fixed_node_size)
        node_colors.append(attractor_color_map[att_id])
        node_alphas.append(0.9)
        if info.get('code') in observed_codes:
            node_edge_colors.append('royalblue')
            node_linewidths.append(3.5)
            drew_observed = True
        else:
            node_edge_colors.append('black')
            node_linewidths.append(1)

    # Draw the network
    node_order = list(G.nodes())
    size_of = dict(zip(node_order, node_sizes))
    width_of = dict(zip(node_order, node_linewidths))

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                           alpha=node_alphas if (pred_nodes or ghost_nodes) else 0.8,
                           edgecolors=node_edge_colors, linewidths=node_linewidths)

    from matplotlib.patches import FancyArrowPatch

    # Length of the separately drawn arrowhead on a dotted route, in points. The
    # connector stops this far short of the rim so the two meet instead of overlap.
    GHOST_HEAD_POINTS = 11.0

    def margin_for(node_id):
        """
        How far an arrow must stop short of a node, in points.

        networkx takes one margin for a whole edgelist, but the nodes here differ
        by a factor of six in area: a margin that clears a depth-2 state leaves the
        arrowhead buried inside the attractor. So edges are grouped by the size of
        the node they point at and each group gets its own margin.
        """
        radius = float(np.sqrt(max(size_of.get(node_id, 300), 1) / np.pi))
        return radius + width_of.get(node_id, 1.0) / 2.0 + 3.0

    def draw_edges(edgelist, extra_target=0.0, **kwargs):
        """
        Draw edges in groups that share a target margin.

        extra_target reserves room at the far end for an arrowhead drawn
        separately; without it the connector runs the whole way to the rim, under
        the head and out through its tip.
        """
        groups = {}
        for u, v in edgelist:
            key = (round(margin_for(u)), round(margin_for(v) + extra_target))
            groups.setdefault(key, []).append((u, v))
        for (source_margin, target_margin), edges in groups.items():
            nx.draw_networkx_edges(G, pos, edgelist=edges,
                                   min_source_margin=source_margin,
                                   min_target_margin=target_margin, **kwargs)

    # Separate self-loops from regular edges
    self_loops = [(u, v) for u, v in G.edges() if u == v]
    regular_edges = [(u, v) for u, v in G.edges()
                     if u != v and u not in pred_nodes and u not in ghost_nodes]
    ghost_edges = [(u, v) for u, v in G.edges() if u in ghost_nodes]
    # Upstream edges carry the same meaning but must not compete with the attractor
    # for attention, so they are drawn thin and grey underneath it.
    pred_edges = [(u, v) for u, v in G.edges() if u != v and u in pred_nodes]

    if ghost_edges:
        # The connector only. A dotted linestyle applies to the arrowhead too and
        # breaks it into chevrons, so the head is drawn separately and solid once
        # the axes scale is known.
        # arrows=True even though this draws no head: with arrows=False networkx
        # falls back to a LineCollection, which ignores min_target_margin outright
        # and runs the dotted line from centre to centre, straight through the head
        # drawn below. arrowstyle='-' keeps it headless while honouring the margins.
        draw_edges(ghost_edges, edge_color='gray', arrows=True, arrowstyle='-',
                   width=1.0, alpha=0.6, style='dotted',
                   extra_target=GHOST_HEAD_POINTS)

    if pred_edges:
        draw_edges(pred_edges, edge_color='gray', arrows=True, arrowsize=11,
                   arrowstyle='-|>', width=0.9, alpha=0.6)

    if regular_edges:
        # A filled head at the old size swamps the short hop between two states of
        # a cycle, which can be barely longer than the heads at each end.
        draw_edges(regular_edges, edge_color='black', arrows=True, arrowsize=10,
                   arrowstyle='-|>', width=1.8, alpha=0.75)

    # The count is the whole point of a placeholder, so it is written on it.
    # Placed just outside the node along its own radius, which keeps it clear of
    # the ring it sits on.
    for node_id, info in ghost_nodes.items():
        x, y = pos[node_id]
        angle = info['angle']
        # Written just beyond the node, continuing the direction it was placed in,
        # so the label never falls back onto the state it hangs off.
        plt.text(x + 0.55 * np.cos(angle), y + 0.55 * np.sin(angle),
                 f"{info['count']:,}", fontsize=8, color='dimgray',
                 ha='center', va='center', zorder=6,
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                           edgecolor='none', alpha=0.85))

    # The legend names the attractors and the blue ring. The rest of what the
    # figure draws — upstream states, hollow entry-point nodes — is labelled on the
    # canvas itself or belongs in the caption, and spelling it all out here crowded
    # the plot more than it explained it.
    legend_elements = []
    total_basin = sum(basin_sizes.values())

    for att_id in sorted(attractors):
        basin_size = basin_sizes.get(att_id, 0)
        percentage = (basin_size / total_basin) * 100 if total_basin > 0 else 0
        att_type = str(df[df['attractor_id'] == att_id].iloc[0]['type']).replace('_', ' ')
        legend_elements.append(
            mpatches.Patch(color=attractor_color_map[att_id],
                           label=f'Attractor {att_id} ({att_type}): '
                                 f'{basin_size:,} states ({percentage:.1f}%)')
        )

    # Only claimed when a ring was actually drawn: with no -b there are no observed
    # states, and a legend row for a symbol absent from the canvas is noise.
    if drew_observed:
        legend_elements.append(
            Line2D([0], [0], marker='o', linestyle='none',
                   markerfacecolor='none', markeredgecolor='royalblue',
                   markeredgewidth=2.2, markersize=11,
                   label='States in the binarized matrix')
        )

    if target_attractor_id is not None:
        log_message(f"Target attractor is {target_attractor_id}", verbose)

    if pred_nodes or ghost_nodes:
        # Moved out from under the drawing. The hollow nodes are placed by
        # clearance from other NODES, which cannot see the legend box, so a legend
        # sitting inside the axes ends up with a node on top of it.
        plt.legend(handles=legend_elements, loc='upper center',
                   bbox_to_anchor=(0.5, -0.01), frameon=True, fontsize=9)
    else:
        plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    plt.title('Boolean Network Attractors', fontsize=14, fontweight='bold')
    plt.axis('off')
    if pred_nodes or ghost_nodes:
        # The layers are concentric rings; without an equal aspect they render as
        # ellipses and "one step back" looks like a different distance up than across.
        plt.gca().set_aspect('equal', adjustable='box')
    
    # Set axis limits to fit content tightly
    if pos:
        x_coords = [p[0] for p in pos.values()]
        y_coords = [p[1] for p in pos.values()]
        if pred_bounds:
            bx0, bx1, by0, by1 = pred_bounds
            plt.xlim(bx0 - PLOT_PAD, bx1 + PLOT_PAD)
            plt.ylim(by0 - PLOT_PAD, by1 + PLOT_PAD)
        else:
            margin = 1
            plt.xlim(min(x_coords) - margin, max(x_coords) + margin)
            # No headroom needed: with predecessors the legend sits below the axes.
            plt.ylim(min(y_coords) - margin, max(y_coords) + margin)
    
    plt.tight_layout()

    # Self-loops last: their size is a node radius, which is given in points, and
    # converting that to data units needs the axes scale — which only exists once
    # the limits and the aspect are settled. Sizing them off a hardcoded divisor
    # beforehand made the loop vanish on wide figures and swamp the node on narrow
    # ones.
    if self_loops:
        ax = plt.gca()
        fig = plt.gcf()
        fig.canvas.draw()
        box = ax.get_window_extent()
        x_lo, x_hi = ax.get_xlim()
        points_per_data = (box.width * 72.0 / fig.dpi) / max(x_hi - x_lo, 1e-9)
        data_per_point = 1.0 / max(points_per_data, 1e-9)

        for u, _ in self_loops:
            x, y = pos[u]
            radius = float(np.sqrt(max(size_of.get(u, 800), 1) / np.pi)) * data_per_point
            # Start and end on the node's own rim, a sixth of a turn apart, with the
            # arc bulging up and away. A negative rad curls the loop back inside the
            # node, where it reads as a scribble rather than a return to self.
            # Both ends sit on the node's rim near the top, and the arc bulges up
            # and over between them. The sign of rad is what decides whether it
            # goes over the node or dips down inside it, where it reads as a
            # scribble; checked against all four combinations before settling here.
            a0, a1 = np.deg2rad(65), np.deg2rad(115)
            start = (x + radius * np.cos(a0), y + radius * np.sin(a0))
            finish = (x + radius * np.cos(a1), y + radius * np.sin(a1))
            ax.add_patch(FancyArrowPatch(
                start, finish, connectionstyle="arc3,rad=1.9",
                arrowstyle='-|>', mutation_scale=13, color='black',
                linewidth=1.8, alpha=0.85, zorder=4))

    # Solid heads for the dotted routes, now that points convert to data units.
    if ghost_edges:
        ax = plt.gca()
        fig = plt.gcf()
        fig.canvas.draw()
        box = ax.get_window_extent()
        x_lo, x_hi = ax.get_xlim()
        points_per_data = (box.width * 72.0 / fig.dpi) / max(x_hi - x_lo, 1e-9)
        data_per_point = 1.0 / max(points_per_data, 1e-9)

        for u, v in ghost_edges:
            ux, uy = pos[u]
            vx, vy = pos[v]
            length = float(np.hypot(vx - ux, vy - uy))
            if length <= 0:
                continue
            dx, dy = (vx - ux) / length, (vy - uy) / length
            stop = margin_for(v) * data_per_point          # the target's rim
            head = GHOST_HEAD_POINTS * data_per_point
            tip = (vx - dx * stop, vy - dy * stop)
            tail = (vx - dx * (stop + head), vy - dy * (stop + head))
            ax.add_patch(FancyArrowPatch(
                tail, tip, arrowstyle='-|>', mutation_scale=11,
                color='gray', linewidth=1.0, alpha=0.75, zorder=4))

    # Save plot
    plt.savefig(f"{output_path}_network.png", dpi=300, bbox_inches='tight')
    if svg_output:
        plt.savefig(f"{output_path}_network.svg", bbox_inches='tight')

    plt.close()

    if svg_output:
        log_message(f"Network diagram saved as {output_path}_network.png and .svg", verbose)
    else:
        log_message(f"Network diagram saved as {output_path}_network.png", verbose)


def visualize_attractors(args):
    """
    Main function to visualize Boolean network attractors
    
    Args:
        args: Parsed command line arguments
    """
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(args.input_file)
    
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read attractors data
        df, gene_cols = read_attractors_file(args.input_file, args.verbose)

        # Get basin sizes from dataframe
        basin_sizes = get_basin_sizes_from_df(df, args.verbose)

        # Create output base path
        input_path = Path(args.input_file)
        if args.output_base:
            output_base = os.path.join(output_dir, args.output_base)
        else:
            mutation_suffix = extract_mutation_suffix(input_path.name)
            if mutation_suffix:
                output_base = os.path.join(output_dir, f"attractors_visualization{mutation_suffix}")
            else:
                output_base = os.path.join(output_dir, "attractors_visualization")

        # Resolve target attractor from binarized matrix (optional)
        target_attractor_id = None
        target_label = None
        gene_rules = None
        if args.binarized_matrix:
            final_state, sample_name = read_final_matrix_state(
                args.binarized_matrix, gene_cols, args.verbose
            )
            if args.rules_file:
                gene_rules = read_rules_simple(args.rules_file, args.verbose)
            attractors_dict = build_attractors_dict(df, gene_cols)
            target_attractor_id = find_target_attractor_id(
                final_state, attractors_dict, gene_cols,
                gene_rules=gene_rules, verbose=args.verbose
            )
            if target_attractor_id is not None:
                target_label = f"{sample_name} → Attractor {target_attractor_id}"
            else:
                print("Warning: could not determine target attractor for final matrix state.", file=sys.stderr)

        # Predecessor states need the rules: the attractors file records where the
        # system ends up, not how any state maps onto the next one.
        transition_table = None
        clamp_mask = None
        pinned_codes = set()
        pred_depth = max(0, int(getattr(args, 'predecessors', 0) or 0))
        if pred_depth:
            if gene_rules is None:
                if not args.rules_file:
                    print("ERROR: --predecessors needs -r/--rules_file to know the "
                          "network's update function.", file=sys.stderr)
                    sys.exit(1)
                gene_rules = read_rules_simple(args.rules_file, args.verbose)
            transition_table = build_transition_table(gene_cols, gene_rules, args.verbose)
            # A constant rule is how a knockout or overexpression is encoded, and
            # step 1 counts only the states that respect it.
            clamp_mask = clamped_state_mask(gene_cols, gene_rules,
                                            transition_table.size)
            if clamp_mask is not None:
                log_message(f"Constant rules clamp the state space to "
                            f"{int(clamp_mask.sum()):,} of {transition_table.size:,} "
                            f"states", args.verbose)


        # Every row of the binarized matrix. These states are measurements rather
        # than possibilities, so they outrank the Hamming selection when upstream
        # states are chosen, and they are the ones the figure rings in blue — which
        # it can do whether or not predecessors are drawn.
        if args.binarized_matrix:
            matrix_df = pd.read_csv(args.binarized_matrix, sep='\t')
            first_col = matrix_df.columns[0]
            if not pd.api.types.is_numeric_dtype(matrix_df[first_col]):
                matrix_df = matrix_df.drop(columns=[first_col])
            matrix_df.columns = sanitize_gene_names(
                matrix_df.columns, os.path.basename(args.binarized_matrix))
            present = set(matrix_df.columns)
            for _, row in matrix_df.iterrows():
                pinned_codes.add(state_to_code(
                    [bool(int(row[g])) if g in present else False for g in gene_cols],
                    len(gene_cols)))
            log_message(f"{len(pinned_codes)} distinct state(s) read from the "
                        f"binarized matrix", args.verbose)

        # Generate visualizations
        if args.heatmap:
            create_attractor_heatmap(df, gene_cols, basin_sizes, output_base, args.verbose, args.svg)

        if args.network:
            create_attractor_network(df, gene_cols, basin_sizes, output_base, args.verbose, args.svg,
                                     target_attractor_id=target_attractor_id, target_label=target_label,
                                     transition_table=transition_table,
                                     pred_depth=pred_depth,
                                     pred_breadth=max(1, int(args.predecessors_per_state)),
                                     pinned_codes=pinned_codes,
                                     clamp_mask=clamp_mask)
        
        # Print summary
        print(f"\n{'='*60}")
        print("VISUALIZATION COMPLETED SUCCESSFULLY")
        print('='*60)
        print(f"Input file: {args.input_file}")
        print(f"Output directory: {output_dir}")
        
        if args.heatmap:
            if args.svg:
                print(f"- {os.path.basename(output_base)}_heatmap.png/svg (color-coded states heatmap)")
            else:
                print(f"- {os.path.basename(output_base)}_heatmap.png (color-coded states heatmap)")
        if args.network:
            if args.svg:
                print(f"- {os.path.basename(output_base)}_network.png/svg (color-coded transition network)")
            else:
                print(f"- {os.path.basename(output_base)}_network.png (color-coded transition network)")
        
        print(f"Attractors analyzed: {len(df['attractor_id'].unique())}")
        print(f"Total states: {len(df)}")
        print(f"Genes: {len(gene_cols)}")
        # Show mutation info if detected
        mutation_suffix = extract_mutation_suffix(Path(args.input_file).name)
        if mutation_suffix:
            print(f"Mutations detected: {mutation_suffix[1:].replace('_', '=', 1).replace('_', ', ')}")
        print('='*60)
        
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Visualize Boolean network attractors from analysis results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 3.BNI3_Visualize_Attractors.py -i attractors.tsv
  python3 3.BNI3_Visualize_Attractors.py -i attractors/attractors.tsv --heatmap --network -v
  python3 3.BNI3_Visualize_Attractors.py -i attractors.tsv --svg -v
  python3 3.BNI3_Visualize_Attractors.py -i attractors.tsv --network-only -ob my_viz

Visualization Types:
  --heatmap    : Create heatmap showing gene states across attractors (default: enabled)
  --network    : Create network diagram showing attractor transitions (default: enabled)
  --both       : Generate both visualizations (default behavior)

Output Options:
  -o           : Specify output directory
  -ob          : Custom base name for output files

Notes:
  - Input file should be the attractors.tsv file generated by 1.BNI3_Attractors.py
  - Generates PNG format by default, add --svg for vector graphics
  - Heatmap: Each row is an attractor state, columns are genes, colors show activity
  - Network: Nodes are states, edges show transitions, size reflects basin size
  - Each attractor has a distinctive color for easy identification
  - Basin size information is read directly from the TSV file
  - For large networks (>20 genes), heatmap is more readable than network diagram
        """
    )
    
    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input_file', type=str, required=True,
                         help='Input attractors TSV file (from 1.BNI3_Attractors.py)')
    
    # Optional parameters
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-o', '--output_dir', type=str, default=None,
                         help='Output directory (default: same as input file)')
    optional.add_argument('-ob', '--output_base', type=str, default=None,
                         help='Base name for output files (default: auto-generated)')
    optional.add_argument('-b', '--binarized_matrix', type=str, default=None,
                          help='Binarized expression matrix TSV — last row is used to mark '
                               'the target attractor in both heatmap and network')
    optional.add_argument('-r', '--rules_file', type=str, default=None,
                          help='Rules file (selected_rules.tsv or full table) — required only '
                               'if the final matrix state is not already an attractor state')
    
    # Visualization options
    viz_group = parser.add_argument_group('Visualization options')
    viz_group.add_argument('--heatmap', action='store_true', default=False,
                          help='Generate heatmap visualization')
    viz_group.add_argument('--network', action='store_true', default=False,
                          help='Generate network diagram visualization')
    viz_group.add_argument('--heatmap-only', action='store_true',
                          help='Generate only heatmap (shortcut for --heatmap)')
    viz_group.add_argument('--network-only', action='store_true',
                          help='Generate only network diagram (shortcut for --network)')
    viz_group.add_argument('-v', '--verbose', action='store_true',
                          help='Show detailed processing information')
    viz_group.add_argument('--predecessors', type=int, default=0, metavar='STEPS',
                           help='Also draw states that lead INTO each attractor, walking '
                                'this many update steps back (0 = off, 2 is a good start). '
                                'Requires -r/--rules_file')
    viz_group.add_argument('--predecessors-per-state', type=int, default=3, metavar='K',
                           help='How many predecessors to keep per state per step '
                                '(default: 3). A 15-gene attractor can have thousands, '
                                'so the closest ones in Hamming distance are kept; any '
                                'state observed in -b/--binarized_matrix is always kept')
    viz_group.add_argument('--svg', action='store_true',
                      help='Also generate SVG format (vector graphics)')
    
    parser.add_argument('--version', action='version', version='Boolean Network Attractor Visualizer v1.0')
    
    args = parser.parse_args()
    
    # Handle visualization options
    if args.heatmap_only:
        args.heatmap = True
        args.network = False
    elif args.network_only:
        args.heatmap = False
        args.network = True
    elif not args.heatmap and not args.network:
        # Default: generate both
        args.heatmap = True
        args.network = True
    
    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Run visualization
    try:
        visualize_attractors(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()