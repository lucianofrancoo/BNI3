#!/usr/bin/env python3
"""
Boolean Network Attractor Visualizer
Creates visualizations of Boolean network attractors from attractor analysis results.
Part of the Boolean Network Inference (BNI) pipeline.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
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

    matrix_genes = list(df.columns)
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
                             verbose=False, svg_output=False,
                             target_attractor_id=None, target_label=None):
    """
    Create a heatmap visualization of attractor states with color-coded attractors.

    Args:
        df (pd.DataFrame): Attractors dataframe
        gene_cols (list): List of gene column names
        basin_sizes (dict): Basin sizes for each attractor
        output_path (str): Output file path (without extension)
        verbose (bool): Enable verbose output
        svg_output (bool): Also save SVG format
        target_attractor_id: Attractor the system converges to (highlighted in blue)
        target_label (str|None): Label for the blue marker in the legend
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
    
    # Check if we have only one attractor with one state (single fixed point)
    single_fixed_point = (len(heatmap_data) == 1 and 
                         len(attractors) == 1 and 
                         df.iloc[0]['type'] == 'fixed_point')
    
    # Calculate figure size based on data
    fig_width = max(10, min(22, n_genes * 0.5 + 6))
    
    # Adjust height for single fixed point case
    if single_fixed_point:
        fig_height = max(4, min(8, n_genes * 0.3 + 2))  # Más compacto para punto fijo único
    else:
        fig_height = max(6, min(16, len(df) * 0.3 + 3))
    
    # Create figure with gridspec for main plot and attractor color bar
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[20, 3], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])
    
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

    # Display the colored matrix with aspect control
    if single_fixed_point:
        # For single fixed point, control aspect ratio to make it square-like
        im = ax.imshow(colored_matrix, aspect='equal', interpolation='nearest')
    else:
        # Normal behavior for multiple states/attractors
        im = ax.imshow(colored_matrix, aspect='auto', interpolation='nearest')
    
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
    
    # Blue border: mark the rows of the target attractor (final matrix state destination)
    if target_attractor_id is not None:
        target_rows = [i for i, a in enumerate(attractor_assignments) if a == target_attractor_id]
        if target_rows:
            row_start = min(target_rows)
            row_end = max(target_rows)
            from matplotlib.patches import Rectangle as _Rect
            border = _Rect((-0.5, row_start - 0.5), n_genes, row_end - row_start + 1,
                           linewidth=3, edgecolor='royalblue', facecolor='none', alpha=0.9)
            ax.add_patch(border)
            blue_label = target_label or f"System → Attractor {target_attractor_id}"
            legend_elements.append(
                mpatches.Patch(color='royalblue', alpha=0.9, label=blue_label)
            )
            ax_legend.legend(handles=legend_elements, loc='upper left',
                             bbox_to_anchor=(0, 1), title='Attractors', title_fontsize=12)
            log_message(f"Highlighted target attractor {target_attractor_id} rows {target_rows} in heatmap", verbose)

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


def create_attractor_network(df, gene_cols, basin_sizes, output_path,
                             verbose=False, svg_output=False,
                             target_attractor_id=None, target_label=None):
    """
    Create a network diagram of attractor transitions with color-coded attractors.

    Args:
        df (pd.DataFrame): Attractors dataframe
        gene_cols (list): List of gene column names
        basin_sizes (dict): Basin sizes for each attractor
        output_path (str): Output file path (without extension)
        verbose (bool): Enable verbose output
        svg_output (bool): Also save SVG format
        target_attractor_id: Attractor the system converges to (highlighted in blue)
        target_label (str|None): Label for the blue marker in the legend
    """
    log_message("Creating attractor transition network with color-coded attractors...", verbose)
    
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
                'cycle_length': 1
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
                    'cycle_length': cycle_length
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
        fig_width = min(16, max(8, x_range + 2))
        fig_height = min(12, max(6, y_range + 2))
        
        plt.figure(figsize=(fig_width, fig_height))
    else:
        plt.figure(figsize=(10, 8))
    
    node_sizes = []
    node_colors = []
    node_edge_colors = []
    node_linewidths = []

    fixed_node_size = 800

    for node_id in G.nodes():
        info = node_info[node_id]
        att_id = info['attractor_id']
        node_sizes.append(fixed_node_size)
        node_colors.append(attractor_color_map[att_id])
        if target_attractor_id is not None and att_id == target_attractor_id:
            node_edge_colors.append('royalblue')
            node_linewidths.append(3.5)
        else:
            node_edge_colors.append('black')
            node_linewidths.append(1)

    # Draw the network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                           alpha=0.8, edgecolors=node_edge_colors, linewidths=node_linewidths)

    # Import for fancy arrows
    from matplotlib.patches import FancyArrowPatch

    # Separate self-loops from regular edges
    self_loops = [(u, v) for u, v in G.edges() if u == v]
    regular_edges = [(u, v) for u, v in G.edges() if u != v]

    # Draw regular edges first
    if regular_edges:
        # Calculate dynamic margins based on node sizes
        base_margin = 15
        margin = base_margin + (fixed_node_size / 2000) * 10
        
        nx.draw_networkx_edges(G, pos, edgelist=regular_edges, edge_color='black', 
                            arrows=True, arrowsize=15, arrowstyle='->', 
                            width=2, alpha=0.7,
                            min_source_margin=margin, min_target_margin=margin)

   # Draw self-loops with FancyArrowPatch for better control
    if self_loops:
        # Check if we only have fixed points (no cycles)
        only_fixed_points = all(df[df['attractor_id'] == att_id].iloc[0]['type'] == 'fixed_point' 
                            for att_id in attractors)
        
        # Scale factor for single fixed point vs multiple attractors
        if len(attractors) == 1 and only_fixed_points:
            scale_factor = 0.5  # Much smaller for single fixed point
        elif only_fixed_points:
            scale_factor = 0.5  # Smaller for multiple fixed points
        else:
            scale_factor = 1.0  # Normal size when there are cycles
        
        for u, v in self_loops:
            x, y = pos[u]
            
            # Calculate node radius with scaling
            node_radius = np.sqrt(fixed_node_size / np.pi) / 100
            
            # Position loop above and to the side of node with adaptive sizing
            loop_radius = node_radius * 0.8 * scale_factor
            start_angle = np.pi/4  # 45 degrees
            end_angle = 3*np.pi/4  # 135 degrees
            
            start_x = x + loop_radius * np.cos(start_angle)
            start_y = y + loop_radius * np.sin(start_angle) + node_radius * 0.3 * scale_factor
            end_x = x + loop_radius * np.cos(end_angle)
            end_y = y + loop_radius * np.sin(end_angle) + node_radius * 0.3 * scale_factor
            
            # Create self-loop with FancyArrowPatch
            loop = FancyArrowPatch(
                (start_x, start_y), (end_x, end_y),
                connectionstyle=f"arc3,rad={3 * scale_factor}",
                arrowstyle='->',
                mutation_scale=15 * scale_factor,
                color='black',
                linewidth=2,
                alpha=0.7
            )
            plt.gca().add_patch(loop)
    
    # Create legend with attractor colors and basin information
    legend_elements = []
    total_basin = sum(basin_sizes.values())

    for att_id in sorted(attractors):
        color = attractor_color_map[att_id]
        basin_size = basin_sizes.get(att_id, 0)
        percentage = (basin_size / total_basin) * 100 if total_basin > 0 else 0
        att_type = df[df['attractor_id'] == att_id].iloc[0]['type']
        legend_elements.append(
            mpatches.Patch(color=color, label=f'A{att_id} ({att_type}): {basin_size} states ({percentage:.1f}%)')
        )

    if target_attractor_id is not None:
        blue_label = target_label or f"System → Attractor {target_attractor_id}"
        legend_elements.append(
            mpatches.Patch(color='royalblue', alpha=0.9, label=blue_label)
        )
        log_message(f"Highlighted target attractor {target_attractor_id} node(s) in network", verbose)

    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    plt.title('Boolean Network Attractor Transitions', fontsize=14, fontweight='bold')
    plt.axis('off')
    
    # Set axis limits to fit content tightly
    if pos:
        x_coords = [p[0] for p in pos.values()]
        y_coords = [p[1] for p in pos.values()]
        margin = 1
        plt.xlim(min(x_coords) - margin, max(x_coords) + margin)
        plt.ylim(min(y_coords) - margin, max(y_coords) + margin)
    
    plt.tight_layout()
    
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
        if args.binarized_matrix:
            final_state, sample_name = read_final_matrix_state(
                args.binarized_matrix, gene_cols, args.verbose
            )
            gene_rules = None
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

        # Generate visualizations
        if args.heatmap:
            create_attractor_heatmap(df, gene_cols, basin_sizes, output_base, args.verbose, args.svg,
                                     target_attractor_id=target_attractor_id, target_label=target_label)

        if args.network:
            create_attractor_network(df, gene_cols, basin_sizes, output_base, args.verbose, args.svg,
                                     target_attractor_id=target_attractor_id, target_label=target_label)
        
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