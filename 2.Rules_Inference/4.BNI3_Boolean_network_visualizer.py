#!/usr/bin/env python3
"""
Boolean Network Visualizer
Creates visualizations of Boolean gene regulatory networks from inferred rules.
"""

import argparse
import os
import sys
import re
import json
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def read_rules_table(table_file, verbose):
    """
    Read the rules table from TSV file
    
    Args:
        table_file (str): Path to the rules table file
        verbose (bool): Enable verbose output
        
    Returns:
        pd.DataFrame: Rules dataframe or None if error
    """
    try:
        log_message(f"Reading rules table from: {table_file}", verbose)
        
        if table_file.endswith('.tsv'):
            df = pd.read_csv(table_file, sep='\t', encoding='utf-8')
        else:
            df = pd.read_csv(table_file, encoding='utf-8')
        
        log_message(f"Table loaded successfully with {len(df)} rules", verbose)
        
        # Validate required columns
        required_cols = ['Gene', 'Position', 'Rule']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
            return None
            
        return df
        
    except Exception as e:
        print(f"ERROR reading file {table_file}: {e}", file=sys.stderr)
        return None


def select_rules_by_criteria(df, criteria, custom_positions, verbose):
    """
    Select rules based on Position values only
    
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
            print("ERROR: custom_positions criteria requires --positions parameter", file=sys.stderr)
            sys.exit(1)
        
        try:
            position_list = [int(x.strip()) for x in custom_positions.split(',')]
            log_message(f"Custom positions: {position_list}", verbose)
            
            if len(position_list) != len(genes):
                print(f"ERROR: Number of positions ({len(position_list)}) must match number of genes ({len(genes)})", file=sys.stderr)
                print(f"Genes found: {', '.join(genes)}", file=sys.stderr)
                sys.exit(1)
                
        except ValueError as e:
            print(f"ERROR: Invalid positions format. Use comma-separated integers: {e}", file=sys.stderr)
            sys.exit(1)
    
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
        
        # Store rule information - convert pandas types to native Python types
        rule_info = {
            'position': int(selected_rule['Position']),
            'rule': str(selected_rule['Rule']),
            'score': str(selected_rule.get('Score', 'N/A')),
            'n_correct': int(selected_rule.get('Correct', selected_rule.get('N_Correct', selected_rule.get('N_Correct_Normalized', 0)))) if pd.notna(selected_rule.get('Correct', selected_rule.get('N_Correct', selected_rule.get('N_Correct_Normalized', 0)))) else 'N/A',
            'n_regulators': int(selected_rule.get('N_Regulators', selected_rule.get('N_regulators', 0))) if pd.notna(selected_rule.get('N_Regulators', selected_rule.get('N_regulators', 0))) else 'N/A',
            'mse': str(selected_rule.get('MSE', 'N/A'))
        }
        
        selected_rules[gene] = rule_info
        log_message(f"{gene}: Position {rule_info['position']} - {rule_info['rule']}", verbose)
    
    log_message(f"Selected {len(selected_rules)} rules total", verbose)
    return selected_rules


def extract_genes_from_rule(rule_str, available_genes):
    """
    Extract genes from rule string and determine regulation type
    
    Args:
        rule_str (str): Boolean rule string
        available_genes (set): Set of available gene names
        
    Returns:
        dict: Dictionary mapping gene names to regulation type ('positive' or 'negative')
    """
    genes_in_rule = {}
    
    # Search for all available genes in the rule
    for gene in available_genes:
        # Search for gene as complete word (avoid sub-strings)
        pattern = r'\b' + re.escape(gene) + r'\b'
        matches = list(re.finditer(pattern, rule_str, re.IGNORECASE))
        
        if matches:
            # Determine if the relationship is positive or negative
            # Search for negation indicators near the gene
            for match in matches:
                start_pos = max(0, match.start() - 20)  # Search 20 characters before
                end_pos = min(len(rule_str), match.end() + 20)  # Search 20 characters after
                context = rule_str[start_pos:end_pos].lower()
                
                # Patterns indicating negative regulation
                negative_patterns = [
                    r'not\s+' + re.escape(gene.lower()),
                    r'¬\s*' + re.escape(gene.lower()),
                    r'!\s*' + re.escape(gene.lower()),
                    r'~\s*' + re.escape(gene.lower()),
                    r'-\s*' + re.escape(gene.lower()),
                    r'neg\s+' + re.escape(gene.lower()),
                    r'without\s+' + re.escape(gene.lower()),
                    r'absence\s+of\s+' + re.escape(gene.lower())
                ]
                
                is_negative = any(re.search(pattern, context) for pattern in negative_patterns)
                
                if is_negative:
                    genes_in_rule[gene] = 'negative'
                else:
                    genes_in_rule[gene] = 'positive'
                break  # Only need to evaluate first occurrence
    
    return genes_in_rule


def build_gene_network(selected_rules, verbose):
    """
    Build gene regulatory network from selected rules
    
    Args:
        selected_rules (dict): Dictionary with selected rules per gene
        verbose (bool): Enable verbose output
        
    Returns:
        tuple: NetworkX DiGraph and connection info dictionary
    """
    log_message("Building gene regulatory network...", verbose)
    
    # Create directed graph
    G = nx.DiGraph()
    
    # Get all available genes
    all_genes = set(selected_rules.keys())
    
    # Connection information
    connection_info = defaultdict(list)
    
    # Add all genes as nodes
    for gene in all_genes:
        G.add_node(gene)
    
    # Process each rule to find connections
    for target_gene, rule_info in selected_rules.items():
        rule = rule_info['rule']
        
        log_message(f"Processing {target_gene}: {rule}", verbose)
        
        # Extract genes that appear in the rule with their regulation type
        genes_with_relation = extract_genes_from_rule(rule, all_genes)
        
        if genes_with_relation:
            positive_genes = [gene for gene, reg_type in genes_with_relation.items() if reg_type == 'positive']
            negative_genes = [gene for gene, reg_type in genes_with_relation.items() if reg_type == 'negative']
            autoregulation_genes = [gene for gene in genes_with_relation.keys() if gene == target_gene]
            
            if positive_genes:
                log_message(f"  Positive regulation: {', '.join(sorted(positive_genes))}", verbose)
            if negative_genes:
                log_message(f"  Negative regulation: {', '.join(sorted(negative_genes))}", verbose)
            if autoregulation_genes:
                auto_type = genes_with_relation[target_gene]
                log_message(f"  {auto_type.upper()} autoregulation: {target_gene}", verbose)
            
            # Add edges with regulation type attribute
            for regulator_gene, regulation_type in genes_with_relation.items():
                G.add_edge(regulator_gene, target_gene, regulation_type=regulation_type)
                connection_info[regulator_gene].append({
                    'target': target_gene,
                    'rule': rule,
                    'score': rule_info['score'],
                    'regulation_type': regulation_type,
                    'is_autoregulation': regulator_gene == target_gene
                })
        else:
            log_message(f"  No regulatory genes found", verbose)
    
    log_message(f"Network built with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges", verbose)
    return G, dict(connection_info)


def show_network_statistics(G, connection_info, verbose):
    """
    Display network statistics
    
    Args:
        G (nx.DiGraph): Gene regulatory network
        connection_info (dict): Connection information
        verbose (bool): Enable verbose output
    """
    print(f"\n{'='*80}")
    print("NETWORK STATISTICS")
    print('='*80)
    
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    
    # Statistics by regulation type
    positive_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('regulation_type') == 'positive']
    negative_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('regulation_type') == 'negative']
    autoregulation_edges = [(u, v) for u, v, d in G.edges(data=True) if u == v]
    
    print(f"Positive regulations: {len(positive_edges)}")
    print(f"Negative regulations: {len(negative_edges)}")
    print(f"Autoregulations: {len(autoregulation_edges)}")
    
    # In-degree and out-degree statistics
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())
    
    print(f"\nGenes with most regulators (in-degree):")
    genes_by_in_degree = sorted(in_degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    for gene, degree in genes_by_in_degree:
        if degree > 0:
            # Count positive and negative regulations
            pos_regs = len([e for e in G.in_edges(gene, data=True) if e[2].get('regulation_type') == 'positive'])
            neg_regs = len([e for e in G.in_edges(gene, data=True) if e[2].get('regulation_type') == 'negative'])
            print(f"  {gene}: {degree} regulators total (+{pos_regs}, -{neg_regs})")
    
    print(f"\nGenes that regulate most genes (out-degree):")
    genes_by_out_degree = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    for gene, degree in genes_by_out_degree:
        if degree > 0:
            # Count positive and negative regulations
            pos_regs = len([e for e in G.out_edges(gene, data=True) if e[2].get('regulation_type') == 'positive'])
            neg_regs = len([e for e in G.out_edges(gene, data=True) if e[2].get('regulation_type') == 'negative'])
            print(f"  {gene}: regulates {degree} genes (+{pos_regs}, -{neg_regs})")


def separate_overlapping_nodes(pos, G, min_dist=0.1, max_iterations=50, edge_clearance=0.08):
    """
    Auxiliary function to separate nodes that are too close to each other
    and avoid edges crossing through other nodes
    """
    import numpy as np
    
    pos_array = np.array(list(pos.values()))
    nodes = list(pos.keys())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    
    def point_to_line_distance(point, line_start, line_end):
        """Calculate minimum distance from point to line segment"""
        line_vec = line_end - line_start
        point_vec = point - line_start
        line_len_sq = np.dot(line_vec, line_vec)
        
        if line_len_sq == 0:
            return np.linalg.norm(point - line_start)
        
        t = max(0, min(1, np.dot(point_vec, line_vec) / line_len_sq))
        closest_point = line_start + t * line_vec
        return np.linalg.norm(point - closest_point)
    
    for iteration in range(max_iterations):
        moved = False
        
        # Step 1: Separate nodes that are too close
        for i in range(len(pos_array)):
            for j in range(i + 1, len(pos_array)):
                diff = pos_array[i] - pos_array[j]
                dist = np.linalg.norm(diff)
                
                if dist < min_dist and dist > 0:
                    separation = diff / dist * (min_dist - dist) / 2
                    pos_array[i] += separation
                    pos_array[j] -= separation
                    moved = True
        
        # Step 2: Move nodes that are crossed by edges
        for u, v in G.edges():
            if u == v:  # Skip self-loops
                continue
                
            u_idx = node_to_idx[u]
            v_idx = node_to_idx[v]
            edge_start = pos_array[u_idx]
            edge_end = pos_array[v_idx]
            
            for node in nodes:
                if node == u or node == v:
                    continue
                    
                node_idx = node_to_idx[node]
                node_pos = pos_array[node_idx]
                dist_to_edge = point_to_line_distance(node_pos, edge_start, edge_end)
                
                if dist_to_edge < edge_clearance:
                    # Move node away from edge
                    line_vec = edge_end - edge_start
                    point_vec = node_pos - edge_start
                    line_len_sq = np.dot(line_vec, line_vec)
                    
                    if line_len_sq > 0:
                        t = max(0, min(1, np.dot(point_vec, line_vec) / line_len_sq))
                        closest_point = edge_start + t * line_vec
                        to_node = node_pos - closest_point
                        to_node_len = np.linalg.norm(to_node)
                        
                        if to_node_len > 0:
                            move_distance = edge_clearance - dist_to_edge + 0.02
                            move_vec = (to_node / to_node_len) * move_distance
                            pos_array[node_idx] += move_vec
                            moved = True
        
        if not moved:
            break
    
    # Update position dictionary
    new_pos = {}
    for i, node in enumerate(nodes):
        new_pos[node] = tuple(pos_array[i])
    
    return new_pos


def visualize_network(G, output_base, show_plot, verbose):
    """
    Create network visualization
    
    Args:
        G (nx.DiGraph): Gene regulatory network
        output_base (str): Base name for output files
        show_plot (bool): Whether to display the plot
        verbose (bool): Enable verbose output
        
    Returns:
        dict: Node positions for use in XGMML export
    """
    log_message("Creating network visualization...", verbose)
    
    # Adjust figure size based on number of nodes
    num_nodes = len(G.nodes())
    if num_nodes <= 20:
        figsize = (16, 14)
        node_size = 3000
        font_size = 12
    elif num_nodes <= 50:
        figsize = (20, 18)
        node_size = 2000
        font_size = 10
    else:
        figsize = (24, 20)
        node_size = 1500
        font_size = 8
    
    plt.figure(figsize=figsize)
    
    # Choose layout based on network size
    if num_nodes <= 15:
        pos = nx.kamada_kawai_layout(G)
    elif num_nodes <= 30:
        pos = nx.spring_layout(G, k=5/np.sqrt(num_nodes), iterations=100, seed=42)
        pos = separate_overlapping_nodes(pos, G, min_dist=0.12, edge_clearance=0.08)
    else:
        pos_spring = nx.spring_layout(G, k=8/np.sqrt(num_nodes), iterations=50, seed=42)
        pos = nx.fruchterman_reingold_layout(G, pos=pos_spring, k=6/np.sqrt(num_nodes), iterations=50)
        pos = separate_overlapping_nodes(pos, G, min_dist=0.15, edge_clearance=0.1)

    # Calculate center for autoregulation orientation
    if pos:
        center_x = sum(x for x, y in pos.values()) / len(pos)
        center_y = sum(y for x, y in pos.values()) / len(pos)
        center = (center_x, center_y)
    else:
        center = (0, 0)

    def calculate_outward_angle(node_pos, center):
        """Calculate angle from center to node (for outward orientation)"""
        dx = node_pos[0] - center[0]
        dy = node_pos[1] - center[1]
        return np.arctan2(dy, dx)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          node_color='lightblue', 
                          node_size=node_size,
                          alpha=0.8,
                          edgecolors='black',
                          linewidths=1.5)
    
    # Separate edges by regulation type and autoregulation
    positive_edges = [(u, v) for u, v, d in G.edges(data=True) 
                      if d.get('regulation_type') == 'positive' and u != v]
    negative_edges = [(u, v) for u, v, d in G.edges(data=True) 
                      if d.get('regulation_type') == 'negative' and u != v]
    positive_autoregulation = [(u, v) for u, v, d in G.edges(data=True) 
                               if u == v and d.get('regulation_type') == 'positive']
    negative_autoregulation = [(u, v) for u, v, d in G.edges(data=True) 
                               if u == v and d.get('regulation_type') == 'negative']

    # Pre-calculate curvatures to avoid overlap
    all_edges = set()
    for u, v in positive_edges + negative_edges:
        all_edges.add((u, v))

    # Identify mutual regulation pairs
    mutual_pairs = set()
    for u, v in all_edges:
        if (v, u) in all_edges:
            normalized_pair = tuple(sorted([u, v]))
            mutual_pairs.add(normalized_pair)

    # Assign curvatures
    curvatures = {}
    curvature_base = 0.3 if num_nodes > 30 else 0.2
    
    for pair in mutual_pairs:
        node1, node2 = pair
        if (node1, node2) in all_edges:
            curvatures[(node1, node2)] = curvature_base
        if (node2, node1) in all_edges:
            curvatures[(node2, node1)] = curvature_base

    for u, v in all_edges:
        if (u, v) not in curvatures:
            curvatures[(u, v)] = 0.0
    
    def get_curvature(u, v):
        return curvatures.get((u, v), 0.0)
    
    # Adjust line width based on network size
    line_width = max(1, 3 - num_nodes/20)
    
    # Draw positive regulation edges (blue arrows)
    if positive_edges:
        for u, v in positive_edges:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            
            dx = x2 - x1
            dy = y2 - y1
            length = (dx**2 + dy**2)**0.5
            
            if length > 0:
                dx_norm = dx / length
                dy_norm = dy / length

                node_radius = np.sqrt(node_size/3000) * 0.08
                start_offset = node_radius
                end_offset = node_radius
                
                x_start = x1 + dx_norm * start_offset
                y_start = y1 + dy_norm * start_offset
                x_end = x2 - dx_norm * end_offset
                y_end = y2 - dy_norm * end_offset

                rad = get_curvature(u, v)

                arrow = FancyArrowPatch((x_start, y_start), (x_end, y_end),
                                    connectionstyle=f"arc3,rad={rad}",
                                    color='#3690C0',
                                    lw=line_width, alpha=0.7,
                                    arrowstyle='-|>', mutation_scale=20 if num_nodes > 30 else 30)
                plt.gca().add_patch(arrow)
    
    # Draw negative regulation edges (red lines with T)
    if negative_edges:
        for u, v in negative_edges:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            
            dx = x2 - x1
            dy = y2 - y1
            length = (dx**2 + dy**2)**0.5
            
            if length > 0:
                dx_norm = dx / length
                dy_norm = dy / length
                
                node_radius = np.sqrt(node_size/3000) * 0.08
                start_offset = node_radius
                end_offset = node_radius
                
                x_start = x1 + dx_norm * start_offset
                y_start = y1 + dy_norm * start_offset
                x_end = x2 - dx_norm * end_offset
                y_end = y2 - dy_norm * end_offset

                rad = get_curvature(u, v)

                path = FancyArrowPatch((x_start, y_start), (x_end, y_end),
                                   connectionstyle=f"arc3,rad={rad}",
                                   color='#D7301F', lw=line_width, alpha=0.8,
                                   arrowstyle='|-|, widthA=0, widthB=0.4',
                                   mutation_scale=20 if num_nodes > 30 else 30)
                plt.gca().add_patch(path)
    
    # Draw positive autoregulation (blue self-loops)
    if positive_autoregulation:
        for u, v in positive_autoregulation:
            x, y = pos[u]
            outward_angle = calculate_outward_angle((x, y), center)
            r_node = np.sqrt(node_size/5000) * 0.08

            offset_angle = np.pi / 12
            angle1 = outward_angle - offset_angle
            angle2 = outward_angle + offset_angle
            start = (x + r_node * np.cos(angle1), y + r_node * np.sin(angle1))
            end = (x + r_node * np.cos(angle2), y + r_node * np.sin(angle2))
            
            loop = FancyArrowPatch(start, end,
                                connectionstyle=f"arc3,rad=3",
                                color='#3690C0', lw=line_width, linestyle='solid',
                                arrowstyle='->', mutation_scale=20 if num_nodes > 30 else 30, alpha=0.8)
            plt.gca().add_patch(loop)
    
    # Draw negative autoregulation (red self-loops)
    if negative_autoregulation:
        for u, v in negative_autoregulation:
            x, y = pos[u]
            outward_angle = calculate_outward_angle((x, y), center)
            r_node = np.sqrt(node_size/3000) * 0.08
            
            offset_angle = np.pi / 12
            angle1 = outward_angle - offset_angle
            angle2 = outward_angle + offset_angle

            start = (x + r_node * np.cos(angle1), y + r_node * np.sin(angle1))
            end = (x + r_node * np.cos(angle2), y + r_node * np.sin(angle2))

            path = FancyArrowPatch(
                start, end,
                connectionstyle=f"arc3,rad=3",
                color='#D7301F', lw=line_width, linestyle='solid',
                arrowstyle='|-|, widthA=0, widthB=0.4',
                mutation_scale=20 if num_nodes > 30 else 30, alpha=0.8
            )
            plt.gca().add_patch(path)
            
    # Draw gene labels
    nx.draw_networkx_labels(G, pos, 
                           font_size=font_size,
                           font_weight='bold',
                           font_color='black')
    
    # Create legend
    from matplotlib.lines import Line2D
    legend_elements = []
    
    if positive_edges or positive_autoregulation:
        total_positive = len(positive_edges) + len(positive_autoregulation)
        legend_elements.append(Line2D([0], [0], color='#3690C0', lw=3,
                                      label=f'Positive Regulation ({total_positive})'))
    if negative_edges or negative_autoregulation:
        total_negative = len(negative_edges) + len(negative_autoregulation)
        legend_elements.append(Line2D([0], [0], color='#D7301F', lw=3, 
                                    label=f'Negative Regulation ({total_negative})'))
    
    if legend_elements:
        plt.legend(handles=legend_elements, loc='upper right', fontsize=max(9, font_size-1))
    
    plt.title("Boolean Gene Regulatory Network", 
              size=max(14, 18-num_nodes//10), weight='bold')
    plt.axis('off')
    plt.tight_layout()
    
    # Save files
    dpi = 400 if num_nodes > 30 else 300
    plt.savefig(f"{output_base}.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(f"{output_base}.svg", bbox_inches='tight')
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    log_message(f"Network visualization saved as {output_base}.png and {output_base}.svg", verbose)
    
    # Return positions for XGMML export
    return pos


def create_xgmml_file(G, selected_rules, connection_info, pos, filename, verbose):
    """
    Create XGMML file with embedded visual properties for Cytoscape
    
    Args:
        G (nx.DiGraph): Gene regulatory network
        selected_rules (dict): Selected rules
        connection_info (dict): Connection information
        pos (dict): Node positions from NetworkX layout
        filename (str): Output filename
        verbose (bool): Enable verbose output
    """
    log_message("Creating XGMML file with visual properties and layout positions...", verbose)
    
    # Create root element
    root = ET.Element("graph", {
        'xmlns:dc': "http://purl.org/dc/elements/1.1/",
        'xmlns:xlink': "http://www.w3.org/1999/xlink", 
        'xmlns:rdf': "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        'xmlns:cy': "http://www.cytoscape.org",
        'xmlns': "http://www.cs.rpi.edu/XGMML",
        'directed': "1",
        'label': "Boolean Gene Regulatory Network"
    })
    
    # Add graph attributes
    ET.SubElement(root, "att", {'name': 'documentVersion', 'value': '1.1'})
    ET.SubElement(root, "att", {'name': 'networkMetadata'})
    ET.SubElement(root, "att", {'name': 'backgroundColor', 'value': '#FFFFFF'})
    
    # Add nodes
    for node_id in G.nodes():
        node_elem = ET.SubElement(root, "node", {
            'id': str(node_id),
            'label': str(node_id)
        })
        
        # Add node attributes
        ET.SubElement(node_elem, "att", {'name': 'name', 'value': str(node_id)})
        ET.SubElement(node_elem, "att", {'name': 'node_type', 'value': 'gene'})
        
        # Add rule information if available
        if node_id in selected_rules:
            rule_info = selected_rules[node_id]
            ET.SubElement(node_elem, "att", {'name': 'position', 'value': str(rule_info['position'])})
            ET.SubElement(node_elem, "att", {'name': 'rule', 'value': str(rule_info['rule'])})
            ET.SubElement(node_elem, "att", {'name': 'score', 'value': str(rule_info['score'])})
            ET.SubElement(node_elem, "att", {'name': 'n_correct', 'value': str(rule_info['n_correct'])})
            ET.SubElement(node_elem, "att", {'name': 'n_regulators', 'value': str(rule_info['n_regulators'])})
            ET.SubElement(node_elem, "att", {'name': 'mse', 'value': str(rule_info['mse'])})
        
        # Calculate scaled coordinates for Cytoscape
        # NetworkX coordinates are typically in range [-1, 1], scale to appropriate size for Cytoscape
        if node_id in pos:
            x_coord = pos[node_id][0] * 300  # Scale factor for better visualization
            y_coord = pos[node_id][1] * 300  # Scale factor for better visualization
        else:
            # Fallback coordinates if node not in position dict
            x_coord = 100.0
            y_coord = 100.0
        
        # Add graphics element with actual coordinates and visual properties
        graphics = ET.SubElement(node_elem, "graphics", {
            'type': 'ELLIPSE',
            'x': f'{x_coord:.1f}',
            'y': f'{y_coord:.1f}',
            'w': '75.0',
            'h': '50.0',
            'fill': '#7BCCC4',
            'outline': '#000000'
        })
        
        # Add border width and font properties inside graphics element
        ET.SubElement(graphics, "att", {'name': 'node_border_width', 'value': '2.0'})
        ET.SubElement(graphics, "att", {'name': 'NODE_LABEL_FONT_FACE', 'value': 'Arial Black'})
        ET.SubElement(graphics, "att", {'name': 'NODE_LABEL_FONT_SIZE', 'value': '12'})
    
    # Add edges
    edge_id = 0
    for source, target, edge_data in G.edges(data=True):
        regulation_type = edge_data.get('regulation_type', 'unknown')
        is_autoregulation = (source == target)
        
        edge_elem = ET.SubElement(root, "edge", {
            'id': str(edge_id),
            'source': str(source),
            'target': str(target),
            'label': f"{source} -> {target}"
        })
        edge_id += 1
        
        # Add edge attributes  
        ET.SubElement(edge_elem, "att", {'name': 'regulation_type', 'value': regulation_type})
        ET.SubElement(edge_elem, "att", {'name': 'is_autoregulation', 'value': str(is_autoregulation)})
        ET.SubElement(edge_elem, "att", {'name': 'interaction_type', 'value': 'activation' if regulation_type == 'positive' else 'inhibition' if regulation_type == 'negative' else 'unknown'})
        
        # Add connection-specific information
        if source in connection_info:
            for conn in connection_info[source]:
                if conn['target'] == target and conn['regulation_type'] == regulation_type:
                    ET.SubElement(edge_elem, "att", {'name': 'rule', 'value': str(conn['rule'])})
                    ET.SubElement(edge_elem, "att", {'name': 'score', 'value': str(conn['score'])})
                    break
        
        # Add visual properties as graphics element
        graphics_attrs = {
            'width': '2.0',
            'style': 'SOLID'
        }
        
        # Set colors and arrow shapes based on regulation type
        if regulation_type == 'positive':
            graphics_attrs.update({
                'fill': '#3690C0',           # Stroke color
                'targetArrow': 'DELTA',      # Target arrow shape
                'targetArrowColor': '#3690C0' # Target arrow color
            })
        elif regulation_type == 'negative':
            graphics_attrs.update({
                'fill': '#D7301F',           # Stroke color
                'targetArrow': 'T',          # Target arrow shape
                'targetArrowColor': '#D7301F' # Target arrow color
            })
        else:
            graphics_attrs.update({
                'fill': '#666666',
                'targetArrow': 'ARROW',
                'targetArrowColor': '#666666'
            })
            
        graphics = ET.SubElement(edge_elem, "graphics", graphics_attrs)
    
    # Write to file
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ", level=0)  # Pretty print
    tree.write(filename, encoding='utf-8', xml_declaration=True)
    
    log_message(f"XGMML file created successfully: {filename}", verbose)


def save_network_formats(G, selected_rules, connection_info, pos, output_base, verbose):
    """
    Save network in multiple formats
    
    Args:
        G (nx.DiGraph): Gene regulatory network
        selected_rules (dict): Selected rules
        connection_info (dict): Connection information
        pos (dict): Node positions from NetworkX layout
        output_base (str): Base name for output files
        verbose (bool): Enable verbose output
    """
    log_message("Saving network in multiple formats...", verbose)
    
    # Save as XGMML format (Cytoscape native with visual properties and layout)
    xgmml_file = f"{output_base}.xgmml"
    create_xgmml_file(G, selected_rules, connection_info, pos, xgmml_file, verbose)
    
    # Save selected rules (JSON)
    rules_file = f"{output_base}_selected_rules.json"
    with open(rules_file, 'w', encoding='utf-8') as f:
        json.dump(selected_rules, f, indent=2, ensure_ascii=False)
    log_message(f"Selected rules saved: {rules_file}", verbose)


def visualize_boolean_network(args):
    """
    Main function to visualize Boolean gene regulatory network
    
    Args:
        args: Parsed command line arguments
    """
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Read rules table
    df = read_rules_table(args.input_file, args.verbose)
    if df is None:
        sys.exit(1)
    
    # Select rules based on criteria
    selected_rules = select_rules_by_criteria(df, args.criteria, args.positions, args.verbose)
    
    if not selected_rules:
        print("ERROR: No rules selected", file=sys.stderr)
        sys.exit(1)
    
    # Build gene network
    G, connection_info = build_gene_network(selected_rules, args.verbose)
    
    # Show network statistics
    show_network_statistics(G, connection_info, args.verbose)
    
    # Create output paths
    input_file_path = Path(args.input_file)
    input_name = input_file_path.stem
    input_dir = input_file_path.parent.absolute()
    
    if args.output:
        # User specified custom output path
        output_base = args.output
        network_dir = os.path.dirname(output_base)
    elif args.output_dir:
        # User specified custom output directory
        network_dir = args.output_dir
        output_base = os.path.join(network_dir, f"{input_name}_network")
    else:
        # Default: create 'network' subdirectory in the same folder as input file
        network_dir = os.path.join(input_dir, "network")
        output_base = os.path.join(network_dir, f"{input_name}_network")
    
    # Create output directory
    log_message(f"Creating output directory: {network_dir}", args.verbose)
    os.makedirs(network_dir, exist_ok=True)
    
    # Visualize network and get positions
    pos = visualize_network(G, output_base, args.show, args.verbose)
    
    # Save network in multiple formats using the calculated positions
    save_network_formats(G, selected_rules, connection_info, pos, output_base, args.verbose)
    
    print(f"\n{'='*80}")
    print("PROCESS COMPLETED SUCCESSFULLY")
    print('='*80)
    print(f"Files generated in: {network_dir}")
    print(f"- {input_name}_network.png (network visualization)")
    print(f"- {input_name}_network.svg (network visualization)")
    print(f"- {input_name}_network.xgmml (Cytoscape native format with embedded visual properties)")
    print(f"- {input_name}_network_selected_rules.json (selected rules)")
    print('='*80)
    print("\nTo use in Cytoscape:")
    print("1. Import the .xgmml file (File > Import > Network)")
    print("2. Visual properties should be applied automatically!")
    print("   - Positive regulations: #3690C0 color with DELTA arrows")
    print("   - Negative regulations: #D7301F color with T arrows")
    print("3. Layout positions are preserved from NetworkX calculation")
    print('='*80)


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Visualize Boolean gene regulatory networks from inferred rules',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 visualize_network.py -i rules_by_gene.tsv
  python3 visualize_network.py -i rules.tsv -c top_position -o my_network -v
  python3 visualize_network.py -i results/rules.tsv -d output/ --show
  python3 visualize_network.py -i rules.tsv -c custom_positions -p "1,1,2,3,1" -v

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

Output Formats:
  - PNG and SVG images for visualization
  - XGMML format for import into Cytoscape with embedded visual properties and layout
  - JSON file with selected rules

Notes:
  - Input file should contain columns: Gene, Position, Rule
  - Optional columns: Score, MSE, N_Correct_Normalized, N_Regulators
  - Network automatically detects positive/negative regulation from rule syntax
  - Self-loops (autoregulation) are supported and visualized differently
  - Layout algorithm adapts to network size for optimal visualization
  - Node positions are preserved between matplotlib and Cytoscape for consistency
  - If specified position not found for a gene, uses minimum available position
        """
    )
    
    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input_file', type=str, required=True,
                         help='Input TSV file containing rules by gene')
    
    # Optional parameters
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-c', '--criteria', type=str, default='top_position',
                         choices=['top_position', 'custom_positions'],
                         help='Rule selection criteria (default: top_position)')
    optional.add_argument('-p', '--positions', type=str, default=None,
                         help='Comma-separated list of positions for custom_positions criteria (e.g., "1,1,2,3,1")')
    optional.add_argument('-o', '--output', type=str, default=None,
                         help='Output base name for files (default: auto-generated from input)')
    optional.add_argument('-d', '--output_dir', type=str, default=None,
                         help='Output directory for generated files (default: network/ subdirectory in input file location)')
    optional.add_argument('--show', action='store_true',
                         help='Display the network plot on screen')
    optional.add_argument('--no-show', dest='show', action='store_false',
                         help='Do not display the network plot (default)')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')
    
    parser.add_argument('--version', action='version', version='Boolean Network Visualizer v1.0')
    parser.set_defaults(show=False)
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Run visualization
    try:
        visualize_boolean_network(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()