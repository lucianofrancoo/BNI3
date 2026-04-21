# 3. Attractors Analysis (`3.Attractors/`)

Once the logical relationships between genes are established and evaluated as a Boolean network, this module computes its eventual steady states or cycles, known as **attractors**. These attractors represent the biological phenotypes, specific distinct cell states, or functional modes the network can settle into.

The module provides tools for exhaustive attractor detection, trajectory (path) simulations, and rich aesthetic visualizations of the network's dynamics.

## Pipeline Scripts and Usage

### 1. Attractors Finder (`1.BNI3_Attractors.py`)
This script analyzes your Boolean network rules to detect all reachable attractors (steady states or limit cycles) and precisely calculates the sizes of their basins of attraction. 

**Basic Usage:**
```bash
python3 1.BNI3_Attractors.py -i ../2.Rules_Inference/Example/rules_by_gene_evaluated.tsv -o Example/
```
*Outputs generated:* `attractors.tsv` and `selected_rules.tsv` (useful for tracking state sizes and exact rules used).

**Perturbation Analysis (Simulating Overexpression or Knockout):**
You can also evaluate the effect of mutations or forced states on the network's attractors landscape by using the `-m` flag. For example, to simulate an overexpression of ABF4 and a knockout of MYB44:
```bash
python3 1.BNI3_Attractors.py -i ../2.Rules_Inference/Example/rules_by_gene_evaluated.tsv -o Example/ -m "ABF4:1,MYB44:0"
```

### 2. Path to Attractors Simulator (`2.BNI3_Path_to_Attractors.py`)
Simulates and tracks the dynamic trajectory in the network starting from a specific initial configuration of active/inactive genes until it converges into a known attractor.

**Basic Usage:**
```bash
python3 2.BNI3_Path_to_Attractors.py -a Example/attractors.tsv -r Example/selected_rules.tsv -s "1010110010110" -o Example/
```
*(You can pass the initial state `-s` as a complete binary string or as a comma-separated list of active genes: `"GENE1,GENE2"`)*

### 3. Attractors Visualizer (`3.BNI3_Visualize_Attractors.py`)
Creates rich graphical diagram representations and heatmaps to visually interpret the detected attractors and state transitions from `attractors.tsv`.

**Basic Usage:**
```bash
python3 3.BNI3_Visualize_Attractors.py -i Example/attractors.tsv --heatmap --network
```
*Outputs generated:* High quality image formats (`.png` / `.svg`) representing the boolean transition networks and expression state heatmaps.
