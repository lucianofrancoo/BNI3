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

Simulates the dynamic trajectory of the Boolean network from a starting gene-activity configuration until it converges into a known attractor. Additionally, **every sample in the binarized expression matrix is mapped to its corresponding attractor basin**, revealing whether each experimental condition already sits inside an attractor or is transitioning toward one.

#### Required inputs

| Flag | Description |
|------|-------------|
| `-a` | `attractors.tsv` — output of `1.BNI3_Attractors.py` |
| `-r` | `selected_rules.tsv` — output of `1.BNI3_Attractors.py` (or full rules table) |
| `-b` | Binarized expression matrix TSV — output of `1.Binarization/BNI3_SSD.py` or `BNI3_WCSS.py` (genes as columns, samples as rows) |

#### Optional inputs

| Flag | Description | Default |
|------|-------------|---------|
| `-s` | Initial state as a binary string (`"1010110010110"`) or comma-separated active gene names (`"ABF3,ABF4"`) | **First row of the binarized matrix** |
| `-o` | Output directory | Same directory as `attractors.tsv` |
| `-ob` | Base name for output files | `trajectory` |
| `--max-steps` | Maximum simulation steps per trajectory | `1000` |
| `-v` | Verbose output | off |

**Basic usage — use first matrix sample as starting point (default):**
```bash
python3 2.BNI3_Path_to_Attractors.py \
    -a Example/attractors.tsv \
    -r Example/selected_rules.tsv \
    -b ../1.Binarization/Example/Counts_lite_binarized_SSD.tsv \
    -o Example/
```

**Override the initial state with an explicit binary string:**
```bash
python3 2.BNI3_Path_to_Attractors.py \
    -a Example/attractors.tsv \
    -r Example/selected_rules.tsv \
    -b ../1.Binarization/Example/Counts_lite_binarized_SSD.tsv \
    -s "1010110010110" -o Example/
```

**Override with a comma-separated list of active genes:**
```bash
python3 2.BNI3_Path_to_Attractors.py \
    -a Example/attractors.tsv \
    -r Example/selected_rules.tsv \
    -b ../1.Binarization/Example/Counts_lite_binarized_SSD.tsv \
    -s "ABF3,ABF4,DREB2A" -v
```

#### Output files

| File | Description |
|------|-------------|
| `<base>.tsv` | Complete step-by-step trajectory from the initial state to the attractor |
| `<base>_trajectory.png` / `.svg` | Heatmap visualization of gene activity along the trajectory |
| `<base>_matrix_attractor_mapping.tsv` | Attractor assignment for **every sample** in the binarized matrix: whether the sample is already inside an attractor (`in_attractor`), transitioning toward one (`transient_to_attractor`), or did not converge within the step limit (`did_not_converge`) |

The mapping table includes: `sample`, `binary_state`, `status`, `attractor_id`, `attractor_type`, `step_in_cycle`, `steps_to_attractor`, `basin_size`, `basin_percentage`.

### 3. Attractors Visualizer (`3.BNI3_Visualize_Attractors.py`)
Creates rich graphical diagram representations and heatmaps to visually interpret the detected attractors and state transitions from `attractors.tsv`.

**Basic Usage:**
```bash
python3 3.BNI3_Visualize_Attractors.py -i Example/attractors.tsv --heatmap --network
```
*Outputs generated:* High quality image formats (`.png` / `.svg`) representing the boolean transition networks and expression state heatmaps.
