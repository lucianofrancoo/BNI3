# BNI3: Boolean Network Inference and Analysis

Welcome to the **BNI3** project! BNI3 is a comprehensive pipeline designed for the binarization of gene expression data, inference of Boolean network rules, and the analysis of their dynamic attractors. This toolkit is well-suited for systems biology research, gene regulatory network (GRN) modeling, and dynamic state analysis.

## Overview

The BNI3 pipeline is structured into three main stages, each corresponding to a specific module in the project:

### 1. Binarization (`1.Binarization/`)
Continuous gene expression data (e.g., RNA-Seq counts) must be discretized into binary states (0 for inactive/unexpressed, 1 for active/expressed) to model them as Boolean networks. 
* **Key Scripts**: 
  * `BNI3_SSD.py`, `BNI3_WCSS.py`: Different algorithms for optimal binarization (e.g., Within-Cluster Sum of Squares).
* **Input/Output Data**: Example continuous data and outputs are provided in `Example/`.

### 2. Boolean Rules Inference (`2.Rules_Inference/`)
Once the data is binarized, this module infers the logical relationships (Boolean rules) between different genes/nodes over time.
* **Key Scripts**:
  * `1.BNI3_Boolean_Rules_Inference.py`: The core algorithm to infer Boolean transition rules from the binarized datasets.
  * `2.BNI3_Analyze_results.py`: Analyzes the output of the inference process.
  * `3.BNI3_Evaluate_rules.py`: Evaluates the accuracy and robustness of the inferred rules.
  * `4.BNI3_Boolean_network_visualizer.py`: Generates visual representations of the inferred network topologies.

### 3. Attractors Analysis (`3.Attractors/`)
Boolean networks eventually settle into steady states or cycles known as *attractors*. These represent biological phenotypes or stable cell states.
* **Key Scripts**:
  * `1.BNI3_Attractors.py`: Calculates the attractors and their corresponding basins of attraction.
  * `2.BNI3_Path_to_Attractors.py`: Simulates the sequence of state transitions (trajectories) from an initial point towards an attractor.
  * `3.BNI3_Visualize_Attractors.py`: Visualizes the attractors (e.g., as networks or heatmaps).

## Project Structure

```text
BNI3/
├── 1.Binarization/          # Stage 1: Discretization of continuous data
├── 2.Rules_Inference/       # Stage 2: Boolean network rule inference
├── Atracttors/              # Stage 3: Attractor mapping and visualization
├── BoolNet_comparison/      # Benchmarks against the BoolNet R package
├── DREAM4/                  # Evaluation using standard DREAM4 in silico network challenges
└── ...                      # Additional datasets and specific network tests (guardCells, etc.)
```

## Getting Started

### Prerequisites
Make sure you have the following installed:
* **Python 3.x**: Required for the core pipeline (inference, python-based binarization, and attractors).
* **R**: Required for running the comparison scripts (BoolNet) and some binarization/validation scripts (`.r` files).

### Usage Pipeline
1. **Binarize your data**: Navigate to the `1.Binarization/` directory and run your continuous data through the WCSS or SSD Python/R scripts.
2. **Infer Rules**: Move the binary output to the `2.Rules_Inference/` module and run `BNI3_Boolean_Rules_Inference.py` to generate the logical rules.
3. **Analyze Post-Inference**: Use the scripts inside `Atracttors/` to find the stable steady states of your inferred network.

## Benchmarks & Comparisons
The pipeline includes modules (`BoolNet_comparison/` and `DREAM4/`) that benchmark the BNI3 inference methodology against established tools like the BoolNet R package and the DREAM4 *in silico* network challenges.

## Contributing
Contributions, issues, and feature requests are welcome!

## License
[Specify License Here]
