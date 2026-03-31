# BNI3: Boolean Network Inference and Analysis

Welcome to the **BNI3** project! BNI3 is a comprehensive pipeline designed for the binarization of gene expression data, inference of Boolean network rules, and the analysis of their dynamic attractors. This toolkit is well-suited for systems biology research, gene regulatory network (GRN) modeling, and dynamic state analysis.

## Overview

The BNI3 pipeline is structured into three main stages, each corresponding to a specific module in the project:

### 1. Binarization (`Binarization/`)
Continuous gene expression data (e.g., RNA-Seq counts) must be discretized into binary states (0 for inactive/unexpressed, 1 for active/expressed) to model them as Boolean networks. 
* **Key Scripts**: 
  * `BNI3_SSD.py`, `BNI3_WCSS.py`: Different algorithms for optimal binarization (e.g., Within-Cluster Sum of Squares).
  * `reducir_matriz.r`, `WCSS.r`: R scripts for matrix reduction and statistical clustering validation.
* **Input**: Continuous expression matrices (`continuous_expression_matrix/`).
* **Output**: Binarized expression matrices (`binary_expression_matrix/`).

### 2. Boolean Rules Inference (`Rules_Inference/`)
Once the data is binarized, this module infers the logical relationships (Boolean rules) between different genes/nodes over time.
* **Key Scripts**:
  * `BNI3_Boolean_Rules_Inference.py`: The core algorithm to infer Boolean transition rules from the binarized datasets.
  * `BNI3_Evaluate_rules.py`: Evaluates the accuracy and robustness of the inferred rules.
  * `BNI3_Analyze_results.py`: Analyzes the output of the inference process.
  * `BNI3_Boolean_network_visualizer.py`: Generates visual representations of the inferred network topologies.

### 3. Attractors Analysis (`Atracttors/`)
Boolean networks eventually settle into steady states or cycles known as *attractors*. These represent biological phenotypes or stable cell states.
* **Key Scripts**:
  * `BNI3_Attractors.py`: Computes and identifies the attractors from the inferred Boolean rules.
  * `BNI3_Path_to_Attractors.py`: Analyzes the state transition graphs and paths leading to specific attractors.
  * `BNI3_Visualize_Attractors.py`: Visualizes the attractor states and the transition landscape.

## Project Structure

```text
BNI3/
├── Binarization/            # Stage 1: Discretization of continuous data
├── Rules_Inference/         # Stage 2: Boolean network rule inference
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
1. **Binarize your data**: Navigate to the `Binarization/` directory and run your continuous data through the WCSS or SSD Python/R scripts.
2. **Infer Rules**: Move the binary output to the `Rules_Inference/` module and run `BNI3_Boolean_Rules_Inference.py` to generate the logical rules.
3. **Analyze Post-Inference**: Use the scripts inside `Atracttors/` to find the stable steady states of your inferred network.

## Benchmarks & Comparisons
The pipeline includes modules (`BoolNet_comparison/` and `DREAM4/`) that benchmark the BNI3 inference methodology against established tools like the BoolNet R package and the DREAM4 *in silico* network challenges.

## Contributing
Contributions, issues, and feature requests are welcome!

## License
[Specify License Here]
