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
├── 3.Attractors/            # Stage 3: Attractor mapping and visualization
├── bni3_launcher.sh         # Interactive CLI wrapper for the complete pipeline
└── README.md                # Project documentation
```

## Getting Started

### Interactive Wrapper (Recommended for standard usage)
We provide an interactive bash launcher designed to simplify the execution of the entire pipeline using robust default settings. 

To start the interactive CLI:
```bash
./bni3_launcher.sh
```
> **Note**: The interactive launcher is a simplified wrapper. If you require highly detailed, granular control over every possible algorithm hyperparameter, please invoke the underlying Python scripts inside each folder directly.

### Manual & Advanced Pipeline
If you need full parameter control, navigate to each module's directory and run the scripts directly from your terminal:
1. **Binarize your data**: Navigate to the `1.Binarization/` directory and run your continuous data through the `BNI3_WCSS.py` or `BNI3_SSD.py` scripts.
2. **Infer Rules**: Feed the binary output to the `2.Rules_Inference/` module and run `1.BNI3_Boolean_Rules_Inference.py` to generate the logical rules.
3. **Analyze Post-Inference**: Use the scripts inside `3.Attractors/` to discover and map the steady states of your inferred network.

## Contributing
Contributions, issues, and feature requests are welcome!

## License
[Specify License Here]
