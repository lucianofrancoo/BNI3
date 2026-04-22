# BNI3: Interactive Usage Vignette

Welcome to the **Boolean Network Inference 3 (BNI3)** official vignette. This guide will walk you through a complete, end-to-end analysis of gene expression data using our interactive bash launcher.

BNI3 is designed to abstract away the complexity of mathematical modeling, providing you with an automated pipeline that takes raw continuous data, binarizes it, infers logical rules using Gene Expression Programming (GEP), and plots the topological attractors representing the stable states of your biological system.

---

## 1. Getting Started

To ensure reproducibility across different operating systems, we highly recommend using `conda` to resolve all the Python dependencies required for BNI3.

```bash
# 1. Clone the repository
git clone https://github.com/lucianofrancoo/BNI3.git
cd BNI3

# 2. Create and activate the Conda environment
conda env create -f environment.yml
conda activate bni3_env

# 3. Launch the Interactive Pipeline
./bni3_launcher.sh
```

Upon executing the launcher, you will be greeted by the Main Menu:

```text
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║                    🧬 BNI3 PIPELINE LAUNCHER                  ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝

Please select the pipeline step you want to execute:

  1) Binarization Module
  2) Rules Inference Module
  3) Attractors Analysis Module
  4) Exit
```

---

## 2. Module 1: Binarization

Continuous gene expression data (e.g., RNA-Seq) needs to be discretized. BNI3 natively supports **SSD** (StepMiner) and **WCSS** (K-Means).

### Input Formatting
Your initial continuous matrix should have **Genes as rows** and **Timepoints/Samples as columns**:

| ID    | T1  | T2  | T3  |
|-------|-----|-----|-----|
| GeneA | 0.5 | 2.1 | 1.8 |
| GeneB | 1.2 | 0.3 | 0.4 |

### Interactive Execution
Select `1` in the Main Menu and follow the prompts. You simply drag and drop the path to your raw `.tsv`. The pipeline will automatically generate the corresponding **Binarized Matrix** and a **Behavior Review Summary** indicating how many unique patterns were found.

---

## 3. Module 2: Rules Inference & Evaluation

This is the core of BNI3. Using Gene Expression Programming, the module searches millions of possible Boolean logic combinations (`AND`, `OR`, `NOT`) to find the mathematical equations that best describe how your genes transition from one timepoint to the next.

### Input Formatting
The inference engine expects the **Binarized Matrix** to be transposed compared to the original input (Genes as Columns, Timepoints as Rows):

| GeneA | GeneB |
|-------|-------|
| 0     | 1     |
| 1     | 0     |

### Interactive Execution
Select `2` in the Main Menu. Provide the path to your *Continuous Matrix* and your newly generated *Binarized Matrix*.

> **Note on Multithreading**: BNI3 automatically detects your CPU hardware through `$(nproc)` and will distribute the evolutionary workload across 100% of your available cores, drastically reducing computational time.

The orchestrator will output the best rules per gene:
`rules_by_gene_evaluated.tsv`

---

## 4. Module 3: Attractors & Visualization

Once the rules are established, we can calculate the **Attractors**. Attractors are the steady states (fixed points or cycles) where the network eventually settles, often representing biological phenotypes like cell death, proliferation, or steady homeostasis.

Select `3` in the Main Menu and choose `Find Attractors` and then `Visualize Attractors Landscape`.

The pipeline will use the `rules_by_gene_evaluated.tsv` generated in the previous step and will render the topological network mapping the basins of attraction.

### Output
BNI3 will automatically generate high-resolution Vector (`.svg`) and Pixel (`.png`) images of your network's topology.
