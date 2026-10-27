# Binarization Module

This directory contains the necessary scripts and tools to discretize continuous gene expression data into binary states (0 or 1). Discretizing the data is the fundamental first step required to model gene interactions as a Boolean network.

We provide two primary algorithms for binarization: **SSD** (Short Series Discretization) and **WCSS** (Within-Cluster Sum of Squares).

## 0. Binarizator (recommended entry point)

`BNI3_Binarizator.py` runs the whole stage from one command: it binarizes a counts matrix with **both** SSD and WCSS, then hands the results to the behavior reviewer. The individual scripts described below still work standalone — this one only saves you from launching all three by hand.

```bash
python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv
```

For an input named `Counts_lite.tsv` this produces:

| File | Content |
|---|---|
| `Counts_lite_binarized_SSD.tsv` | Binarized matrix, SSD |
| `Counts_lite_binarized_WCSS.tsv` | Binarized matrix, WCSS |
| `Counts_lite_binarized_*_pattern.tsv` | Per-gene patterns, one per method |
| `Counts_lite_reviewer_summary.tsv` | SSD vs WCSS comparison in a single table |

Note the summary is now **one table covering both methods**, which is what makes the head-to-head comparison readable; running the reviewer separately per method produced one summary each.

**Options:**

| Flag | Effect |
|---|---|
| `-O`, `--output_dir` | Directory for all generated files (default: next to the input matrix) |
| `--methods` | Comma-separated methods to run (default: `SSD,WCSS`) |
| `-p`, `--processors` | Forwarded to each method (default: 1) |
| `--no_review` | Binarize only, skip the behavior reviewer |
| `-s`, `--summary` | Custom path for the summary TSV |
| `-v`, `--verbose` | Detailed processing information |

Methods run one after another, each using all `-p` processors.

**Two caveats worth knowing:**

- The behavior reviewer scans the output directory **recursively**, so binarized matrices left by earlier runs — including ones sitting in subdirectories — are also picked up and included in the comparison. Give each run its own `-O` directory when that matters.
- A relative `-O` resolves against your current directory, **not** against the input matrix. Running from one folder with an input in another puts the results next to neither. The script prints the resolved absolute path in its header and in the final summary; trust that over the string you typed.

**Run log and the next command.** Each run writes `<input_name>_binarizator_log.txt` next to its outputs, recording the inputs, the methods, the settings and every file produced. It ends with the rule-inference command for what it just made, with both paths already filled in — this script was given the counts matrix and it knows which binarized matrix it wrote:

```
Next step - rule inference
------------------------------------------------------------
python3 /path/to/2.Rules_Inference/1.BNI3_Boolean_Rules_Inference.py \
    -i /path/to/Counts_lite.tsv \
    -i_binary /path/to/Counts_lite_binarized_SSD.tsv \
    -o /path/to/output/Boolean_Rules_Inference/
```

The same command is printed to the console at the end of the run. Paths are absolute so it can be pasted from any directory. It names the SSD matrix when both methods ran; swap in the WCSS one if the behavior reviewer favours it.

One thing the command cannot know: the evaluator later enumerates 2^N states, so the network has to be cut to roughly 15–20 genes first. If you have not selected genes yet, point `-i_binary` at that subset instead. Both the log and the console say so.

## 1. SSD (Short Series Discretization)
The **SSD** algorithm (`BNI3_SSD.py`) binarizes gene expression data by iteratively evaluating the distance between contiguous expression values to find the most natural biological separation point.

**How it works:**
1. **Sorting**: For each gene, its expression values across all samples are sorted in ascending order.
2. **Distance Analysis**: A distance matrix is computed between all ordered values.
3. **Clustering by Isolation**: The algorithm repeatedly finds and removes the maximum distances between values. This isolates groups of values and extreme single measurements.
4. **Group Assignment**: The isolated values and the remaining connected clusters are assigned to either a "low expression" group (Group 0) or a "high expression" group (Group 1) based on their relative magnitude.
5. **Proximity Assignment**: Any intermediate values left unassigned are assigned to Group 0 or Group 1 depending on which cluster boundary they are geometrically closest to.
6. **Binarization**: Values in Group 0 are converted to `0`, and values in Group 1 are converted to `1`.

This algorithm is highly conservative and excels when the transition between active/inactive states has sharp thresholds and outliers.

## 2. WCSS (Within-Cluster Sum of Squares)
The **WCSS** algorithm (`BNI3_WCSS.py`) binarizes the expression matrix by finding the optimal partition of values into exactly two clusters that minimize internal variance. This is mathematically equivalent to 1-dimensional k-means clustering with $k=2$.

**How it works:**
1. **Sorting**: For each gene, its continuous expression values are sorted.
2. **Split Evaluation**: The algorithm iterates over every possible splitting index to divide the array into two clusters: Cluster 1 (lower values) and Cluster 2 (higher values).
3. **WCSS Calculation**: For each split, it calculates the Within-Cluster Sum of Squares (WCSS). WCSS is the sum of squared Euclidean distances from each data point to its respective cluster's mean.
   $$WCSS = \sum (x_i - \mu_1)^2 + \sum (y_i - \mu_2)^2$$
4. **Optimization**: The algorithm identifies the splitting index that produces the lowest overall WCSS.
5. **Binarization**: All values belonging to Cluster 1 are discretized to `0`, while those in Cluster 2 become `1`.

WCSS approaches the problem statistically and ensures that the separation minimizes noise within the "on" and "off" states, maximizing the difference between the two states.

## 3. Behavior Reviewer
The **Behavior Reviewer** (`BNI3_behavior_reviewer.py`) is an analytical tool to evaluate the output of the binarization algorithms. It analyzes binarized expression matrices to extract distinct gene expression patterns and their frequencies, allowing you to quickly compare the robustness of methods like SSD vs WCSS.

**How it works:**
1. **Pattern Extraction**: It scans all binarized TSV files in a given directory, ignoring metadata columns (like `ID`).
2. **Frequency Count**: It concatenates the binary sequence for each gene across all samples to form a "pattern" string (e.g., `00110`), and calculates how often each unique pattern appears.
3. **Summary Generation**: It creates a comprehensive summary table detailing total genes, unique patterns, distinct patterns, and max frequency for each method. 
4. **Method Comparison**: It provides statistical summaries and head-to-head comparisons between binarization algorithms directly in the terminal output.

**Usage:**
Run the script from the command line by providing an input directory containing the binarized matrices (`bin_*.tsv`) and an output summary file path:
```bash
python3 BNI3_behavior_reviewer.py -i ./Example -o ./summary_patterns.tsv
```
This command will process all binarized datasets in `./Example`, save individual gene patterns for each dataset, and output an aggregate summary to `./summary_patterns.tsv` alongside detailed terminal analytics.
