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

### 0. Attractorator (`BNI3_Attractorator.py`) — the whole stage in one command

Runs the three scripts below in order: attractors, then the path the observed data takes into them, then both figures.

```bash
python3 BNI3_Attractorator.py \
    -i rules_by_gene_evaluated.tsv \
    -b binarized_matrix.tsv \
    -O results/ --predecessors 2
```

Everything lands in `results/attractors/`. **The `attractors` subdirectory is not optional** — `1.BNI3_Attractors.py` appends it to whatever `-o` it is given, so the Attractorator follows it there and prints the resolved absolute path up front rather than leaving the destination to be discovered.

| Flag | Effect |
|---|---|
| `-i` | rules table (required) |
| `-b` | binarized matrix — **required for the path step**; without it that step is skipped and the other two still run |
| `-O` | parent directory (default: the rules file's own directory) |
| `-m "ABF3:1,MYB44:0"` | mutations, forwarded to every step; the suffix propagates into all eight filenames |
| `--no_path`, `--no_visualization` | cut the chain at either point |
| `--predecessors N`, `--predecessors-per-state K` | forwarded to step 3 |
| `-s`, `-ob`, `-c`, `-p`, `-n`, `--max-iter`, `--svg`, `--heatmap-only`, `--network-only` | forwarded to whichever step owns them |

Steps stream their output live with no timeout, since attractor enumeration is Θ(2^N) and legitimately runs for minutes. After step 1 the script checks that the two files it expects actually exist and aborts naming the path if they do not, rather than letting a naming mismatch surface as a confusing failure two steps later. If a later step fails the run continues, and the summary lists what was produced and what failed.

Each of the three scripts still works standalone.

### 3. Attractors Visualizer (`3.BNI3_Visualize_Attractors.py`)
Creates rich graphical diagram representations and heatmaps to visually interpret the detected attractors and state transitions from `attractors.tsv`.

**Basic Usage:**
```bash
python3 3.BNI3_Visualize_Attractors.py -i Example/attractors.tsv --heatmap --network
```
*Outputs generated:* High quality image formats (`.png` / `.svg`) representing the boolean transition networks and expression state heatmaps.

## Changelog

#### Gene-name sanitization in the matrix readers (2026-09-25)

`2.BNI3_Path_to_Attractors.py` and `3.BNI3_Visualize_Attractors.py` read the binarized matrix by looking up each network gene among the matrix columns:

```python
state = [bool(int(row[g])) if g in matrix_genes else False for g in gene_cols]
```

`gene_cols` comes from the attractors file, where names are already sanitized to valid Python identifiers (`SnRK2_8`), while the binarized matrix keeps the original symbols (`SnRK2.8`). The lookup therefore failed and the gene **defaulted to 0** for every sample. The only warning was behind `-v`, and in `3.BNI3_Visualize_Attractors.py` there was none at all.

The effect was one wrong bit per affected gene in every state read from the matrix. On the Arabidopsis networks this corrupted `trajectory_matrix_attractor_mapping.tsv` for three of eight samples and suppressed the blue "final matrix state" rectangle in the Control trajectory plot, because the corrupted final state no longer matched any state on the simulated trajectory. It also made the Control network look as if it disagreed with the data in `SnRK2_8` when in fact it reproduces the final state exactly.

Both readers now apply the same `re.sub(r'\W|^(?=\d)', '_', name)` used by the inference and evaluation scripts, refuse sanitization collisions (`GEN.1` and `GEN_1` both becoming `GEN_1`), and print the "defaulted to 0" warning to stderr unconditionally.

Reruns of any attractor analysis performed before this date are worth repeating if the matrix contained gene symbols with `.`, `-`, spaces or parentheses.

#### Upstream states around the attractors (2026-09-25)

`3.BNI3_Visualize_Attractors.py --predecessors N` adds states that lead **into** each attractor, in the spirit of a BoolNet state transition graph but without drawing all 2^N nodes.

```bash
python3 3.BNI3_Visualize_Attractors.py \
    -i attractors.tsv -r selected_rules.tsv -b binarized_matrix.tsv \
    --network --predecessors 2 --predecessors-per-state 3
```

**The fan-in, not the cost, is what forces a selection.** Building the full transition table is cheap: each rule is evaluated once against numpy boolean arrays spanning the whole state space, so a 15-gene network takes a second. But the Arabidopsis Control fixed point has **2,559 direct predecessors**, 10,240 within two steps and 24,576 within three — essentially the entire space. Drawing "two steps back" literally is the unreadable figure the option exists to avoid.

So predecessors are selected, `--predecessors-per-state` of them per state per step:

1. **states observed in the binarized matrix are always kept.** They are the only upstream states that are measurements rather than possibilities, and keeping them draws the observed trajectory inside the figure, ringed in blue.
2. the remaining slots go to the states **closest in Hamming distance** to their successor, which read as "flip these few genes and the system still returns here".

Ties break on the packed state code, so the figure is deterministic.

**Layout.** Each attractor gets concentric rings: step *d* back sits at radius `base + d × 1.5`, and every state's angular wedge is subdivided among its own predecessors, so a child stays visibly attached to the state it feeds. Fixed points spread over the full circle; each node of a cycle gets a wedge of `2π / cycle_length` pointing outward. Nodes shrink and fade with depth so the attractor stays dominant, and the aspect ratio is locked to equal so one step back is the same distance in every direction.

**Stating its own scale, in numbers that add up.** A figure showing six states out of 32,768 would, left alone, imply that six is all there is. So **every drawn state that the rest of the basin flows through gets its own hollow node beside it**, on a short straight dotted line, labelled with how many states arrive that way and sized by that count on a log scale.

On the Arabidopsis Control network, three of the seven drawn states receive anything:

```
28,668  into the attractor
 3,325  into Sample_3
   768  into Sample_2
```

The other four are Garden-of-Eden leaves with nothing behind them. Worth noting which ones receive: the two that are not the attractor are both **observed samples**. States picked for Hamming proximity tend to be dead ends, while a measured state sits on a real path and carries its whole history behind it.

The arithmetic closes twice over:

```
28,668 + 3,325 + 768              = 32,761   all the undrawn states
32,761 + 1 attractor + 6 drawn    = 32,768   the basin
```

This works because forward dynamics are deterministic: every state has exactly one forward path, that path must reach the attractor, and the attractor is drawn — so every undrawn state has exactly one *first* drawn state it lands on. Grouping by it partitions the remainder with no overlap. The counts were cross-checked against a brute-force forward simulation of all 32,768 states, one at a time, and agree exactly.

Two earlier attempts are worth recording because they failed in instructive ways. A `+N` beside each drawn state counted its own undrawn *predecessors*; those sets nest inside one another, so they overshot the basin (`32,761 + 2,556 + 1,021` against 32,768). A single shared hollow node with curved routes to every receiver partitioned correctly but needed long curves across the figure, which crossed the rings and each other.

**Placement.** A hollow node gets a **reserved slot in the ring, beside its target's own predecessors**, inside the same angular wedge. The receivers are therefore computed *before* anything is positioned, since reserving that slot changes how the fan is spread. A receiver with no drawn predecessors of its own takes a slot one ring further out, still in its own wedge.

Two placement strategies were tried first and both failed on cycles. Putting the node in the widest free gap around its target sends it *inside* the cycle, because the fan points outward and the cycle edges run tangentially — six hollow nodes then pile on the centre. Searching all directions for the one whose dotted line stays furthest from every other node does not help either: a line coming from outside a cycle state has to cross that state's own fan no matter which way it approaches. Reserving a slot inside the fan is the only placement that is clean by construction rather than by search, and it keeps every dotted line short.

The legend also breaks each basin into its backward layers, which partition it too:

```
A1 (fixed_point): 32,768 states (100.0%)
    = 1 in attractor + 2,559 at 1 back + 7,680 at 2 back + 22,528 deeper
```

Two counts that look similar are not: `basin_layer_sizes()` expands *every* state of each layer, while the per-step counts inside `build_predecessor_layers()` expand only the handful the figure kept, so its "available at step 2" means "predecessors of the three states I drew" (1,024 here), not "states two steps from the attractor" (7,680). Only the first kind appears in the figure.

**Mutants.** A knockout or overexpression reaches this script as a constant rule (`MYB44 -> 0`), and `1.BNI3_Attractors.py` enumerates only the states that respect it — 2^12 rather than 2^13 for one clamped gene. Both the backward walk and the predecessor selection apply the same restriction, so the figure never shows a state where the knocked-out gene is ON, and its basin sizes match the attractors file exactly (1,536 / 512 / 2,048 on the bundled `MYB44_0` example).

**Arrows.** Edges are grouped by the size of the node they point at and each group gets its own margin, because networkx takes one margin per edgelist while the nodes here differ by a factor of six in area — a margin that clears a depth-2 state leaves the arrowhead buried inside the attractor. The dotted routes carry arrowheads too, drawn separately and solid: a dotted linestyle applies to the head as well and breaks it into chevrons. Their connector stops one arrowhead short of the rim rather than at it, or the dotted line runs underneath the head and out through its tip. It is drawn with `arrows=True` and `arrowstyle='-'` even though it has no head of its own: with `arrows=False` networkx falls back to a `LineCollection`, which ignores `min_target_margin` outright and runs the line from centre to centre no matter what margin is passed.

A self-loop is sized from the node's radius, which is given in points, so it is drawn last — the conversion to data units needs the axes scale, and that only exists once the limits and the aspect are settled. Both ends sit on the rim near the top and the arc bulges over; the sign of `rad` decides whether it goes over the node or dips inside it, where it reads as a scribble.

**Node size follows the drawing's scale.** Node area is in points, which are physical, while the layout is in data units, and the predecessor rings stretch the drawing to tens of data units against a capped figure width. A constant size therefore covers a few hundredths of the plot in one figure and half a ring in another — which is what let the states of a long cycle overlap and swallowed the arrows between them. The size is now a fixed fraction of the ring spacing, and a cycle's radius is set by arc length per state rather than by a flat multiple of its length.

**One attractor, one fixed-size box.** The axes used to be set from wherever the nodes happened to land, and the width and height were chosen independently before the aspect was locked to equal. Two runs of the same pipeline therefore came out on wildly different canvases — the one-attractor Control network at 2370 x 4282, the two-attractor Drought network at 3394 x 1597 — which is unusable for panels meant to sit side by side.

Each attractor now gets a nominal box of its own, and the axes are set from those boxes rather than from the node positions. The figure is then sized at a fixed `DATA_UNIT_INCHES` per data unit in both directions, so the scale never drifts and the node size follows from it as a constant instead of being estimated from the figure's extent and clamped.

| | before | after |
|---|---|---|
| Control (1 attractor) | 2370 x 4282 | 1523 x 1692 |
| Drought (2 attractors) | 3394 x 1597 | 3072 x 1689 |

Same height, width proportional to the number of attractors, and the attractor node measures 49–50 px in both. Measured across one, two and three attractors the node stays within a couple of pixels; the leftover variation is the black rim eating into the fill, not the scale.

The cost is deliberate whitespace: the box is square around each attractor, and a sparse fan does not fill it. That is what keeps an attractor the same size whether it is alone or one of five.

`DATA_UNIT_INCHES` and `PLOT_PAD` sit at the top of the script if a figure needs to be scaled as a whole. Runs without `--predecessors` keep the old sizing untouched.

**The heatmap keeps one cell size too.** It had the same problem from the other direction: the figure height grew with the number of attractor states, and the aspect switched from `equal` to `auto` as soon as there was more than one row, so the cells stretched vertically to fill whatever height had been chosen. A one-row network came out as a strip of squares, a two-row network as a strip of tall rectangles, and the two could not be stacked in one figure.

The grid is now laid out in absolute inches — `HEATMAP_CELL_INCHES` per cell in both directions — inside axes placed with `add_axes` rather than a stretched gridspec, and the aspect stays `equal` always. The canvas grows with the number of genes and states instead of the cells growing to fill it. Measured across four networks (1, 2, 8 and 8 states; 10 to 15 genes) a cell is 98 x 98 px in every one.

Margins are sized from the text they have to hold: the row labels set the left margin, the rotated gene names the bottom one, and the legend's longest entry the right one. The legend gets its own vertical band, because with one or two rows it is taller than the grid, and the grid is centred against it.

Its legend carries two things and no more: `Inactive gene`, and one row per attractor with its basin size and share. The blue box that used to mark where the last row of the binarized matrix ends up was dropped from this figure — the network and trajectory plots already make that their subject, and here it put an extra border and an extra legend row on a panel whose only job is to let the gene states be read off and compared between conditions.

**The legend names the attractors and the blue ring** — `Attractor 1 (fixed point): 28,672 states (87.5%)` and `States in the binarized matrix` — and the title is just `Boolean Network Attractors`. The rest is labelled on the canvas itself (each hollow node carries its own count) or belongs in the figure caption; spelling it all out in the legend crowded the plot more than it explained it.

A royal-blue ring means one thing anywhere in the figure: **that state was read from the binarized matrix**. An attractor state that was measured gets one, an upstream state that was measured gets one, and a state that was only inferred never does. It previously meant two things at once — "observed" on an upstream node and "this is where the system ends up" on an attractor node — which is why it could not be given a legend row without the row being false for half the nodes it described. The destination of the trajectory is the subject of the trajectory figure, and the inferred one is in `trajectory_matrix_attractor_mapping.tsv`.

The ring needs only `-b/--binarized_matrix`, not `--predecessors`: in the plain diagram it marks whichever attractors were themselves observed. It is read against an attractor's own fill, so it is sharpest on the warm hues and softest on the blue one.

The legend sits below the axes whenever predecessors are drawn: slots are reserved against other *nodes*, which cannot see a legend box, so a legend inside the axes ends up with a node on top of it.

Verified across four networks and seven attractors, cycles and fixed points, wild type and mutant: in every case the hollow counts plus the drawn states sum to the basin.

**Outputs.** Alongside the figure, `*_network_predecessors.tsv` lists every drawn state with its distance from the attractor, whether it was observed, and its gene values — the figure's claim in a form that can be checked.

**Limits.** Needs `-r/--rules_file`, since the attractors file records where the system ends up, not how states map onto one another. Refuses networks above 24 genes, where the 2^N table stops being cheap. With `-v`, the per-step line reports how many predecessors were drawn out of how many existed, which is the honest measure of how much the figure leaves out:

```
A1 step 1 back: 3 drawn of 2,559 available predecessors
A1 step 2 back: 3 drawn of 1,024 available predecessors
```

Worth knowing about these networks: **32,716 of the 32,768 states have no predecessor at all** (Garden-of-Eden states), and mean in-degree is 1.0 while the maximum is 2,560. The state graph is a very shallow, very wide funnel, so branches often stop short of the requested depth — that is the network's structure, not a truncation.
