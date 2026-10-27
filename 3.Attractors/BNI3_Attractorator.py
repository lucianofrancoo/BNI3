#!/usr/bin/env python3
"""
BNI3 Attractorator
Runs the whole attractor stage from a single command: finds the attractors of a
rule set, traces the path the observed data takes into them, and draws both
figures.

Each step is still available as a standalone script; this one only orchestrates
them so the three do not have to be launched by hand.

Author: Luciano
"""

import argparse
import os
import subprocess
import sys
import time


ATTRACTORS_SCRIPT = '1.BNI3_Attractors.py'
PATH_SCRIPT = '2.BNI3_Path_to_Attractors.py'
VISUALIZE_SCRIPT = '3.BNI3_Visualize_Attractors.py'


def resolve_script(script_name):
    """
    Locate a pipeline script next to this file.

    Returns the absolute path, or None if it is missing.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, script_name)
    return path if os.path.exists(path) else None


def mutation_suffix(mutations_str):
    """
    Rebuild the filename suffix 1.BNI3_Attractors.py appends for mutations.

    Mirrors generate_mutation_suffix() there: genes sorted by name, each written
    as "_GENE_VALUE". Kept in step with that function so this script knows what
    the attractors file will be called without having to guess or glob for it.
    """
    if not mutations_str:
        return ""

    parts = []
    for mutation in mutations_str.split(','):
        if not mutation.strip():
            continue
        gene, _, value = mutation.partition(':')
        parts.append((gene.strip(), value.strip()))

    if not parts:
        return ""

    return "_" + "_".join(f"{gene}_{value}" for gene, value in sorted(parts))


def run_step(title, cmd):
    """
    Run one pipeline step, streaming its output.

    Output is streamed rather than captured so progress stays visible, and no
    timeout is applied: attractor enumeration is Theta(2^N) and legitimately
    runs for minutes on larger networks.

    Returns:
        bool: True if the step completed successfully
    """
    print(f"\n{'='*60}")
    print(f">> {title}")
    print('='*60)

    try:
        result = subprocess.run(cmd)
        if result.returncode == 0:
            return True
        print(f"\n{title} failed with return code {result.returncode}",
              file=sys.stderr)
        return False

    except KeyboardInterrupt:
        print(f"\n{title} interrupted by user.", file=sys.stderr)
        raise
    except Exception as exc:
        print(f"\nError running {title}: {exc}", file=sys.stderr)
        return False


def run_attractorator(args):
    """Run the full attractor stage and report what was produced."""
    start_time = time.time()

    # Resolved to absolute paths so every path printed below is unambiguous: a
    # relative -O resolves against the current directory, not against the input,
    # which is easy to misread when the input lives somewhere else.
    base_dir = os.path.abspath(args.output_dir) if args.output_dir \
        else os.path.dirname(os.path.abspath(args.input))

    # 1.BNI3_Attractors.py always appends "attractors" to the -o it is given and
    # writes there, so that subdirectory — not -O itself — is where steps 2 and 3
    # must read from and write to. Printed in full below so the destination is
    # never a surprise.
    output_dir = os.path.join(base_dir, "attractors")
    os.makedirs(output_dir, exist_ok=True)

    matrix = os.path.abspath(args.binarized_matrix) if args.binarized_matrix else None

    suffix = mutation_suffix(args.mutations)
    attractors_file = os.path.join(output_dir, f"attractors{suffix}.tsv")
    rules_file = os.path.join(output_dir, f"selected_rules{suffix}.tsv")

    print('='*60)
    print("BNI3 ATTRACTORATOR")
    print('='*60)
    print(f"Rules file:       {os.path.abspath(args.input)}")
    print(f"Output directory: {output_dir}")
    print(f"                  (step 1 appends 'attractors' to -O)")
    print(f"Binarized matrix: {matrix if matrix else 'not given'}")
    if args.mutations:
        print(f"Mutations:        {args.mutations}")
    print(f"Path to attractors: {'skipped' if args.no_path else 'enabled'}")
    print(f"Visualization:      {'skipped' if args.no_visualization else 'enabled'}")

    # The path step needs the observed data; without it there is no trajectory
    # to trace. Say so once, up front, rather than failing three steps in.
    if not args.no_path and matrix is None:
        print("\nNOTE: no binarized matrix given (-b), so the path-to-attractors "
              "step will be skipped. The attractors and the figures do not need it, "
              "but the observed trajectory cannot be traced without it.")

    produced = []
    failed = []

    # ---- Step 1: attractors -------------------------------------------------
    script = resolve_script(ATTRACTORS_SCRIPT)
    if script is None:
        print(f"ERROR: {ATTRACTORS_SCRIPT} not found next to this script.",
              file=sys.stderr)
        sys.exit(1)

    cmd = ["python3", script, "-i", os.path.abspath(args.input),
           "-o", base_dir, "--max-iter", str(args.max_iterations)]
    if args.processes is not None:
        cmd.extend(["-n", str(args.processes)])
    if args.mutations:
        cmd.extend(["-m", args.mutations])
    if args.criteria:
        cmd.extend(["-c", args.criteria])
    if args.positions:
        cmd.extend(["-p", args.positions])
    if args.verbose:
        cmd.append("-v")

    if not run_step("Step 1/3 — finding attractors", cmd):
        print("\nSTATUS: ATTRACTOR SEARCH FAILED — nothing downstream can run.",
              file=sys.stderr)
        sys.exit(1)

    # Everything after this reads these two files, so a wrong guess about their
    # names would surface as a confusing "file not found" three steps later.
    for path, label in ((attractors_file, 'attractors'), (rules_file, 'selected rules')):
        if not os.path.exists(path):
            print(f"\nERROR: expected the {label} file at {path}, but it is not "
                  f"there. Step 1 may name its output differently than this "
                  f"script expects.", file=sys.stderr)
            sys.exit(1)
    produced.extend([attractors_file, rules_file])

    # ---- Step 2: path to attractors ----------------------------------------
    if args.no_path or matrix is None:
        if args.no_path:
            print("\nPath-to-attractors skipped. You can run it manually with:")
            print(f"python3 {PATH_SCRIPT} -a {attractors_file} "
                  f"-r {rules_file} -b <binarized_matrix.tsv> -o {output_dir}")
    else:
        script = resolve_script(PATH_SCRIPT)
        if script is None:
            print(f"WARNING: {PATH_SCRIPT} not found next to this script.")
            failed.append('path to attractors')
        else:
            cmd = ["python3", script, "-a", attractors_file, "-r", rules_file,
                   "-b", matrix, "-o", output_dir,
                   "--max-steps", str(args.max_iterations)]
            if args.output_base:
                cmd.extend(["-ob", args.output_base])
            if args.initial_state:
                cmd.extend(["-s", args.initial_state])
            if args.verbose:
                cmd.append("-v")

            if run_step("Step 2/3 — tracing the path into the attractors", cmd):
                base = args.output_base or 'trajectory'
                produced.extend([
                    os.path.join(output_dir, f"{base}.tsv"),
                    os.path.join(output_dir, f"{base}_matrix_attractor_mapping.tsv"),
                    os.path.join(output_dir, f"{base}_trajectory.png"),
                ])
            else:
                failed.append('path to attractors')

    # ---- Step 3: visualization ---------------------------------------------
    if args.no_visualization:
        print("\nVisualization skipped. You can run it manually with:")
        print(f"python3 {VISUALIZE_SCRIPT} -i {attractors_file} "
              f"-r {rules_file} -o {output_dir}")
    else:
        script = resolve_script(VISUALIZE_SCRIPT)
        if script is None:
            print(f"WARNING: {VISUALIZE_SCRIPT} not found next to this script.")
            failed.append('visualization')
        else:
            cmd = ["python3", script, "-i", attractors_file, "-r", rules_file,
                   "-o", output_dir]
            if matrix:
                cmd.extend(["-b", matrix])
            if args.predecessors:
                cmd.extend(["--predecessors", str(args.predecessors),
                            "--predecessors-per-state", str(args.predecessors_per_state)])
            if args.svg:
                cmd.append("--svg")
            if args.heatmap_only:
                cmd.append("--heatmap-only")
            elif args.network_only:
                cmd.append("--network-only")
            if args.verbose:
                cmd.append("-v")

            if run_step("Step 3/3 — drawing the attractor figures", cmd):
                base = f"attractors_visualization{suffix}"
                if not args.network_only:
                    produced.append(os.path.join(output_dir, f"{base}_heatmap.png"))
                if not args.heatmap_only:
                    produced.append(os.path.join(output_dir, f"{base}_network.png"))
                    if args.predecessors:
                        produced.append(os.path.join(
                            output_dir, f"{base}_network_predecessors.tsv"))
            else:
                failed.append('visualization')

    # ---- Summary ------------------------------------------------------------
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print("ATTRACTOR STAGE SUMMARY")
    print('='*60)
    for path in produced:
        mark = ' ' if os.path.exists(path) else '?'
        print(f" {mark} {path}")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    print(f"  Total time: {total_time:.2f} seconds")
    print('='*60)

    if failed:
        print(f"STATUS: {len(failed)} STEP(S) FAILED")
        sys.exit(1)

    print("STATUS: ATTRACTOR STAGE COMPLETED SUCCESSFULLY")


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Run the full BNI3 attractor stage: attractors + path + visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Attractors, path and both figures
  python3 BNI3_Attractorator.py -i rules_by_gene_evaluated.tsv \\
      -b binarized_matrix.tsv -O attractors/

  # Add upstream states to the network figure
  python3 BNI3_Attractorator.py -i rules_by_gene_evaluated.tsv \\
      -b binarized_matrix.tsv -O attractors/ --predecessors 2

  # Knockout / overexpression
  python3 BNI3_Attractorator.py -i rules.tsv -b matrix.tsv -m "ABF3:1,MYB44:0"

  # Attractors and figures only, no trajectory
  python3 BNI3_Attractorator.py -i rules.tsv --no_path

Output files (for a run without mutations):
  attractors.tsv                             attractor states, cycles and basins
  selected_rules.tsv                         the rule set that was analysed
  trajectory.tsv                             simulated path into the attractor
  trajectory_matrix_attractor_mapping.tsv    attractor assignment per sample
  trajectory_trajectory.png                  trajectory heatmap
  attractors_visualization_heatmap.png       attractor states heatmap
  attractors_visualization_network.png       attractor transition diagram

With -m the mutation suffix appears in every name, e.g. attractors_ABF3_1.tsv.

Notes:
  - -b is required for the path step. Without it that step is skipped and the
    other two still run.
  - Each underlying script (1.BNI3_Attractors.py, 2.BNI3_Path_to_Attractors.py,
    3.BNI3_Visualize_Attractors.py) still works standalone.
        """
    )

    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input', type=str, required=True,
                         help='Boolean rules TSV (rules_by_gene_evaluated.tsv, '
                              'rules_by_gene.tsv or selected_rules.tsv)')

    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-b', '--binarized_matrix', type=str, default=None,
                         help='Binarized expression matrix TSV. Required for the '
                              'path step; also marks the observed final state in '
                              'the figures')
    optional.add_argument('-O', '--output_dir', type=str, default=None,
                         help='Directory for all generated files '
                              '(default: the directory holding the rules file)')
    optional.add_argument('-m', '--mutations', type=str, default=None,
                         help='Gene mutations as "GENE1:1,GENE2:0" '
                              '(1=overexpression, 0=knockout). Forwarded to every step')
    optional.add_argument('-n', '--processes', type=int, default=None,
                         help='Worker processes for the attractor search '
                              '(default: auto-detect)')
    optional.add_argument('--max-iter', '--max_iterations', dest='max_iterations',
                         type=int, default=1000,
                         help='Maximum simulation steps per trajectory (default: 1000)')
    optional.add_argument('-c', '--criteria', type=str, default=None,
                         choices=['top_position', 'custom_positions'],
                         help='Rule selection criteria (default: top_position)')
    optional.add_argument('-p', '--positions', type=str, default=None,
                         help='Comma-separated positions for custom_positions, '
                              'e.g. "1,1,2,3,1"')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')

    path_group = parser.add_argument_group('Path to attractors parameters')
    path_group.add_argument('--no_path', action='store_true',
                           help='Skip the path-to-attractors step')
    path_group.add_argument('-s', '--initial_state', type=str, default=None,
                           help='Initial state as a binary string or comma-separated '
                                'active gene names (default: first row of the matrix)')
    path_group.add_argument('-ob', '--output_base', type=str, default=None,
                           help='Base name for the trajectory files (default: trajectory)')

    viz_group = parser.add_argument_group('Visualization parameters')
    viz_group.add_argument('--no_visualization', action='store_true',
                          help='Skip the figures')
    viz_group.add_argument('--predecessors', type=int, default=0, metavar='STEPS',
                          help='Draw states leading INTO each attractor, walking this '
                               'many update steps back (0 = off, 2 is a good start)')
    viz_group.add_argument('--predecessors-per-state', type=int, default=3, metavar='K',
                          help='Predecessors kept per state per step (default: 3)')
    viz_group.add_argument('--heatmap-only', action='store_true',
                          help='Draw only the attractor states heatmap')
    viz_group.add_argument('--network-only', action='store_true',
                          help='Draw only the attractor transition diagram')
    viz_group.add_argument('--svg', action='store_true',
                          help='Also save the figures as SVG')

    parser.add_argument('--version', action='version', version='BNI3 Attractorator v1.0')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: Rules file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)

    if args.binarized_matrix and not os.path.exists(args.binarized_matrix):
        print(f"ERROR: Binarized matrix '{args.binarized_matrix}' does not exist.",
              file=sys.stderr)
        sys.exit(1)

    if args.heatmap_only and args.network_only:
        print("ERROR: --heatmap-only and --network-only are mutually exclusive.",
              file=sys.stderr)
        sys.exit(1)

    try:
        run_attractorator(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
