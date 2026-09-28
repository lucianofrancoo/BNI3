#!/usr/bin/env python3
"""
BNI3 Binarizator
Runs the whole binarization stage from a single command: binarizes one counts
matrix with SSD and WCSS, then reviews both results with the behavior reviewer.

Each step is still available as a standalone script; this one only orchestrates
them so the three do not have to be launched by hand.

Author: Luciano
"""

import argparse
import os
import subprocess
import sys
import time


# Binarization methods this script can run: name -> script filename.
# Both scripts share the same CLI contract (-i input, -o output, -v, -p).
METHODS = {
    'SSD':  'BNI3_SSD.py',
    'WCSS': 'BNI3_WCSS.py',
}

REVIEWER_SCRIPT = 'BNI3_behavior_reviewer.py'


def log_message(message, verbose):
    """Print message only if verbose is enabled"""
    if verbose:
        print(message)


def resolve_script(script_name):
    """
    Locate a pipeline script next to this file.

    Returns the absolute path, or None if it is missing.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, script_name)
    return path if os.path.exists(path) else None


def run_binarization(method, input_file, output_file, processors, verbose):
    """
    Run one binarization method on the counts matrix.

    Output is streamed rather than captured so progress stays visible on large
    matrices, and no timeout is applied.

    Args:
        method (str): Key of METHODS ('SSD' or 'WCSS')
        input_file (str): Counts matrix TSV
        output_file (str): Where to write the binarized matrix
        processors (int): Forwarded to the method as -p
        verbose (bool): Forwarded to the method as -v

    Returns:
        bool: True if the method completed successfully
    """
    script = resolve_script(METHODS[method])

    if script is None:
        print(f"WARNING: {method} script not found at "
              f"{os.path.join(os.path.dirname(os.path.abspath(__file__)), METHODS[method])}")
        return False

    cmd = ["python3", script, "-i", input_file, "-o", output_file,
           "-p", str(processors)]
    if verbose:
        cmd.append("-v")

    try:
        print(f"\n{'='*60}")
        print(f">> Running {method} binarization...")
        print('='*60)
        result = subprocess.run(cmd)

        if result.returncode == 0:
            print(f"{method} completed -> {output_file}")
            return True

        print(f"{method} failed with return code {result.returncode}")
        return False

    except KeyboardInterrupt:
        print(f"\n{method} interrupted by user.")
        raise
    except Exception as e:
        print(f"Error running {method}: {str(e)}")
        return False


def run_behavior_review(review_path, summary_file, verbose):
    """
    Run the behavior reviewer over the binarized matrices.

    Args:
        review_path (str): Directory holding the binarized matrices
        summary_file (str): Where to write the comparison summary
        verbose (bool): Enable verbose logging in this script

    Returns:
        bool: True if the reviewer completed successfully
    """
    script = resolve_script(REVIEWER_SCRIPT)

    if script is None:
        print(f"WARNING: Behavior reviewer not found at "
              f"{os.path.join(os.path.dirname(os.path.abspath(__file__)), REVIEWER_SCRIPT)}")
        print(f"Please run it manually with: "
              f"python3 {REVIEWER_SCRIPT} -i {review_path} -o {summary_file}")
        return False

    cmd = ["python3", script, "-i", review_path, "-o", summary_file]

    try:
        print(f"\n{'='*60}")
        print(f">> Running behavior reviewer...")
        print('='*60)
        # The reviewer scans the whole directory for binarized matrices, so any
        # file from an earlier run that is still there is included in the
        # comparison. Use -O to isolate a run in its own directory.
        log_message(f"Reviewing binarized matrices in: {review_path}", verbose)
        result = subprocess.run(cmd)

        if result.returncode == 0:
            return True

        print(f"Behavior reviewer failed with return code {result.returncode}")
        return False

    except KeyboardInterrupt:
        print("\nBehavior reviewer interrupted by user.")
        raise
    except Exception as e:
        print(f"Error running behavior reviewer: {str(e)}")
        return False


def sibling_script(*parts):
    """Absolute path to another pipeline script, relative to this one."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(root, *parts)


def next_step_command(counts_matrix, binarized_matrix, output_dir):
    """
    The rule-inference command for what this run just produced.

    Written with absolute paths so it can be pasted from any directory, and with
    the two inputs already filled in: this script was given the counts matrix and
    it knows which binarized matrix it wrote.
    """
    return (
        f"python3 {sibling_script('2.Rules_Inference', '1.BNI3_Boolean_Rules_Inference.py')} \\\n"
        f"    -i {os.path.abspath(counts_matrix)} \\\n"
        f"    -i_binary {binarized_matrix} \\\n"
        f"    -o {os.path.join(output_dir, 'Boolean_Rules_Inference')}/"
    )


def write_run_log(path, args, output_dir, produced, failed, summary_file,
                  review_ok, total_time, command):
    """Record what this run did and how to continue it."""
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write("BNI3 Binarizator run log\n")
        handle.write("=" * 60 + "\n")
        handle.write(f"generated       : {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        handle.write(f"counts matrix   : {os.path.abspath(args.input)}\n")
        handle.write(f"output directory: {output_dir}\n")
        handle.write(f"methods         : {', '.join(produced) if produced else 'none'}\n")
        handle.write(f"processors      : {args.processors}\n")
        handle.write(f"behavior review : {'skipped' if args.no_review else 'enabled'}\n")
        handle.write(f"total time      : {total_time:.2f} s\n")
        handle.write("\nFiles produced\n")
        handle.write("-" * 60 + "\n")
        for method, produced_path in produced.items():
            handle.write(f"  {method:<6} -> {produced_path}\n")
        if review_ok:
            handle.write(f"  review -> {summary_file}\n")
            handle.write("            (plus one *_pattern.tsv per binarized matrix)\n")
        if failed:
            handle.write(f"  failed: {', '.join(failed)}\n")
        handle.write("\nNext step - rule inference\n")
        handle.write("-" * 60 + "\n")
        handle.write(command + "\n")
        handle.write(
            "\nNote: the evaluator enumerates 2^N states, so the network has to be\n"
            "cut down to roughly 15-20 genes first. Point -i_binary at the selected\n"
            "subset rather than the full matrix if you have not already.\n")


def run_binarizator(args):
    """Run the full binarization stage and report what was produced."""
    start_time = time.time()

    input_stem = os.path.splitext(os.path.basename(args.input))[0]

    # Default output directory: alongside the input matrix.
    # Resolved to an absolute path so every path printed below is unambiguous —
    # a relative -O resolves against the current directory, not against the input,
    # which is easy to misread when the input lives somewhere else.
    output_dir = os.path.abspath(args.output_dir) if args.output_dir \
        else os.path.dirname(os.path.abspath(args.input))
    os.makedirs(output_dir, exist_ok=True)

    methods = [m.strip().upper() for m in args.methods.split(',') if m.strip()]
    unknown = [m for m in methods if m not in METHODS]
    if unknown:
        print(f"ERROR: Unknown method(s): {', '.join(unknown)}. "
              f"Available: {', '.join(METHODS)}", file=sys.stderr)
        sys.exit(1)

    print('='*60)
    print("BNI3 BINARIZATOR")
    print('='*60)
    print(f"Input matrix:     {args.input}")
    print(f"Output directory: {output_dir}")
    print(f"Methods:          {', '.join(methods)}")
    print(f"Processors:       {args.processors}")
    print(f"Behavior review:  {'skipped' if args.no_review else 'enabled'}")

    # Binarize with each requested method
    produced = {}
    failed = []
    for method in methods:
        output_file = os.path.join(
            output_dir, f"{input_stem}_binarized_{method}.tsv")
        if run_binarization(method, args.input, output_file,
                            args.processors, args.verbose):
            produced[method] = output_file
        else:
            failed.append(method)

    if not produced:
        print(f"\n{'='*60}")
        print("STATUS: ALL BINARIZATION METHODS FAILED")
        print('='*60)
        sys.exit(1)

    # Review the binarized matrices unless disabled
    summary_file = os.path.abspath(args.summary) if args.summary else os.path.join(
        output_dir, f"{input_stem}_reviewer_summary.tsv")
    review_ok = False

    if args.no_review:
        print("\nBehavior review skipped. You can run it manually with:")
        print(f"python3 {REVIEWER_SCRIPT} -i {output_dir} -o {summary_file}")
    else:
        review_ok = run_behavior_review(output_dir, summary_file, args.verbose)
        if not review_ok:
            print("\nBehavior review failed. You can run it manually with:")
            print(f"python3 {REVIEWER_SCRIPT} -i {output_dir} -o {summary_file}")

    # Final summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print("BINARIZATION SUMMARY")
    print('='*60)
    for method, path in produced.items():
        print(f"  {method:<6} -> {path}")
    if review_ok:
        print(f"  Review -> {summary_file}")
        print(f"            (plus one *_pattern.tsv per binarized matrix)")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    print(f"  Total time: {total_time:.2f} seconds")
    print('='*60)

    if failed:
        print(f"STATUS: {len(failed)} METHOD(S) FAILED")
        sys.exit(1)

    print("STATUS: BINARIZATION COMPLETED SUCCESSFULLY")

    command = next_step_command(
        args.input,
        produced.get('SSD', list(produced.values())[0]),
        output_dir)

    log_path = os.path.join(output_dir, f"{input_stem}_binarizator_log.txt")
    write_run_log(log_path, args, output_dir, produced, failed, summary_file,
                  review_ok, total_time, command)
    print(f"  Log    -> {log_path}")

    print(f"\nNext step — rule inference:")
    print(command)
    print("\nNote: the evaluator enumerates 2^N states, so cut the network down to")
    print("roughly 15-20 genes first and point -i_binary at that subset.")


def main():
    """Main function with argument parsing"""
    parser = argparse.ArgumentParser(
        description='Run the full BNI3 binarization stage: SSD + WCSS + behavior review',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Binarize with both methods and review the results
  python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv

  # Write everything to a separate directory
  python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv -O results/

  # Only one method
  python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv --methods SSD

  # Binarize only, no behavior review
  python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv --no_review

  # Use 8 processors, verbose
  python3 BNI3_Binarizator.py -i Example/Counts_lite.tsv -p 8 -v

Output files (for an input named Counts_lite.tsv):
  Counts_lite_binarized_SSD.tsv            binarized matrix, SSD
  Counts_lite_binarized_WCSS.tsv           binarized matrix, WCSS
  Counts_lite_binarized_*_pattern.tsv      per-gene patterns (behavior reviewer)
  Counts_lite_reviewer_summary.tsv         SSD vs WCSS comparison

Notes:
  - Methods run one after another, each using all -p processors.
  - The behavior reviewer scans the whole output directory, so binarized
    matrices left there by earlier runs also appear in the comparison.
    Use -O to give a run its own directory.
  - Each underlying script (BNI3_SSD.py, BNI3_WCSS.py,
    BNI3_behavior_reviewer.py) still works standalone.
        """
    )

    # Required parameters
    required = parser.add_argument_group('Required parameters')
    required.add_argument('-i', '--input', type=str, required=True,
                         help='Input counts matrix (expression matrix in TSV format)')

    # Optional parameters
    optional = parser.add_argument_group('Optional parameters')
    optional.add_argument('-O', '--output_dir', type=str, default=None,
                         help='Directory for all generated files '
                              '(default: the directory holding the input matrix)')
    optional.add_argument('--methods', type=str, default='SSD,WCSS',
                         help='Comma-separated binarization methods to run '
                              '(default: SSD,WCSS)')
    optional.add_argument('-p', '--processors', type=int, default=1,
                         help='Number of processors forwarded to each method (default: 1)')
    optional.add_argument('-v', '--verbose', action='store_true',
                         help='Show detailed processing information')

    # Behavior review
    review = parser.add_argument_group('Behavior review parameters')
    review.add_argument('--no_review', action='store_true',
                       help='Binarize only, skipping the behavior reviewer')
    review.add_argument('-s', '--summary', type=str, default=None,
                       help='Path for the reviewer summary TSV '
                            '(default: <input_name>_reviewer_summary.tsv in the output directory)')

    parser.add_argument('--version', action='version', version='BNI3 Binarizator v1.0')

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"ERROR: Input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)

    try:
        run_binarizator(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
