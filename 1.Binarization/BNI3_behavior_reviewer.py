#!/usr/bin/env python3
"""
Behavior Reviewer script
Analyzes binarized expression matrices to extract expression patterns and frequencies.
"""

import os
import argparse
import sys
import pandas as pd

def process_binarized_file(input_file):
    """
    Processes a binarized file and returns patterns and their frequencies.
    """
    df = pd.read_csv(input_file, sep="\t")

    # Detect ID column if it exists
    gene_cols = [c for c in df.columns if c.lower() != "id"]

    patterns = {}
    for gene in gene_cols:
        sequence = ''.join(df[gene].astype(str).tolist())
        patterns[gene] = sequence

    pattern_counts = pd.Series(patterns).value_counts()

    output = pd.DataFrame([
        {
            "Gene": gene,
            "Pattern": pat,
            "Frequency": pattern_counts[pat]
        }
        for gene, pat in patterns.items()
    ])

    output = output.sort_values(
        by=["Frequency", "Pattern"],
        ascending=[True, True]
    )

    base, _ = os.path.splitext(input_file)
    output_file = f"{base}_pattern.tsv"
    output.to_csv(output_file, sep="\t", index=False)

    return output, os.path.basename(input_file)


def main():
    parser = argparse.ArgumentParser(
        description="Analyzes behavior and patterns of binarized gene expression matrices.",
        epilog="Example:\n  python3 BNI3_behavior_reviewer.py -i ./binary_expression_matrix -o ./summary_patterns.tsv",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', '--input_path', required=True, help="File or directory containing binarized matrices")
    parser.add_argument('-o', '--output', required=True, help="Output file for the summary (TSV format)")
    
    args = parser.parse_args()

    input_path = args.input_path
    output_summary = args.output

    all_files = []
    if os.path.isfile(input_path):
        all_files.append(input_path)
    elif os.path.isdir(input_path):
        # === Find binarized files ===
        for root, _, files in os.walk(input_path):
            for file in files:
                file_lower = file.lower()
                if "bin" in file_lower and file_lower.endswith(".tsv") and "_pattern" not in file_lower:
                    all_files.append(os.path.join(root, file))
    else:
        print(f"ERROR: Input '{input_path}' does not exist.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(all_files)} binarized files\n")

    if not all_files:
        print("No binarized files found in the specified directory.")
        sys.exit(0)

    # === Process files ===
    summary_data = []

    for file_path in sorted(all_files):
        filename = os.path.basename(file_path)
        print(f"Processing: {filename}")

        try:
            patterns_df, fname = process_binarized_file(file_path)

            # === ROBUST name parsing ===
            name = fname.replace(".tsv", "")
            method_used = "SSD" if "SSD" in name.upper() else ("WCSS" if "WCSS" in name.upper() else "Unknown")

            # === Statistics ===
            n_genes = len(patterns_df)
            n_unique_patterns = (patterns_df["Frequency"] == 1).sum()
            n_distinct_patterns = patterns_df["Pattern"].nunique()
            max_frequency = patterns_df["Frequency"].max()

            summary_data.append({
                "File": fname,
                "Method": method_used,
                "Total_Genes": n_genes,
                "Unique_Patterns": n_unique_patterns,
                "Distinct_Patterns": n_distinct_patterns,
                "Max_Frequency": max_frequency,
                "Unique_Ratio": round(n_unique_patterns / n_genes, 3) if n_genes > 0 else 0
            })

            print(
                f"  [+] Genes: {n_genes}, "
                f"Unique: {n_unique_patterns}, "
                f"Distinct patterns: {n_distinct_patterns}"
            )

        except Exception as e:
            print(f"  [-] Error processing {filename}: {e}")

        print()

    # === Create summary ===
    summary_df = pd.DataFrame(summary_data)

    if summary_df.empty:
        raise RuntimeError(
            "Summary generation failed. "
            "Please check the format and content of the binarized files."
        )

    summary_df = summary_df.sort_values(by=["File"])

    summary_df.to_csv(output_summary, sep="\t", index=False)

    print(f"Summary saved to: {output_summary}\n")
    print(summary_df.to_string(index=False))

if __name__ == "__main__":
    main()
