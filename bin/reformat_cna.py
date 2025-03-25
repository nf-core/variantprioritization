#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import sys

def compute_major_minor(cn):
    """
    Compute nMajor and nMinor based on total copy number (cn).
    """
    if cn % 2 == 0:
        return cn // 2, cn // 2
    else:
        return (cn + 1) // 2, cn // 2

def reformat_cna(cna_file, sample):
    """
    Reformat CNVkit .call.cns file for PCGR allele-specific input.
    :param cna_file: CNVkit .call.cns file.
    :param sample: Sample name for output (not used directly but kept for compatibility).
    """
    # Read input (tab-delimited CNVkit .call.cns)
    df = pd.read_csv(cna_file, sep='\t')
    #tool = guess_tool(df)


    # Validate required columns
    required_cols = {"chromosome", "start", "end", "cn"}
    if not required_cols.issubset(df.columns):
        print(f"Error: Input file must contain columns: {', '.join(required_cols)}")
        sys.exit(1)

    # Ensure CN is rounded and integer
    df['cn'] = df['cn'].round().astype(int)

    # Compute nMajor and nMinor
    df[['nMajor', 'nMinor']] = df['cn'].apply(lambda cn: pd.Series(compute_major_minor(cn)))

    # Convert start to 1-based for PCGR
    df['start'] = df['start'] + 1

    # Prepare output with required PCGR columns
    pcgr_df = df[['chromosome', 'start', 'end', 'nMajor', 'nMinor']]

    # Output filename
    out_file = sample + ".allele_specific.pcgr.tsv"
    pcgr_df.to_csv(out_file, sep='\t', index=False)

    print(f"PCGR allele-specific CNA file written to: {out_file}")

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Convert CNVkit .call.cns to PCGR allele-specific format.")
    parser.add_argument("-i", "--input", required=True, help="CNVkit .call.cns input file (tab-delimited).")
    parser.add_argument("-s", "--sample", required=True, help="Sample name (used for output filename).")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    # Process CNA file
    reformat_cna(args.input, args.sample)

if __name__ == "__main__":
    main()
