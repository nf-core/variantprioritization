#!/usr/bin/env python3

import argparse
from pysam import VariantFile
import os
import sys

def reformat_pon(pon_file, out):
    with VariantFile(pon_file) as fr:
        header = fr.header
        header.info.add("PANEL_OF_NORMALS", number=0, type="Flag",
                        description="Overlap with germline call among panel of normals")
        with VariantFile("tmp_.vcf", "w", header=header) as fw:
            for record in fr:
                record.info["set"] = "PANEL_OF_NORMALS"
                fw.write(record)

    os.system(f"mv tmp_.vcf {out}")
    os.system(f"bgzip -f {out}")
    os.system(f"tabix -f {out}.gz")

def main():
    parser = argparse.ArgumentParser(description="Reformat Panel of Normals (PoN) VCF file and annotate with PANEL_OF_NORMALS flag.")
    parser.add_argument("-i", "--input", required=True, help="Input PoN VCF file (can be .vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file name (without .gz)")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist.")
        sys.exit(1)

    reformat_pon(args.input, args.output)

if __name__ == "__main__":
    main()
