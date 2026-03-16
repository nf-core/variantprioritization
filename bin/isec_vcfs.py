#!/usr/bin/env python3

import os
import glob
from pathlib import Path
import pandas as pd
import argparse
import subprocess


def intersect_variants(input_dir, output_file):

    input_path = Path(input_dir)
    if not input_path.exists() or not input_path.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist or is not a directory: {input_dir}"
        )

    sample_files = sorted(
        file.name
        for file in input_path.iterdir()
        if file.is_file()
        and (file.name.endswith(".vcf") or file.name.endswith(".vcf.gz"))
    )
    print(sample_files)

    if not sample_files:
        raise ValueError(f"No VCF files found in input directory: {input_dir}")

    tbi_files = {
        file.name
        for file in input_path.iterdir()
        if file.is_file() and file.name.endswith(".tbi")
    }
    for vcf_name in sample_files:
        if vcf_name.endswith(".vcf.gz") and f"{vcf_name}.tbi" not in tbi_files:
            raise FileNotFoundError(
                f"Missing index for {vcf_name}: expected {vcf_name}.tbi in {input_dir}"
            )

    tool_names = {}
    for idx, file in enumerate(sample_files):
        tool = file.split(".")[2]  # change this if you change prefix
        if tool.split("_")[0] == "strelka":
            tool = "strelka"
        tool_names[idx] = tool

    if len(sample_files) > 1:
        for idx, _x in enumerate(sample_files):
            idx = idx + 1  # cant use 0 for -n
            subprocess.run(
                [
                    "bcftools",
                    "isec",
                    "-c",
                    "none",
                    f"-n={idx}",
                    "-p",
                    str(idx),
                    *sample_files,
                ],
                cwd=input_dir,
                check=True,
            )

        pattern = f"{input_dir}/**/sites.txt"
        fn_size = {}
        file_list = glob.glob(pattern, recursive=True)
        for file in file_list:
            file_size = os.stat(file).st_size
            fn_size[file] = file_size
        # remove sites.txt files that are empty
        fn_size = {key: val for key, val in fn_size.items() if val != 0}
        file_list = sorted(fn_size.keys())

        if not file_list:
            raise ValueError(
                "No non-empty sites.txt files produced by bcftools isec. "
                "Check staged VCF/TBI pairs in --input_dir."
            )

        li = []
        for filename in file_list:
            df = pd.read_table(
                filename, sep="\t", header=None, converters={4: str}
            )  # preserve leading zeros
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)

        # loop over the 4th column containing bytes '0101' etc
        convert_column = []
        for byte_str in frame[4]:
            # print(byte_str)
            # convert 0101 to itemized list
            code = [x for x in byte_str]
            # init list to match 1's in bytestring to corresponding tool name
            grab_index = []
            for idx, val in enumerate(code):
                if val != "0":
                    grab_index.append(int(idx))
            bytes_2_tal = {k: tool_names[k] for k in grab_index if k in tool_names}
            bytes_2_tal = ",".join(
                sorted(bytes_2_tal.values())
            )  # Sort tool names alphabetically
            convert_column.append(bytes_2_tal)

        assert len(convert_column) == len(frame), (
            "bytes to TAL section failed - length of list != length DF"
        )
        frame[4] = convert_column
        # I noticed duplicate rows in the output file during testing. Worrying as I'm not sure how they got there...
        # chr1    3866080 C       T       freebayes
        # chr1    3866080 C       T       freebayes
        frame = frame.drop_duplicates()
        # Sort by chromosome and position for reproducible output
        frame = frame.sort_values(by=[0, 1]).reset_index(drop=True)
        frame.to_csv(output_file, sep="\t", index=None, header=None)

    else:
        single_vcf = input_path / sample_files[0]
        result = subprocess.run(
            ["bcftools", "view", str(single_vcf), "-G", "-H"],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        with open(output_file, "w") as out_fh:
            for line in result.stdout.splitlines():
                fields = line.split("\t")
                out_fh.write(
                    f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}\t{tool_names[0]}\n"
                )


def main():
    # Argument parsing using argparse
    parser = argparse.ArgumentParser(
        description="Reformat somatic CNA files for PCGR input."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing staged input VCF and TBI files.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output key mapping filename."
    )

    args = parser.parse_args()

    # Call reformat_cna function with arguments
    intersect_variants(args.input_dir, args.output)


if __name__ == "__main__":
    main()
