#!/usr/bin/env python3

import os
import glob
import re
import pandas as pd
import argparse


def intersect_variants(sample):

    sample_id = sample
    suffixes = (".cns", ".tbi")
    r = re.compile(f"{sample_id}*")
    sample_files = os.listdir("./")
    sample_files = list(filter(r.match, sample_files))
    sample_files = [file for file in sample_files if not file.endswith(suffixes)]
    print(sample_files)

    tool_names = {}
    for idx, file in enumerate(sample_files):
        tool = file.split(".")[2]  # change this if you change prefix
        tool_names[idx] = tool

    if len(sample_files) > 1:

        for idx, x in enumerate(sample_files):
            idx = idx + 1  # cant use 0 for -n
            os.system(f'bcftools isec -c none -n={idx} -p {idx} {" ".join(str(x) for x in sample_files)}')

        pattern = "./**/sites.txt"
        fn_size = {}
        file_list = glob.glob(pattern, recursive=True)
        for file in file_list:
            file_size = os.stat(file).st_size
            fn_size[file] = file_size
        # remove sites.txt files that are empty
        fn_size = {key: val for key, val in fn_size.items() if val != 0}
        file_list = list(fn_size.keys())

        li = []
        for filename in file_list:
            df = pd.read_table(filename, sep="\t", header=None, converters={4: str})  # preserve leading zeros
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
            bytes_2_tal = ",".join(bytes_2_tal.values())
            convert_column.append(bytes_2_tal)

        assert len(convert_column) == len(frame), f"bytes to TAL section failed - length of list != length DF"
        frame[4] = convert_column
        # I noticed duplicate rows in the output file during testing. Worrying as I'm not sure how they got there...
        # chr1    3866080 C       T       freebayes
        # chr1    3866080 C       T       freebayes
        frame = frame.drop_duplicates()
        frame.to_csv(f"{sample}_keys.txt", sep="\t", index=None, header=None)

    else:

        os.system(
            f'bcftools view {sample_files[0]} -G -H | awk -v OFS="\t" \'{{print $1, $2, $4, $5, "{tool_names[0]}"}}\' > {sample}_keys.txt'
        )


def main():
    # Argument parsing using argparse
    parser = argparse.ArgumentParser(description="Reformat somatic CNA files for PCGR input.")
    parser.add_argument("-s", "--sample", required=True, help="Sample name (meta.id) for the output.")

    args = parser.parse_args()


    # Call reformat_cna function with arguments
    intersect_variants(args.sample)


if __name__ == "__main__":
    main()
