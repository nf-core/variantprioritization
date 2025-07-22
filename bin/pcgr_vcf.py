#!/usr/bin/env python3

import os
import re
import pandas as pd
from pysam import VariantFile
import argparse

# same as before, just make two similar functions for tumor-normal, tumor-only.

def tumor_normal(sample):
    # pass meta.id to the script to remove .command* etc files from file search:
    sample_id = sample
    suffixes = (".cns", ".tbi")
    r = re.compile(f"{sample_id}*")
    sample_files = os.listdir("./")
    sample_files = list(filter(r.match, sample_files))
    sample_files = [file for file in sample_files if not file.endswith(suffixes)]

    # there will be a sample_keys.txt file in the work directory
    # load it as key_Df and remove from the VCF list before running vcf2tsvpy
    key = "".join([file for file in sample_files if file.endswith(".txt")])
    key_df = pd.read_table(key, header=None, sep="\t")
    key_df.index = pd.MultiIndex.from_arrays(
        key_df.values.T[
            (0, 1, 2, 3),
        ]
    )
    sample_files.remove(key)

    # stage fields we are interested in capturing
    fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "NDP", "ADN", "NAF", "TDP", "ADT", "TAF"]

    # create dict of dataframes containing dataframes of the abpve fields.
    # strategy is to average TDP,NDP, TAF,NAF and append TAL from $meta.id_keys.txt

    vcf_dict = {}
    for idx, x in enumerate(sample_files):
        idx = idx + 1
        os.system(f"vcf2tsvpy --input_vcf {x} --out_tsv {idx}.tmp --skip_genotype_data")
        os.system(f"tail -n +2 {idx}.tmp > {idx}.tsv && rm {idx}.tmp")
        df = pd.read_table(f"{idx}.tsv", usecols=fields, low_memory=True, sep="\t")
        df.index = pd.MultiIndex.from_arrays(
            df.values.T[
                (0, 1, 3, 4),
            ]
        )
        vcf_dict[idx] = df

    # Originally I averaged TAF, NAF etc. but now I have decided to take the max, and the max allelic depths (eg: 32,0) for the REF and ALT
    # such that they add up and make sense. I'm aware that Strelka does not add up perfectly. This is due to their tier 1 and tier 2 malarky.
    # Strelka does not give a genotype either, which means the genotypes output in thid script are dummies.
    # IIRC, PCGR does not use GT information in their HTML ? double check this.
    avg_df = (
        pd.concat(list(vcf_dict.values()))
        .assign(QUAL=lambda d: d["QUAL"].apply(lambda s: 0 if s == "." else float(s)))
        .groupby(level=[0, 1, 2, 3])
        .agg(
            {
                "TAF": "max",
                "TDP": "max",
                "NAF": "max",
                "NDP": "max",
                "QUAL": "max",
            }
        )
    )
    master_df = pd.concat(list(vcf_dict.values()))

    # I want to keep rs snp ids if they are present
    # I want to take the highest allelic depths found by a somatic caller.
    master_df = master_df.sort_values(["ID", "ADT", "ADN"], ascending=[False, False, False])
    master_df = master_df[~master_df.index.duplicated(keep="first")]

    assert len(avg_df) == len(
        master_df
    ), "averaged values dataframe does not match master dataframe with duplicate index vals removed"

    # add the rest to avg df. revert QUAL back to object
    avg_df = avg_df.assign(QUAL=lambda d: d["QUAL"].apply(lambda s: "." if s == 0.0 else float(s)))
    avg_df["ID"] = master_df["ID"]
    avg_df["FILTER"] = master_df["FILTER"]
    avg_df["ADT"] = master_df["ADT"]
    avg_df["ADN"] = master_df["ADN"]  # dont actually need these in INFO field, theyre in FORMAT and SAMPLE. pop later

    # append TAL from keys
    key_df = key_df[~key_df.index.duplicated(keep="first")]
    avg_df["TAL"] = key_df[4]

    # convert index back to proper names
    avg_df = avg_df.reset_index()
    avg_df.rename(columns={"level_0": "#CHROM", "level_1": "POS", "level_2": "REF", "level_3": "ALT"}, inplace=True)
    fields.append("TAL")
    fields[0] = "#CHROM"
    avg_df = avg_df[fields]

    # these columns are strings for the FORMAT field.
    # dont touch em
    avg_df["FMT_GT"] = "GT"
    avg_df["FMT_DPN"] = "DPN"
    avg_df["FMT_DPT"] = "DPT"
    avg_df["FMT_ADN"] = "ADN"
    avg_df["FMT_ADT"] = "ADT"
    avg_df["FMT_AL"] = "AL"

    # convert TAL to AL
    algo_dict = {
        "1": "mutect2",
        "2": "freebayes",
        "3": "strelka",
        "1,2": "freebayes,mutect2",
        "1,3": "freebayes,strelka",
        "2,3": "mutect2,strelka",
        "1,2,3": "freebayes,mutect2,strelka",
    }

    # sort so the above key:values are always valid
    avg_df["TAL"] = avg_df.TAL.apply(lambda x: ",".join(sorted(x.split(","))))

    # genotype info! objects are ok, int are not for apply() later.
    avg_df["SAMPLE_GT"] = "0/1"
    avg_df["SAMPLE_DPT"] = avg_df["TDP"]  # int
    avg_df["SAMPLE_DPT"] = avg_df["SAMPLE_DPT"].map(str)
    avg_df["SAMPLE_DPN"] = avg_df["NDP"]  # int
    avg_df["SAMPLE_DPN"] = avg_df["SAMPLE_DPN"].map(str)
    avg_df["SAMPLE_ADT"] = avg_df["ADT"]
    avg_df["SAMPLE_ADN"] = avg_df["ADN"]
    # things got weird
    avg_df["TAL"] = avg_df["TAL"].astype("str")
    avg_df["SAMPLE_AL"] = avg_df["TAL"].replace({v: k for k, v in algo_dict.items()})

    # format INFO i.e TAF=0.034;NAF=0;...
    avg_df["TAF"] = avg_df["TAF"].apply(lambda x: "{}{}".format("TAF=", x))
    avg_df["TDP"] = avg_df["TDP"].apply(lambda x: "{}{}".format("TDP=", x))
    avg_df["NAF"] = avg_df["NAF"].apply(lambda x: "{}{}".format("NAF=", x))
    avg_df["NDP"] = avg_df["NDP"].apply(lambda x: "{}{}".format("NDP=", x))
    avg_df["TAL"] = avg_df["TAL"].apply(lambda x: "{}{}".format("TAL=", x))
    # avg_df['ADN'] = avg_df['ADN'].apply(lambda x: "{}{}".format('ADN=', x))
    # avg_df['ADT'] = avg_df['ADT'].apply(lambda x: "{}{}".format('ADT=', x))

    # Add INFO column
    avg_df["INFO"] = avg_df[["NDP", "NAF", "TDP", "TAF", "TAL"]].apply(lambda x: ";".join(x[x.notnull()]), axis=1)
    avg_df = avg_df.drop(
        ["NDP", "NAF", "TDP", "TAF", "TAL", "ADT", "ADN"], axis=1
    )  # drop ADT ADN without using in INFO

    # ADd FORMAT column
    avg_df["FORMAT"] = avg_df[["FMT_GT", "FMT_DPN", "FMT_DPT", "FMT_ADN", "FMT_ADT", "FMT_AL"]].apply(
        lambda x: ":".join(x[x.notnull()]), axis=1
    )
    avg_df = avg_df.drop(["FMT_GT", "FMT_DPN", "FMT_DPT", "FMT_ADN", "FMT_ADT", "FMT_AL"], axis=1)

    # Add sample column
    avg_df[f"{sample_id}"] = avg_df[
        ["SAMPLE_GT", "SAMPLE_DPN", "SAMPLE_DPT", "SAMPLE_ADN", "SAMPLE_ADT", "SAMPLE_AL"]
    ].apply(lambda x: ":".join(x[x.notnull()]), axis=1)
    avg_df = avg_df.drop(["SAMPLE_GT", "SAMPLE_DPN", "SAMPLE_DPT", "SAMPLE_ADN", "SAMPLE_ADT", "SAMPLE_AL"], axis=1)

    avg_df.to_csv("tmp.vcf", sep="\t", index=None, header=True)
    # os.system(f'bcftools view -h {sample_files[0]} | sed \'$d\' > header.txt')
    os.system(f"cat pcgr_header.txt tmp.vcf > {sample_id}.vcf")
    os.system(f"bgzip {sample_id}.vcf")
    os.system(f"tabix {sample_id}.vcf.gz")


def tumor_only(sample):
    # pass meta.id to the script to remove .command* etc files from file search:
    sample_id = sample
    suffixes = (".cns", ".tbi")
    r = re.compile(f"{sample_id}*")
    sample_files = os.listdir("./")
    sample_files = list(filter(r.match, sample_files))
    sample_files = [file for file in sample_files if not file.endswith(suffixes)]

    # there will be a sample_keys.txt file in the work directory
    # load it as key_Df and remove from the VCF list before running vcf2tsvpy
    key = "".join([file for file in sample_files if file.endswith(".txt")])
    key_df = pd.read_table(key, header=None, sep="\t")
    key_df.index = pd.MultiIndex.from_arrays(
        key_df.values.T[
            (0, 1, 2, 3),
        ]
    )
    sample_files.remove(key)

    # stage fields we are interested in capturing
    fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "TDP", "ADT", "TAF"]

    # create dict of dataframes containing dataframes of the abpve fields.
    # strategy is to average TDP,NDP, TAF,NAF and append TAL from $meta.id_keys.txt

    vcf_dict = {}
    for idx, x in enumerate(sample_files):
        idx = idx + 1
        os.system(f"vcf2tsvpy --input_vcf {x} --out_tsv {idx}.tmp --skip_genotype_data")
        os.system(f"tail -n +2 {idx}.tmp > {idx}.tsv && rm {idx}.tmp")
        df = pd.read_table(f"{idx}.tsv", usecols=fields, low_memory=True, sep="\t")
        df.index = pd.MultiIndex.from_arrays(
            df.values.T[
                (0, 1, 3, 4),
            ]
        )
        vcf_dict[idx] = df

    # Originally I averaged TAF, NAF etc. but now I have decided to take the max, and the max allelic depths (eg: 32,0) for the REF and ALT
    # such that they add up and make sense. I'm aware that Strelka does not add up perfectly. This is due to their tier 1 and tier 2 malarky.
    # Strelka does not give a genotype either, which means the genotypes output in thid script are dummies.
    # IIRC, PCGR does not use GT information in their HTML ? double check this.
    avg_df = (
        pd.concat(list(vcf_dict.values()))
        .assign(QUAL=lambda d: d["QUAL"].apply(lambda s: 0 if s == "." else float(s)))
        .groupby(level=[0, 1, 2, 3])
        .agg(
            {
                "TAF": "max",
                "TDP": "max",
                "QUAL": "max",
            }
        )
    )
    master_df = pd.concat(list(vcf_dict.values()))

    # I want to keep rs snp ids if they are present
    # I want to take the highest allelic depths found by a somatic caller.
    master_df = master_df.sort_values(["ID", "ADT"], ascending=[False, False])
    master_df = master_df[~master_df.index.duplicated(keep="first")]

    assert len(avg_df) == len(
        master_df
    ), "averaged values dataframe does not match master dataframe with duplicate index vals removed"

    # add the rest to avg df. revert QUAL back to object
    avg_df = avg_df.assign(QUAL=lambda d: d["QUAL"].apply(lambda s: "." if s == 0.0 else float(s)))
    avg_df["ID"] = master_df["ID"]
    avg_df["FILTER"] = master_df["FILTER"]
    avg_df["ADT"] = master_df["ADT"]

    # append TAL from keys
    key_df = key_df[~key_df.index.duplicated(keep="first")]
    avg_df["TAL"] = key_df[4]

    # convert index back to proper names
    avg_df = avg_df.reset_index()
    avg_df.rename(columns={"level_0": "#CHROM", "level_1": "POS", "level_2": "REF", "level_3": "ALT"}, inplace=True)
    fields.append("TAL")
    fields[0] = "#CHROM"
    avg_df = avg_df[fields]

    # these columns are strings for the FORMAT field.
    # dont touch em
    avg_df["FMT_GT"] = "GT"
    avg_df["FMT_DPT"] = "DPT"
    avg_df["FMT_ADT"] = "ADT"
    avg_df["FMT_AL"] = "AL"

    # convert TAL to AL
    algo_dict = {
        "1": "mutect2",
        "2": "freebayes",
        "3": "strelka",
        "1,2": "freebayes,mutect2",
        "1,3": "freebayes,strelka",
        "2,3": "mutect2,strelka",
        "1,2,3": "freebayes,mutect2,strelka",
    }

    # sort so the above key:values are always valid
    avg_df["TAL"] = avg_df.TAL.apply(lambda x: ",".join(sorted(x.split(","))))

    # genotype info! objects are ok, int are not for apply() later.
    avg_df["SAMPLE_GT"] = "0/1"
    avg_df["SAMPLE_DPT"] = avg_df["TDP"]  # int
    avg_df["SAMPLE_DPT"] = avg_df["SAMPLE_DPT"].map(str)
    avg_df["SAMPLE_ADT"] = avg_df["ADT"]
    # things got weird
    avg_df["TAL"] = avg_df["TAL"].astype("str")
    avg_df["SAMPLE_AL"] = avg_df["TAL"].replace({v: k for k, v in algo_dict.items()})

    # format INFO i.e TAF=0.034;NAF=0;...
    avg_df["TAF"] = avg_df["TAF"].apply(lambda x: "{}{}".format("TAF=", x))
    avg_df["TDP"] = avg_df["TDP"].apply(lambda x: "{}{}".format("TDP=", x))
    avg_df["TAL"] = avg_df["TAL"].apply(lambda x: "{}{}".format("TAL=", x))
    # avg_df['ADN'] = avg_df['ADN'].apply(lambda x: "{}{}".format('ADN=', x))
    # avg_df['ADT'] = avg_df['ADT'].apply(lambda x: "{}{}".format('ADT=', x))

    # Add INFO column
    avg_df["INFO"] = avg_df[["TDP", "TAF", "TAL"]].apply(lambda x: ";".join(x[x.notnull()]), axis=1)
    avg_df = avg_df.drop(
        ["TDP", "TAF", "TAL", "ADT"], axis=1
    )  # drop ADT ADN without using in INFO

    # ADd FORMAT column
    avg_df["FORMAT"] = avg_df[["FMT_GT", "FMT_DPT", "FMT_ADT", "FMT_AL"]].apply(
        lambda x: ":".join(x[x.notnull()]), axis=1
    )
    avg_df = avg_df.drop(["FMT_GT", "FMT_DPT", "FMT_ADT", "FMT_AL"], axis=1)

    # Add sample column
    avg_df[f"{sample_id}"] = avg_df[
        ["SAMPLE_GT", "SAMPLE_DPT", "SAMPLE_ADT", "SAMPLE_AL"]
    ].apply(lambda x: ":".join(x[x.notnull()]), axis=1)
    avg_df = avg_df.drop(["SAMPLE_GT", "SAMPLE_DPT", "SAMPLE_ADT", "SAMPLE_AL"], axis=1)

    avg_df.to_csv("tmp.vcf", sep="\t", index=None, header=True)
    # os.system(f'bcftools view -h {sample_files[0]} | sed \'$d\' > header.txt')
    os.system(f"cat pcgr_header.txt tmp.vcf > {sample_id}.vcf")
    os.system(f"bgzip {sample_id}.vcf")
    os.system(f"tabix {sample_id}.vcf.gz")

def pcgr_ready_vcf(sample):
    sample_id = sample
    suffixes = (".cns", ".tbi")
    r = re.compile(f"{sample_id}*")
    sample_files = os.listdir("./")
    sample_files = list(filter(r.match, sample_files))
    sample_files = [file for file in sample_files if not file.endswith(suffixes)]
    key = "".join([file for file in sample_files if file.endswith(".txt")])
    sample_files.remove(key)

    with VariantFile(sample_files[0]) as fr:
        header = fr.header
        info = list(header.info)
        if 'NAF' in info:
            tumor_normal(sample)
        else:
            tumor_only(sample)

def main():
    # Argument parsing using argparse
    parser = argparse.ArgumentParser(description="Reformat somatic CNA files for PCGR input.")
    parser.add_argument("-s", "--sample", required=True, help="Sample name (meta.id) for the output.")

    args = parser.parse_args()


    # Call reformat_cna function with arguments
    pcgr_ready_vcf(args.sample)


if __name__ == "__main__":
    main()
