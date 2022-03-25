#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CombineAndSplitSummaryStatsForColoc
# @created     : Tuesday Jul 27, 2021 13:58:56 CDT
#
# @description : tab delimited summary stats files contain stats for all genes,
# and in a hypr coloc R script, it might be convenient for debugging and
# parralellization purposes to be able to just read in summary stats for all
# phenotypes for one gene at a time. So here I'll make a script to
# split/combine summary stats files into smaller summary stats files for each
# gwas accession
######################################################################


import sys
import gzip
import pathlib
import os
import pandas as pd

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter, parse_args with hardcoded args below
if hasattr(sys, "ps1"):
    sys.argv = [
        "",
        "scratch/testcolocgwasfiles/",
        "0.01",
        "scratch/test_summarystats/PermutationPassForColoc.txt.gz",
        "scratch/test_summarystats/nominal/eQTL.part.aa",
        "scratch/test_summarystats/PermutationPassForColoc.txt.gz",
        "scratch/test_summarystats/nominal/eQTL.part.ab",
        "scratch/test_summarystats/PermutationPassForColoc.txt.gz",
        "scratch/test_summarystats/nominal/eQTL.part.ac",
        "scratch/test_summarystats/PermutationPassForColoc.txt.gz",
        "scratch/test_summarystats/nominal/eQTL.part.ad",
    ]

FilesOutPrefix = sys.argv[1]
NominalPThreshold = float(sys.argv[2])
SummaryStatFilesIn = sys.argv[3:]
PairedSummaryStatFilesIn = tuple(
    zip(SummaryStatFilesIn[0::2], SummaryStatFilesIn[1::2])
)

def return_file_handle(input_file, open_mode="rt"):
    """
    Handles compressed and uncompressed files. Accepts open modes r/w/w+
    """
    # compressed
    if input_file.endswith(".gz"):
        with gzip.open(input_file, open_mode) as gzipped_file_handle:
            yield gzipped_file_handle
    else:
        with open(input_file, open_mode) as normal_fh:
            yield normal_fh


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


FilesOutPrefixDir = os.path.dirname(FilesOutPrefix)
if FilesOutPrefixDir:
    pathlib.Path(FilesOutPrefixDir).mkdir(parents=True, exist_ok=True)

FilesOpened = set([])

for PermutationPassFileIn, NominalPassFileIn in PairedSummaryStatFilesIn:
    permutation_df = pd.read_csv(PermutationPassFileIn, sep=" ")
    PhenotypesToInclude = set(
        permutation_df.loc[permutation_df["adj_emp_pval"] < NominalPThreshold]["phe_id"]
    )
    OpenFileOut = None
    if is_gz_file(NominalPassFileIn):
        fh = gzip.open(NominalPassFileIn, 'rt')
    else:
        fh = open(NominalPassFileIn, 'rt')
    for i, line in enumerate(fh):
        if i >= 0:
            l = line.strip("\n").split(" ")
            if l[0] in PhenotypesToInclude:
                phenotype, gwas_locus = l[0].rsplit(":",1)
                gwas_accession = gwas_locus.split("_")[-1]
                FileOut = FilesOutPrefix + gwas_accession + ".txt.gz"
                if FileOut != OpenFileOut:
                    try:
                        fh_out.close()
                    except:
                        pass
                    if FileOut not in FilesOpened:
                        with gzip.open(FileOut, "wt") as fh_out:
                            _ = fh_out.write("")
                        FilesOpened.add(FileOut)
                    fh_out = gzip.open(FileOut, "at")
                _ = fh_out.write(
                    #phenotype\tsource_file\tgwas_locus\tsnp\tbeta\tbeta_se\tp\n
                    f"{phenotype}\t{NominalPassFileIn}\t{gwas_locus}\t{l[7]}\t{l[13]}\t{l[14]}\t{l[11]}\n"
                )
                OpenFileOut = FileOut
    fh.close()
    fh_out.close()
