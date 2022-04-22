#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
"""
GetGWASLeadVariantWindows

main function returns a dataframe and optionally writes it out of lead snp rows
of summary stats for harmonized summary stats from gwas catalog as described in
Pheonix's paper

usage: python GetGWASLeadVariantWindows.py <Input.h.tsv.gz> [<FileOut>
        <PvalueThreshold>]
"""

import sys
import pandas as pd

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, "ps1"):
    sys.argv = [
        "",
        # "/scratch/midway2/bjf79/27618452-GCST006258-EFO_0006336.h.tsv.gz",
        "gwas_summary_stats/NonGwasCatalogBeds/IMSGC2019.bed.gz",
    ]


df = pd.read_csv(sys.argv[1], sep='\t', usecols=range(0,4), header=0, names=["chrom", "start", "stop", "P"])

def main(
    GWASCatalogHarmonizedSummaryStats_f,
    FileOut=None,
    threshold=5e-8,
):
    """main script function used to read in a file downloaded from GWAS
    catalog of harmonized summary stats, and return a pandas data frame of the
    lead SNPs for each genomewide significant locus"""
    # Big files sometimes throw errors when guessing column types that aren't easy
    # to guess by looking at first bunch of rows. Solution: Read first 1000 rows to
    # guess column types, then manually reassign the ones that I know might be
    # finnicky and read_csv with all rows
    df = pd.read_csv(sys.argv[1], sep='\t', usecols=range(0,4), header=0, names=["hm_chrom", "hm_pos", "stop", "p_value"])

    # Filter rows for significant P-value and no problematic NAs
    df = df.loc[
        (df["p_value"] < float(threshold))
        & pd.notna(df["hm_chrom"])
        & (pd.notna(df["p_value"]))
    ]

    LeadSNPList = []
    Iterations = 1

    while not df.empty:
        # print(Iterations)
        df.sort_values("p_value", inplace=True)
        LeadSNPList.append(df.iloc[0])
        LeadSNPChr = df.iloc[0]["hm_chrom"]
        LeadSNPPos = df.iloc[0]["hm_pos"]
        df = df.loc[
            ~(
                (df["hm_chrom"] == LeadSNPChr)
                & (
                    df["hm_pos"].between(
                        LeadSNPPos - 5e5, LeadSNPPos + 5e5
                    )
                )
            )
        ]
        Iterations += 1

    LeadSnps = pd.DataFrame(LeadSNPList)
    if FileOut:
        LeadSnps.to_csv(FileOut, sep="\t", index=False)
    return LeadSnps


if __name__ == "__main__":
    main(*sys.argv[1:])
