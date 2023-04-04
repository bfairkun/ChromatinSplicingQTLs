#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CalculatePi1_GetAscertainmentP
# @created     : Thursday Jun 23, 2022 12:15:33 CDT
#
# @description : 
######################################################################

import sys
import pysam
import pandas as pd
import re
import numpy as np
import glob
from operator import itemgetter

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    # sys.argv = ["", "scratch/PairwisePi1Traits.1.txt.gz" ,"scratch/PairwisePi1Traits.P.1.txt.gz"] + "QTLs/QTLTools/chRNA.Expression_cheRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_eRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_lncRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_ncRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_snoRNA/NominalPassForColoc.txt.tabix.gz".split(' ')
    sys.argv = "scripts/CalculatePi1_GetAscertainmentP_AllPairs.py pi1/PairwiseTraitsToCompare/32.txt.gz pi1/PairwiseTraitsToCompare/P.32.txt.gz QTLs/QTLTools/Expression.Splicing/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/Expression.Splicing.Subset_YRI/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression.Splicing/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/MetabolicLabelled.30min/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/MetabolicLabelled.60min/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/H3K27AC/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/H3K4ME3/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/H3K4ME1/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/H3K36ME3/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/APA_Nuclear/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/APA_Total/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/polyA.Splicing/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/polyA.Splicing.Subset_YRI/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Splicing/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/ProCap/NominalPassForColoc.txt.tabix.gz".split(' ')
    # sys.argv = ["", "scratch/PairwisePi1Traits.1.txt.gz" ,"scratch/PairwisePi1Traits.P.1.txt.gz"] + glob.glob("QTLs/QTLTools/*/NominalPassForColoc.txt.tabix.gz")

_, f_in, f_out = sys.argv[:3]
tabix_f_in_list = sys.argv[3:]

def GetSummaryStatsFromOtherTabixFiles(row, TabixFilesDict, ColumnIndexes):
    """
    from row of df
    """
    fetch_region = row["singletrait_topvar_chr.x"], row["singletrait_topvar_pos.x"]-1, row["singletrait_topvar_pos.x"]
    try:
        tbx_fetch = TabixFilesDict[row["PC2"]].fetch(*fetch_region)
        ToReturn = pd.Series([None for i in ColumnIndexes]) 
        for line in tbx_fetch:
            if (row["singletrait_topvar.x"] == line[7]) and f'{row["P2"]}:{row["GeneLocus"]}' == line[0]:
                ToReturn = pd.Series([float(i) for i in itemgetter(*ColumnIndexes)(line)])
        return ToReturn
    except ValueError:
        return pd.Series([None for i in ColumnIndexes])
# GetSummaryStatsFromOtherTabixFiles(dfhead.loc[1], TabixFilesDict, [11,13,14])

df = pd.read_csv(f_in, sep='\t')

TabixFilesDict = {re.search("QTLs/QTLTools/(.+?)/NominalPassForColoc.txt.tabix.gz", fn).group(1):pysam.TabixFile(fn, parser=pysam.asTuple()) for fn in tabix_f_in_list}

# dfhead = df.head(100)
# dfhead[['trait.x.p.in.y', 'x.beta.in.y', 'x.beta_se.in.y']] = dfhead.apply(GetSummaryStatsFromOtherTabixFiles, axis=1, TabixFilesDict=TabixFilesDict, ColumnIndexes = [11, 13,14])
# dfhead['trait.x.p.in.y'] = dfhead.apply(GetAscertainmentSNP_P, axis=1, TabixFilesDict=TabixFilesDict )
# dfhead.to_csv(f_out, sep='\t', index=False)

df[['trait.x.p.in.y', 'x.beta.in.y', 'x.beta_se.in.y']] = df.apply(GetSummaryStatsFromOtherTabixFiles, axis=1, TabixFilesDict=TabixFilesDict, ColumnIndexes = [11, 13,14])

df.to_csv(f_out, sep='\t', index=False)

