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

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "scratch/PairwisePi1Traits.1.txt.gz" ,"scratch/PairwisePi1Traits.P.1.txt.gz"] + "QTLs/QTLTools/chRNA.Expression_cheRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_eRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_lncRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_ncRNA/NominalPassForColoc.txt.tabix.gz QTLs/QTLTools/chRNA.Expression_snoRNA/NominalPassForColoc.txt.tabix.gz".split(' ')

_, f_in, f_out = sys.argv[:3]
tabix_f_in_list = sys.argv[3:]


def GetAscertainmentSNP_P(row, TabixFilesDict):
    """
    from row of df
    """
    fetch_region = row["singletrait_topvar_chr.x"], row["singletrait_topvar_pos.x"]-1, row["singletrait_topvar_pos.x"]
    try:
        tbx_fetch = TabixFilesDict[row["PC2"]].fetch(*fetch_region)
        for line in tbx_fetch:
            if (row["singletrait_topvar.x"] == line[7]) and f'{row["P2"]}:{row["GeneLocus"]}' == line[0]:
                return float(line[11])
    except ValueError:
        return None

def GetAscertainmentSNP_beta(row, TabixFilesDict):
    """
    from row of df
    """
    fetch_region = row["singletrait_topvar_chr.x"], row["singletrait_topvar_pos.x"]-1, row["singletrait_topvar_pos.x"]
    try:
        tbx_fetch = TabixFilesDict[row["PC2"]].fetch(*fetch_region)
        for line in tbx_fetch:
            if (row["singletrait_topvar.x"] == line[7]) and f'{row["P2"]}:{row["GeneLocus"]}' == line[0]:
                return float(line[13])
    except ValueError:
        return None

df = pd.read_csv(f_in, sep='\t')

TabixFilesDict = {re.search("QTLs/QTLTools/(.+?)/NominalPassForColoc.txt.tabix.gz", fn).group(1):pysam.TabixFile(fn, parser=pysam.asTuple()) for fn in tabix_f_in_list}

# dfhead = df.head(100)
# dfhead['trait.x.p.in.y'] = dfhead.apply(GetAscertainmentSNP_P, axis=1, TabixFilesDict=TabixFilesDict )
# dfhead.to_csv(f_out, sep='\t', index=False)



df['trait.x.p.in.y'] = df.apply(GetAscertainmentSNP_P, axis=1, TabixFilesDict=TabixFilesDict )
df['trait.x.beta.in.y'] = df.apply(GetAscertainmentSNP_beta, axis=1, TabixFilesDict=TabixFilesDict )
df.to_csv(f_out, sep='\t', index=False)

