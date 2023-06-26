#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Tidy_GTEx_SummaryStats
# @created     : Friday Jun 23, 2023 11:17:38 CDT
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
    sys.argv = ["", "scratch/Tidy.SummaryStats.eQTLsQTLhQTL.tsv.gz" ,"../output/eQTLs_FromSplicingVsChromatin_AcrossGTEx.ForCarlos.tsv.gz"]

_, f_out, f_in = sys.argv

df = pd.read_csv(f_in, sep='\t')
df[['chrom', 'pos', 'ref', 'alt']] = df['TopSNP'].str.split(':', 4, expand=True)
df[['chrom']] = 'chr' + df[['chrom']]
df[['b38']] = 'b38'
df['GtexSnpName'] = df[['chrom', 'pos', 'ref', 'alt', 'b38']].agg('_'.join, axis=1)
df[['pos']] = df[['pos']].astype(int)
 


TabixFilesDict = {re.search("/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/QTLs/(.+?)/NominalPass.cpm.txt.tabix.gz", fn).group(1):pysam.TabixFile(fn, parser=pysam.asTuple()) for fn in glob.glob("/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/QTLs/*/NominalPass.cpm.txt.tabix.gz")} 

out = []
ColumnIndexes = [11,13,14]
for index, row in df.iterrows():
    if index >= 0:
        print(index)
        fetch_region = row["chrom"], row["pos"]-1, row["pos"]-0
        for tissue, tabix_f in TabixFilesDict.items():
            tbx_fetch = tabix_f.fetch(*fetch_region)
            ToReturn = pd.Series([None for i in ColumnIndexes]) 
            for line in tbx_fetch:
                if (row["eGene"] == line[0]) and row['GtexSnpName'] == line[7]:
                    ToReturn = [float(i) for i in itemgetter(*ColumnIndexes)(line)] + [tissue] + row.tolist() 
                    out.append(ToReturn)
                # print(ToReturn)
out_df = pd.DataFrame(out, columns=['nom_pval', 'beta_GTEx', 'betaSE_GTEx', 'tissue'] + df.columns.tolist())

out_df.drop(columns=['chrom', 'pos', 'ref', 'alt', 'b38', 'GtexSnpName']).to_csv(f_out, sep='\t', index=False)

