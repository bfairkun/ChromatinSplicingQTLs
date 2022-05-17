#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : GatherTopSNPPvalsCrossTraits
# @created     : Friday May 13, 2022 11:38:21 CDT
#
# @description : 
######################################################################

import sys
import pysam
import os
import pandas as pd
import statsmodels.stats.multitest as multi
import gzip
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "scratch/TestPvalsForPi.txt.gz", "0.1" ,"MetabolicLabelled.30min", "MetabolicLabelled.60min", "Expression.Splicing.Subset_YRI"]

f_outgz = sys.argv[1]
fdr = float(sys.argv[2])
Phenotypes = sys.argv[3:]

PermutationPassFns = [(p, f"QTLs/QTLTools/{p}/PermutationPassForColoc.txt.gz") for p in Phenotypes ]
NominalPassTbxs = [(p, pysam.TabixFile(f"QTLs/QTLTools/{p}/NominalPassForColoc.txt.tabix.gz")) for p in Phenotypes ]


# Read in permutation pass results to df and calculate fdr adjusted p
li = []
for p, fn in PermutationPassFns:
    df = pd.read_csv(fn, sep=' ')
    df['phenotype_class'] = p
    df['q'] = multi.fdrcorrection(df['adj_beta_pval'])[1]
    li.append(df)
df = pd.concat(li, axis=0, ignore_index=True)

# Get TopSNPs for significant QTLs
TopSnps = df.loc[df['q']<fdr]

TopSnps.groupby('phenotype_class').size()

# Write out Pvals for those top SNPs across all phenotype classes
with gzip.open(f_outgz, mode='wt') as f:
    _ = f.write("PhenotypeID\tDiscoveryPhenotypeClass\tDiscoveryPhenotypePermutationP\tDiscoveryPhenotypePermutationQ\tTestPhenotypeClass\tTestPhenotypeP\n")
    for row in TopSnps.itertuples():
        for p, tbx in NominalPassTbxs:
            try:
                for line in tbx.fetch(row.var_chr, row.var_from-1, row.var_to+1):
                    l = line.strip('\n').split('\t')
                    TestP = l[11]
                    TestPhenotypeID = l[0]
                    TestVarID = l[7]
                    if (TestPhenotypeID == row.phe_id) and (TestVarID == row.var_id):
                        _ = f.write(f"{row.phe_id}\t{row.phenotype_class}\t{row.adj_beta_pval}\t{row.q}\t{p}\t{TestP}\n")
            except ValueError: pass


