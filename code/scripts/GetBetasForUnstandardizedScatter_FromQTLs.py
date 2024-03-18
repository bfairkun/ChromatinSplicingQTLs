#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : GetBetasForUnstandardizedScatter_FromQTLs
# @created     : Tuesday Mar 05, 2024 12:54:48 CST
#
# @description : 
######################################################################

import sys
import subprocess, shlex
import pandas as pd
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Arg1" ,"Arg2"]

_, MyArg1, MyArg2 = sys.argv

def main():
    pass

if __name__ == '__main__':
    main()

df = pd.read_csv("../output/TraitPairsForUnnormalizedScatter.ForReviewer.tsv.gz",low_memory=False, sep='\t')

len(df['P2'].unique())

df.columns
df[['singletrait_topvar_chr.x', 'singletrait_topvar_pos.x']].drop_duplicates().to_csv("scratch/includePositions.txt", sep=' ', header=False, index=False)

df[['P2']].drop_duplicates().to_csv("scratch/includePhenotypes.txt", sep=' ', header=False, index=False)

df[['singletrait_topvar.x']].drop_duplicates().to_csv("scratch/includeVariants.txt", sep=' ', header=False, index=False)
