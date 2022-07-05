#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CalculateLD_PerGeneWindow_Chunks
# @created     : Monday Jun 20, 2022 21:25:31 CDT
#
# @description : 
######################################################################

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'scripts/')
import CalculateLD_PerGeneWindow
import pandas as pd
import numpy as np
import glob
import re

if hasattr(sys, 'ps1'):
    sys.argv = ["", "hyprcoloc/LociWiseSummaryStatsInput/ForColoc/", "scratch/test_LD_out.", "1" ,"10000"]

LociwiseSummaryStatsFolder, PrefixOut, n, N = sys.argv[1:]
n = int(n)
N = int(N)

f_ins = glob.glob(LociwiseSummaryStatsFolder + "*.txt.gz")
loci = [re.search(LociwiseSummaryStatsFolder + "(.+?).txt.gz", fn).group(1) for fn in f_ins]

df = pd.DataFrame(list(zip(f_ins, loci)), columns=["f_in", "loci"])

df_split = np.array_split(df, N)

for i, row in df_split[n].iterrows():
    print(row['f_in'])
    print(row['loci'])
    CalculateLD_PerGeneWindow.main(row['f_in'], f"{PrefixOut}{row['loci']}.txt.gz")


