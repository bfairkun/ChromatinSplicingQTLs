#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : MergeSummaryStatChunks
# @created     : Friday Mar 25, 2022 12:14:28 CDT
#
# @description : 
######################################################################

import sys
import os
import glob
from collections import defaultdict
import shutil
import gzip
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below

if hasattr(sys, 'ps1'):
    sys.argv = ["", "scratch/dirtest_merged/" ,"scratch/dirtest/Dir1/", "scratch/dirtest/Dir2/", "scratch/dirtest/Dir3/"]

DestinationDir = sys.argv[1].rstrip('/')
SourceDirs = [i.rstrip('/') for i in sys.argv[2:]]

# # Make test data
# Dict = {
#         'scratch/dirtest/Dir1':[str(i) for i in range(0,10)],
#         'scratch/dirtest/Dir2':[str(i) for i in range(5, 15)],
#         'scratch/dirtest/Dir3':[str(i) for i in range(10, 20)]
#         }
# for dirname, fn_list in Dict.items():
#     os.makedirs(dirname, exist_ok=True)
#     for fn in fn_list:
#         with open(dirname + '/' + fn + '.txt', 'wt') as fh:
#             _ = fh.write(f"{fn}\n")


#Dict mapping of all filepaths needed to concatenate, with basenames as keys
FilesToMerge = defaultdict(list)
for dirname in SourceDirs:
    for fn in glob.glob(dirname + '/*.txt.gz'):
        FilesToMerge[os.path.basename(fn)].append(fn)
# print(FilesToMerge)


# mkdir for output
os.makedirs(DestinationDir, exist_ok=True)
# Loop thru dictionary, merging files to fdst_fn in DestinationDir
for fdst_fn, fsrc_fnlist in FilesToMerge.items():
    with gzip.open(DestinationDir + "/" + fdst_fn,'wt') as fdst:
        _ = fdst.write("phenotype\tsource_file\tgwas_locus\tsnp\tbeta\tbeta_se\tp\n")
        for fsrc_fn in fsrc_fnlist:
            with gzip.open(fsrc_fn,'rt') as fsrc:
                shutil.copyfileobj(fsrc, fdst)
