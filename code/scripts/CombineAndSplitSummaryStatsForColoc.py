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
# gene
######################################################################


import sys
import gzip
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter, parse_args with hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "scratch/testcolocfiles/" ,"scratch/FileOut.test.gz"]

FilesOutPrefix = sys.argv[1]
SummaryStatFilesIn = sys.argv[2:]
print(SummaryStatFilesIn)

# FilesOpened = set([])

# OpenFile = None
# for FileIn in SummaryStatFilesIn:
#     with gzip.open(FileIn, 'rt') as fh:
#         for i, line in enumerate(fh):
#             if i>= 1:
#                 l = line.strip('\n').split('\t')
#                 gene = l[0].split(':')[1]
#                 FileOut = FilesOutPrefix + gene + ".txt.gz"
#                 if FileOut not in FilesOpened:
#                     with gzip.open(FileOut, 'wt') as fh_out:
#                         _ = fh_out.write("phenotype\tsnp\tbeta\tbeta_se\tp\n")
#                     FilesOpened.add(FileOut)
#                 if FileOut != OpenFile:
#                     try: fh_out.close()
#                     except: pass
#                     fh_out = gzip.open(FileOut, 'at')
#                 _ = fh_out.write(line)
#                 OpenFile = FileOut
#         fh_out.close()


