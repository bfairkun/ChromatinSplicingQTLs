#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : SmallMoleculeParseCassetteExons
# @created     : Wednesday Mar 15, 2023 12:31:50 CDT
#
# @description : 
######################################################################

import sys
from Bio.Seq import Seq
import pysam
import pandas as pd
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "SmallMolecule/CassetteExons/FlankingBases.tsv", "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz"]

_, FlankBases_fn, gtf_fn = sys.argv

flank_df = pd.read_csv(FlankBases_fn, sep='\t')

gtf = pysam.TabixFile(gtf_fn)
for index, row in flank_df.iterrows():
    # print(row['DownstreamFlankBaseEnd'], row['SkipJuncName'])
    entries = set()
    for entry in gtf.fetch(row['Chrom'], row['UpstreamFlankBaseStart'], row['UpstreamFlankBaseEnd'], parser=pysam.asGTF()):
        if entry.feature == 'CDS':
            entries.add((entry.contig, entry.feature, entry.start, entry.end, entry.frame))
    if len(entries) == 2:
        print(row['SkipJuncName'], row['UpstreamFlankBaseStart'])
        print(entries)
# print(entries)



def main():
    pass

if __name__ == '__main__':
    main()

