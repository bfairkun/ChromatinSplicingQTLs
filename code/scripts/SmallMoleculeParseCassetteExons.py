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
    sys.argv = "scripts/SmallMoleculeParseCassetteExons.py SmallMolecule/CassetteExons/ExonsToTranslate.Translated.tsv.gz SmallMolecule/CassetteExons/ExonsToTranslate.tsv.gz ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa".split(' ')
    # sys.argv = ["", "scratch/test.out.tsv.gz","SmallMolecule/CassetteExons/ExonsToTranslate.tsv.gz", "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz", "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"]

_, f_out, FlankBases_fn, gtf_fn, fa_fn = sys.argv

flank_df = pd.read_csv(FlankBases_fn, sep='\t')

ref = pysam.FastaFile(fa_fn)

# Define a function to get the sequence of an exon from the reference genome
def get_exon_sequence(chrom, start, end, strand):
    exon_seq = Seq(ref.fetch(chrom, start, end))
    if strand == "-":
        exon_seq = Seq(exon_seq).reverse_complement()
    return exon_seq

# flank_df.loc[1]
# gtf = pysam.TabixFile(gtf_fn)
# for entry in gtf.fetch("chr1", 8422726, 8422807, parser=pysam.asGTF()):
#     if entry.feature == "CDS" and "ENST00000377464.5" in entry.attributes:
#         # ref.fetch(entry.contig, entry.start, entry.end)
#         get_exon_sequence(entry.contig, entry.start, entry.end, entry.strand)
#         # print(entry)

gtf = pysam.TabixFile(gtf_fn)
for index, row in flank_df.iterrows():
    if index>=0:
        # Initialize some variables
        us_cds = None
        cassette_cds = None
        ds_cds = None
        # Get upstream exon cds and frame
        for entry in gtf.fetch(row['chrom'], row['ExonStart_Upstream'], row['ExonStop_Upstream'], parser=pysam.asGTF()):
            if entry.feature == "CDS" and row['transcript'] in entry.attributes:
                us_cds = get_exon_sequence(entry.contig, entry.start, entry.end, entry.strand)
                us_cds_frame = entry.frame
                break
        # if the upstream exon is coding, get cassette exon and downstream exon cds
        if us_cds:
            cassette_cds = get_exon_sequence(row['chrom'], row['ExonStart_Cassette'], row['ExonStop_Cassette'], row['Strand'])
            ds_cds = get_exon_sequence(row['chrom'], row['ExonStart_Downstream'], row['ExonStop_Downstream'], row['Strand'])
            Skipped = (us_cds + ds_cds)[us_cds_frame:].translate().split('*')[0]
            Included = (us_cds + cassette_cds + ds_cds)[us_cds_frame:].translate().split('*')[0]
            # print(len(Included) - len(Skipped))
            flank_df.loc[index,'IncludedTranslation'] = str(Included)
            flank_df.loc[index,'SkippedTranslation'] = str(Skipped)

flank_df.to_csv(f_out, sep='\t', index=False)


