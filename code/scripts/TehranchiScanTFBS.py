#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TehranchiScanTFBS
# @created     : Thursday Apr 14, 2022 15:36:47 CDT
#
# @description : 
######################################################################

import sys
import pysam
import numpy as np
import gzip
from Bio import motifs
from Bio.Seq import Seq
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "../data/20220414203249_JASPAR2022_combined_matrices_25818_jaspar.txt" ,"ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa", "Tehranchi/TableS2.hg38.bed", "scratch/ScanTFBS.2.txt.gz"]

_, jaspar_TF_motifs_f, RefFasta, TehranchiBed, f_out = sys.argv

motifs_to_search = list()
jaspar_motifs_fh = open(jaspar_TF_motifs_f)
for m in motifs.parse(jaspar_motifs_fh, "jaspar"):
    motifs_to_search.append(m)
NumFlankingBasesToSearch = max([len(m) for m in motifs_to_search])

for m in motifs_to_search:
    print(m.alphabet)
    print(m)
    # print(str(m).partition('\n')[0])

pssms = [m.counts.normalize(pseudocounts=0.5).log_odds() for m in motifs_to_search]
pssm_names = [str(m).partition('\n')[0][8:] for m in motifs_to_search]
pssm_thresholds = [pssm.distribution().threshold_patser() for pssm in pssms]
# pssm_thresholds = [pssm.distribution().threshold_fpr(0.0001) for pssm in pssms]

pssms_dict = dict(zip(pssm_names, zip(pssms, pssm_thresholds)))

for motif, (pssm, threshold) in pssms_dict.items():
    print (motif, pssm.degenerate_consensus, pssm.max, threshold)

genome = pysam.FastaFile(RefFasta)

genome.fetch()

with open(TehranchiBed, 'r') as f, gzip.open(f_out, 'wt') as fh_out:
    _ = fh_out.write(f'chrom\tstart\tstop\tTF\tP\tmotif\tDelta\tMaxAnyValue\thg38refIsHighbind\n')
    for i, line in enumerate(f):
        if i <= 500:
        # if i >= 0:
            chrom, start, stop, TF, hg19alt, hg19ref, highbind, P = line.strip('\n').split('\t')
            start = int(start) -1
            stop = int(stop) -1
            try:
                seq = Seq(genome.fetch(chrom, start-NumFlankingBasesToSearch, stop+NumFlankingBasesToSearch))
            except KeyError:
                continue
            hg38ref = seq[NumFlankingBasesToSearch]
            if str(hg38ref) == hg19ref:
                hg38alt = hg19alt
            elif str(hg38ref) == hg19alt:
                hg38alt = hg19ref
            else:
                continue
            RefSearchStr = Seq(seq[0:NumFlankingBasesToSearch] + hg38ref + seq[NumFlankingBasesToSearch+1:])
            AltSearchStr = Seq(seq[0:NumFlankingBasesToSearch] + hg38alt + seq[NumFlankingBasesToSearch+1:])
            # print(TF, hg19ref, hg19alt, seq, hg38ref, RefSearchStr, AltSearchStr)
            for motif, (pssm, threshold) in pssms_dict.items():
                RefFwd = pssm.calculate(RefSearchStr)
                RefRev = pssm.reverse_complement().calculate(RefSearchStr)
                AltFwd = pssm.calculate(AltSearchStr)
                AltRev = pssm.reverse_complement().calculate(AltSearchStr)
                RefScores = np.array([RefFwd, RefRev])
                AltScores = np.array([AltFwd, AltRev])
                MaxRefValue = np.amax(RefScores)
                MaxAltValue = np.amax(AltScores)
                MaxAnyValue = max([MaxRefValue, MaxAltValue])
                if MaxAnyValue >= threshold:
                    # print("motif  hit")
                    if MaxRefValue > MaxAltValue:
                        # print(AltScores[np.where(RefScores == MaxRefValue)])
                        Delta = MaxRefValue - AltScores[np.where(RefScores == MaxRefValue)]
                    elif MaxRefValue < MaxAltValue:
                        # print(RefScores[np.where(AltScores == MaxAltValue)])
                        Delta = RefScores[np.where(AltScores == MaxAltValue)] - MaxAltValue
                    elif MaxRefValue == MaxAltValue:
                        Delta = [0]
                    else:
                        print("Error")
                        raise
                    _ = fh_out.write(f'{chrom}\t{start}\t{stop}\t{TF}\t{P}\t{motif}\t{Delta[0]}\t{MaxAnyValue}\t{str(hg38ref)==highbind}\n')
