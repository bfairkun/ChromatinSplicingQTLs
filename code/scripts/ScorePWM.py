#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : ScorePWM
# @created     : Tuesday Jan 25, 2022 11:14:03 CST
#
# @description : Input tab delimited file (output of bedtools getfasta -tab
# -name) and output the same file with a column added for the score based on a
# pwm from all sequences.
######################################################################

import sys
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab" ,"SplicingAnalysis/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz"]

_, MyArg1, MyArg2, MyArg3 = sys.argv

df = pd.read_csv(MyArg1, delimiter="\t", names=['name', 'seq'])
df['seq'] = df['seq'].apply(Seq)

m  = motifs.create(df['seq'])
m.weblogo(MyArg3)
pwm = m.counts.normalize()
pssm = pwm.log_odds()

# pssm.calculate(m.consensus)
df['score'] = df['seq'].apply(pssm.calculate)

df.to_csv(MyArg2, sep='\t', header=False, index=False)


