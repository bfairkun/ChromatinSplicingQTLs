#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : cnajar (cnajar@midway2-login1.rcc.local)
# @file        : ScorePWM
# @created     : Tuesday Feb 03, 2023 23:46:03 CST
#
# @description : Input tab delimited file (output of bedtools getfasta -tab
# -name) and output the same file with a column added for the MaxEntScan score
# for each sequence.
######################################################################

import sys
from maxentpy import maxent 
import pandas as pd

_, MyArg1, MyArg2 = sys.argv

if 'FivePrime' in MyArg1:
    score_ = maxent.score5
elif 'ThreePrime' in MyArg1:
    score_ = maxent.score3
else:
    raise Exception('Tab file must contain FivePrime or ThreePrime on its name.')

df_ss = pd.read_csv(MyArg1, sep='\t', names = ['name', 'sequence'])


maxent_score = []
for idx, row in df_ss.iterrows():
    maxent_score.append(score_(row.sequence))

df_ss['MaxEntScan'] = maxent_score

df_ss.to_csv(MyArg2, sep='\t', header=False, index=False)


