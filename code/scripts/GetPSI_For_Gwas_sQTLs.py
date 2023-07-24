#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : GetPSI_For_Gwas_sQTLs
# @created     : Tuesday Jul 18, 2023 05:47:49 CDT
#
# @description : 
######################################################################

import sys
import pandas as pd
import os
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "../output/sQTLsThatColocWithGWAS.tsv.gz" ,"SplicingAnalysis/NormalizedPsiTables/CountsOverClusterMax_Expression.Splicing.bed.gz", "Genotypes/1KG_GRCh38/Autosomes.vcf.gz", "scratch/test.PSI.tsv.gz"]

_, InputTable_f, Bedgz_f, VCF_f, Output_f = sys.argv

DF = pd.read_csv(InputTable_f, sep="\t")

DF.columns

if os.path.isfile(Output_f):
    os.remove(Output_f)

for index, row in DF.iterrows():
    # print(row['Intron'], row['TopCandidateSNP'])
    print(index)
    if index>=0:
        FeatureName = row['Trait'].split(';')[1]
        SnpPos = ':'.join(row['TopCandidateSNP'].split(':')[0:2])
        cmd = f"python scripts/GenometracksByGenotype/ExtractPhenotypeBedByGenotype.py --Bed {Bedgz_f} --VCF {VCF_f} -vvv --SnpName {row['TopCandidateSNP']} --SnpPos chr{SnpPos} --Long  --NoHeader --FeatureName {FeatureName} | gzip - >> {Output_f}"
        os.system(cmd)



#python ExtractPhenotypeBedByGenotype.py --Bed /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/QTLs/QTLTools/Expression.Splicing/OnlyFirstRepsUnstandardizedForColoc.sorted.qqnorm.bed.gz --VCF /project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38/Autosomes.vcf.gz -vvv --SnpName 14:105172594:T:C --SnpPos chr14:105172594 --Long >> ../../scratch/test.tsv
