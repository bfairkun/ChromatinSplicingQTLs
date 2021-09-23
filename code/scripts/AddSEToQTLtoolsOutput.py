#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
######################################################################
# @author      : benfair (benfair@Bens-MacBook-Air.local)
# @file        : cli
# @created     : Saturday Dec 26, 2020 16:47:20 EST
#
# @description : add SE to QTLtools nominal pass output from P
######################################################################
"""
This script adds a column to QTLtools output for the SE of the beta, based on the degrees of freedom (determined by sample size and number of covariates) and P value
"""

import sys
import os
import argparse
import vcf
from scipy.stats import t
import gzip

def cmdline_args(Args=None):
    p = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("QTLsIn",
            metavar="<QTLtoolsNominalOutput.gz>")
    p.add_argument("Covariates",
            metavar="<space delimited QTLtools covariates file>")
    p.add_argument("vcf",
            metavar="<vcf.gz>"),
    p.add_argument("FileOut",
            metavar="<FileOut.gz>")
    return p.parse_args(Args)

def main(some_args):
    pass

if __name__ == '__main__':
    # I like to script and debug with an interactive interpreter. If using
    # interactive interpreter, parse_args with hardcoded args below
    if hasattr(sys, 'ps1'):
        Args = ["QTLs/QTLTools/H3K4ME3/NominalPass_ForColoc.txt.gz", "QTLs/QTLTools/H3K4ME3/OnlyFirstReps.sorted.qqnorm.bed.pca", "QTLs/QTLTools/H3K4ME3/Genotypes/WholeGenome.vcf.gz", "scratch/FileOut.test.gz"]
        args = cmdline_args(Args=Args)
    else:
        args = cmdline_args()
    print(vars(args))
    main(args)

CovariateCount=0
with open(args.Covariates, "rt") as f:
    for i, line in enumerate(f):
        if i == 0:
            covariates_samples = set(line.strip().split(' ')[1:])
        else:
            CovariateCount += 1
print(CovariateCount)

vcf_in = vcf.Reader(filename=args.vcf)
vcf_samples = set(vcf_in.samples)

len(vcf_samples)
len(covariates_samples)
N = len(vcf_samples & covariates_samples)
print("{} samples intersect both vcf and covariates file".format(N))
df = N - CovariateCount - 1

with gzip.open(args.QTLsIn, "rt") as f_in, gzip.open(args.FileOut, "wt") as f_out:
    # f_out.write("phenotype\tsnp\tbeta\tbeta_se\tp\n")
    for i,line in enumerate(f_in):
        if i >= 0:
            l = line.split(' ')
            phenotype, snp, p, beta = [l[0], l[7], l[11], l[12]]
            try:
                tval = t.ppf(float(p)/2, df=df)
                beta_se = abs(float(beta)/tval)
                f_out.write("{}\t{}\t{}\t{:.6}\t{}\n".format(phenotype, snp, beta, beta_se, p))
            except ValueError:
                print(i)

