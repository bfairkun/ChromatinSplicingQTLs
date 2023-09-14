#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : HandpickedColocExampleScatterPlots.GetLD
# @created     : Thursday Sep 07, 2023 21:10:01 CDT
#
# @description : 
######################################################################

import sys
import pandas as pd
import pysam
import re
import numpy as np
import glob
import vcf as vcf
from operator import itemgetter
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Misc/HandpickedColocPlotsToHighlight.dat.NoLD.tsv.gz" ,"Misc/HandpickedColocPlotsToHighlight.dat.WithLD.tsv.gz"]

_, MyArg1, MyArg2 = sys.argv

def SNP_to_tuple(SNP, sep=":"):
    chrom, pos, ref, alt = SNP.split(sep)
    pos = int(pos)
    return (chrom, pos, ref, alt)

def tupleSNP_to_str(tupleSNP, sep=":"):
    listSNP = list(tupleSNP)
    listSNP[1] = str(tupleSNP[1])
    return(sep.join(listSNP))

def preprocess_genotype(genotype_string):
    """
    Preprocess the genotype string from the VCF file to obtain numeric genotypes.

    Parameters:
        genotype_string (str): Genotype string from the VCF file (e.g., '0|1', '1|0', etc.).

    Returns:
        int: The numeric genotype value (0, 1, or 2).
    """
    alleles = genotype_string.split('|')
    return sum(int(allele) for allele in alleles)

def get_genotypes_from_vcf(vcf_file, chrom, pos1, pos2,name1=None, name2=None):
    """
    Get genotype data for two SNPs from an indexed VCF file using pyVCF.
    Parameters:
        vcf_file (str): Path to the indexed VCF file.
        chrom (str): Chromosome of the SNPs (e.g., 'chr1', 'chrX', etc.).
        pos1 (int): Position of the first SNP.
        pos2 (int): Position of the second SNP.
    Returns:
        tuple: Two numpy arrays containing the genotype data for SNP1 and SNP2.
    """
    try:
        snp1_genotypes = []
        snp2_genotypes = []
        snp1_samples = []
        vcf_reader = vcf.Reader(filename=vcf_file)
        for record in vcf_reader.fetch(chrom, pos1 - 1, pos1):
            if name1 and name1 != record.ID:
                next
            snp1_genotypes = [i.gt_type for i in record.samples]
            snp1_samples = [i.sample for i in record.samples]
        for record in vcf_reader.fetch(chrom, pos2 - 1, pos2):
            if name2 and name2 != record.ID:
                next
            snp2_genotypes = [i.gt_type for i in record.samples]
            snp2_samples = [i.sample for i in record.samples]
        snp1_genotypes = np.array(snp1_genotypes)
        snp2_genotypes = np.array(snp2_genotypes)
        if snp2_samples == snp1_samples:
            return snp1_genotypes, snp2_genotypes
        else:
            return "Samples didn't match"
    except:
        return np.nan


def calculate_correlation(snp1_genotypes, snp2_genotypes):
    """
    Calculate the correlation between two sets of SNP genotypes.

    Parameters:
        snp1_genotypes (numpy.ndarray): Genotype data for SNP1.
        snp2_genotypes (numpy.ndarray): Genotype data for SNP2.

    Returns:
        float: The correlation coefficient between SNP1 and SNP2.
    """
    try:
        return np.corrcoef(snp1_genotypes, snp2_genotypes)[0, 1]
    except:
        return np.nan

def GetR(row, vcf, Snp1ColumnName, Snp2ColumnName):
    SNP1 = SNP_to_tuple(row[Snp1ColumnName])
    SNP2 = SNP_to_tuple(row[Snp2ColumnName])
    gts = get_genotypes_from_vcf(vcf,"chr" +  SNP1[0], SNP1[1], SNP2[1], tupleSNP_to_str(SNP1), tupleSNP_to_str(SNP2))
    try:
        return calculate_correlation(*gts)
    except:
        return np.nan

df = pd.read_csv(MyArg1, sep='\t')

df['R'] = df.apply(GetR, axis=1, vcf="QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz", Snp1ColumnName='snp', Snp2ColumnName='TopSNP')

df.to_csv(MyArg2, sep='\t')

# df.head(150).apply(GetR, axis=1, vcf="QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz", Snp1ColumnName='snp', Snp2ColumnName='TopSNP')

