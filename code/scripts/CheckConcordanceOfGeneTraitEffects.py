#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Tidy_GTEx_SummaryStats
# @created     : Friday Jun 23, 2023 11:17:38 CDT
#
# @description : 
######################################################################

import sys
import pysam
import pandas as pd
import re
import numpy as np
import glob
import vcf as vcf
from operator import itemgetter
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    # sys.argv = ["", "../output/GwasLoci_with_sQTL_and_hQTL.tsv" ,"../output/GwasLociConcordance.tsv"]
    sys.argv = ["", "../output/GwasLoci_with_sQTL_and_hQTL.tsv" ,"../output/GwasLociConcordanceWithNonEqtls.tsv"]

_, f_in, f_out = sys.argv

df = pd.read_csv(f_in, sep='\t')

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
    snp1_genotypes = []
    snp2_genotypes = []
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

def GetValue(row, TabixFilesDict, SnpColumn, ColumnName):
    """
    from row of df
    """
    SNP = SNP_to_tuple(row[SnpColumn])
    header = TabixFilesDict[row['GWAS.accession']].header[0].split('\t')
    tbx_fetch = TabixFilesDict[row['GWAS.accession']].fetch("chr"+SNP[0], SNP[1], SNP[1]+1, parser=pysam.asTuple())
    try:
        for line in tbx_fetch:
            # print(SNP)
            # print(line)
            gwas_snp = SNP_to_tuple(line[header.index("hm_variant_id")], sep="_")
            if SNP == gwas_snp:

                out = line[header.index(ColumnName)]
                return out
    except:
        return np.nan

def GetR(row, vcf, Snp1ColumnName, Snp2ColumnName):
    SNP1 = SNP_to_tuple(row[Snp1ColumnName])
    SNP2 = SNP_to_tuple(row[Snp2ColumnName])
    gts = get_genotypes_from_vcf(vcf, "chr" + SNP1[0], SNP1[1], SNP2[1], tupleSNP_to_str(SNP1), tupleSNP_to_str(SNP2))
    return calculate_correlation(*gts)


TabixFilesDict = {gwas:pysam.TabixFile(f"gwas_summary_stats/sorted_index_summarystat_hg38beds/{gwas}.bed.gz", parser=pysam.asTuple()) for gwas in df['GWAS.accession'].unique()} 

df['sQTL.Top.GWAS_beta'] = df.apply(GetValue, axis=1, TabixFilesDict=TabixFilesDict, SnpColumn='sQTL.Top.Name', ColumnName="hm_beta")
df['hQTL.Top.GWAS_beta'] = df.apply(GetValue, axis=1, TabixFilesDict=TabixFilesDict, SnpColumn='hQTL.Top.Name', ColumnName="hm_beta")
df['R'] = df.apply(GetR, axis=1, vcf="QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz", Snp1ColumnName='sQTL.Top.Name', Snp2ColumnName='hQTL.Top.Name')

df.to_csv(f_out, sep='\t', index=False)
