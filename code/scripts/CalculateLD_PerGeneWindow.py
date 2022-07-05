#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : CalculateLD_PerGeneWindow
# @created     : Thursday Jun 16, 2022 13:59:31 CDT
#
# @description : 
######################################################################

import sys
import pandas as pd
import pysam
import numpy as np
import os
from operator import itemgetter


def matrix_to_xy(df, columns=None, reset_index=False):
    bool_index = np.triu(np.ones(df.shape)).astype(bool)
    xy = (
        df.where(bool_index).stack().reset_index()
        if reset_index
        else df.where(bool_index).stack()
    )
    if reset_index:
        xy.columns = columns or ["row", "col", "val"]
    return xy

def main(GenewiseSummaryStats_f_in, MeltedLD_mat_f_out):
    """
    GenewiseSummaryStats_f_in: input genewise summary stats filepath.
    MeltedLD_mat_f_out: name output filepath of melted LD matrix
    """
    snps = pd.read_csv(GenewiseSummaryStats_f_in, sep='\t')['snp'].unique()
    min_snp = min([int(snp.split(':')[1]) for snp in snps])
    max_snp = max([int(snp.split(':')[1]) for snp in snps])
    chrom = snps[0].split(':')[0]
    vcf_f = f'Genotypes/1KG_GRCh38/{chrom}.vcf.gz'
    bcf_in = pysam.VariantFile(vcf_f)


    samples_in_project = pd.read_csv("config/samples.tsv", sep='\t', comment='#')['IndID'].unique()
    ThousandGenomesSamplesMetadata = pd.read_csv("../data/igsr_samples.tsv.gz", sep='\t')
    YRI_list = ThousandGenomesSamplesMetadata.loc[ThousandGenomesSamplesMetadata['Population code']=="YRI"]['Sample name'].tolist()

    # get list of samples to get genotypes for and calculate LD mat
    SampleList = [s for s in samples_in_project if (s in YRI_list) and (s in list(bcf_in.header.samples) ) ]

    GT_list = []
    snp_names_list = []
    for i, record in enumerate(bcf_in.fetch(f'chr{chrom}', min_snp, max_snp)):
        if i >=0:
            if record.id in snps:
                # print(record.id)
                GTs = [sum(s['GT']) for s in itemgetter(*SampleList)(record.samples)]
                GT_list.append(GTs)
                snp_names_list.append(record.id)
    GT_array = np.array(GT_list)
    # np.corrcoef might produce RuntimeWarning: https://stackoverflow.com/questions/45897003/python-numpy-corrcoef-runtimewarning-invalid-value-encountered-in-true-divide
    df = pd.DataFrame(np.corrcoef(GT_array), columns=snp_names_list , index=snp_names_list)
    keep = np.triu(np.ones(df.shape)).astype('bool').reshape(df.size)
    matrix_to_xy(df).rename_axis(['snpA', 'snpB']).to_csv(MeltedLD_mat_f_out, sep='\t', float_format='%.3f', header=['R'])

if __name__ == "__main__":
    # I like to script and debug with an interactive interpreter in the same
    # working dir as this script.. If using interactive interpreter, can step
    # thru the script with args defined below
    try:
        if sys.ps1:
            args = ["hyprcoloc/LociWiseSummaryStatsInput/ForColoc/ENSG00000000419.12.txt.gz", "scratch/test_LDMat.txt.gz"]
    except:
        args = sys.argv[1:3]
    main(*args)

