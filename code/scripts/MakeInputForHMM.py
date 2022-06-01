################################################################
# Perpare count tables for each chromosome, in bins of 200 bp. #
# 5/31/22                                                      #
#  Carlos F Buen Abad Najar, cnajar@uchicago.edu               #
################################################################

import numpy as np
import pandas as pd
import os
import argparse

def MakeCountsPerAssay(counts_dir, assay, chrom, strand):
    assay_dir = counts_dir + '/' + assay + '/'
    assay_samples = [x for x in os.listdir(assay_dir) if x.split('.')[0] == chrom]
    assay_samples = [x for x in assay_samples if x.split('.')[2] == strand]
    
    df = pd.DataFrame()
    
    for assay_sample in assay_samples:
        sample = assay_sample.split('.')[1]
        assay_counts = pd.read_csv(assay_dir + assay_sample, sep='\t',
                    names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'counts'])
        
        counts = assay_counts.counts
        
        if ((assay == 'ProSeq') and (strand == 'minus')):
            counts = -counts
        
        df[assay + '_' + sample] = counts
        
        idx = df.chrom + '_' + df.start.astype(str) + '_' + df.end.astype(str) + '_' + strand
        df.index = idx
        
    return df
        

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--counts_dir', type=str, required=True)
parser.add_argument('--strand', type=str, required=False)
parser.add_argument('--output', type=str, required=False)


if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    counts_dir = args.counts_dir
    strand = args.strand
    output = args.output

    assay_list = os.listdir(counts_dir)
    
    df = pd.DataFrame()
    for assay in assay_list:
        df_assay = MakeCountsPerAssay(counts_dir, assay, chrom, strand)
        df = pd.concat([df, df_assay], axis=1)
        
    df.to_csv(output, sep='\t', index=True, header=True)
        
        
        