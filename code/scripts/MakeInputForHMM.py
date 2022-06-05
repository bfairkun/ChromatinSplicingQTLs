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
    
    df = list()
    
    counter = 1
    
    for assay_sample in assay_samples:
        sample = assay_sample.split('.')[1]
        assay_counts = pd.read_csv(assay_dir + assay_sample, sep='\t',
                    names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'counts'])
        
        counts = assay_counts.counts
        
        if assay == 'ProSeq':
            counts = assay_counts.name
            counts = list(counts.replace('.', '0'))
            counts = np.array([int(x) for x in counts])
            if strand == 'minus':
                counts = -counts
        
        df.append(list(counts))
        chrom = list(assay_counts.chrom)
        start = list(assay_counts.start.astype(str))
        end = list(assay_counts.end.astype(str))
        #df[assay + '_' + sample] = counts
        
        print(sample + ', ' + str(counter) + '/' + str(len(assay_samples)))
        
        counter += 1
        
    print(str(len(df)) + ' samples')
    
    print(str(len(df[0])) + ' bins')

    df = pd.DataFrame(df).T
    print(df.shape)
    df.columns = [assay + '_' + sample.split('.')[1] for sample in assay_samples]
    idx = [chrom[i] + '_' + start[i] + '_' + end[i] + '_' + strand for i in range(len(chrom))]
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
    print('Merging data from assays')
    df = pd.DataFrame()
    for assay in assay_list:
        print('Merging ' + assay + ' data')
        df_assay = MakeCountsPerAssay(counts_dir, assay, chrom, strand)
        df = pd.concat([df, df_assay], axis=1)
        
    df.to_csv(output, sep='\t', index=True, header=True)
        
        
        