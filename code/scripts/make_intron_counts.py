import numpy as np
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--phenotype', type=str, required=True)
parser.add_argument('--output', type=str, required=True)


if __name__ == '__main__':
    
    args = parser.parse_args()
    phenotype = args.phenotype
    output = args.output
    
    samples = [x.split('.')[0] for x in os.listdir('SplicingAnalysis/SPLICEq/{phenotype}.IER/'.format(phenotype=phenotype))]
    
    if phenotype == 'chRNA.IER':
        samples = [x for x in samples if x != 'NA18855']
    
    ix_list = []
    for sample in samples:
        print(sample)
        intron_df = pd.read_csv('SplicingAnalysis/SPLICEq/{phenotype}.IER/{sample}.spliceq.tab'.format(phenotype=phenotype, 
                                                                                                   sample=sample), 
                                sep='\t')
        intron_df['length'] = intron_df.IEnd - intron_df.IStart

        ix = intron_df.groupby(['chr', 'IStart', 'IEnd', 'strand', 'gene_ID'])[['length', 'intron_cov']].first()
        ix.columns = [sample + '_length', sample + '_intron_cov']

        ix_list.append(ix)

    X = pd.concat(ix_list, axis=1)
    
    X_length = X[[i for i in X.columns if i.split('_')[1] == 'length']]
    X_intron = X[[i for i in X.columns if i.split('_')[1] == 'intron']]

    X_intron = X_intron.fillna(0)
    X_length = X_length.fillna(0)
    X_intron.columns = [x.split('_')[0] for x in X_intron]
    
    IntronLength = X_length.max(axis=1).groupby('gene_ID').sum()
    IntronCounts = X_intron.groupby('gene_ID').sum()
    
    IntronCounts['length'] = IntronLength
    
    
    bed_idx = pd.DataFrame()
    bed_idx['chrom'] = ['chr' + str(x[0]) for x in X_intron.index]
    bed_idx['start'] = [x[1] for x in X_intron.index]
    bed_idx['end'] = [x[2] for x in X_intron.index]
    bed_idx['strand'] = [x[3] for x in X_intron.index]
    bed_idx['gene_ID'] = [x[4] for x in X_intron.index]

    start = bed_idx.groupby('gene_ID').start.min()
    end = bed_idx.groupby('gene_ID').end.max()
    chrom = bed_idx.groupby('gene_ID').chrom.first()
    strand = bed_idx.groupby('gene_ID').strand.first()

    IntronCounts['Chr'] = chrom
    IntronCounts['start'] = start
    IntronCounts['end'] = end
    IntronCounts['pid'] = IntronCounts.index
    IntronCounts['gid'] = IntronCounts.index
    IntronCounts['strand'] = strand

    IntronCounts = IntronCounts[['Chr', 'start', 'end', 'pid', 'gid', 'strand', 'length'] + samples]
    
    IntronCounts.to_csv(output, sep='\t', index=False, header=True)


