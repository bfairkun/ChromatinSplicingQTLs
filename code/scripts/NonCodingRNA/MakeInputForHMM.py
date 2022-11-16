################################################################
# Perpare count tables for each chromosome, in bins of 200 bp. #
# 5/31/22                                                      #
#  Carlos F Buen Abad Najar, cnajar@uchicago.edu               #
################################################################

import numpy as np
import pandas as pd
import os
import argparse

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

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
# parser.add_argument('--merge', action='store_true', required=False)


if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    counts_dir = args.counts_dir
    strand = args.strand
    output = args.output
#     merge = args.merge

    assay_list = os.listdir(counts_dir)
    print('Merging data from assays')
    df = pd.DataFrame()
    for assay in assay_list:
        print('Merging ' + assay + ' data')
        df_assay = MakeCountsPerAssay(counts_dir, assay, chrom, strand)
        df = pd.concat([df, df_assay], axis=1)
        
#     if merge:
        
        
        
        # THIS IS FOR 10 GROUPS of 8 and 9
#         df_merged = pd.DataFrame()
#         groups_4 = list(chunks(df.columns[:32], 8))
#         groups_5 = list(chunks(df.columns[32:], 9))
        
#         for i in range(len(groups_4)):
#             col_name = 'chRNA_' + str(i + 1)
#             df_merged[col_name] = df[groups_4[i]].sum(axis=1)
#         for i in range(len(groups_5)):
#             col_name = 'chRNA_' + str(i + 5)
#             df_merged[col_name] = df[groups_5[i]].sum(axis=1)





        # THIS IS FOR 43 Groups of 2
#         df_merged = pd.DataFrame()
#         groups_2 = list(chunks(df.columns, 2))
        
#         for i in range(43):
#             col_name = 'chRNA_' + str(i + 1)
#             df_merged[col_name] = df[groups_2[i]].sum(axis=1)




#         # THIS IS FOR 30 GROUPS of 3 and 2 THIS IS WHAT I WILL MOST LIKELY USE
#         df_merged = pd.DataFrame()
#         groups_4 = list(chunks(df.columns[:85], 3))
#         groups_5 = df.columns[85:]
        
#         for i in range(29):
#             col_name = 'chRNA_' + str(i + 1)
#             df_merged[col_name] = df[groups_4[i]].sum(axis=1)
#         col_name = 'chRNA_30'
#         df_merged[col_name] = df[groups_5].sum(axis=1)

    
        # THIS IS FOR 5 GROUPS of 18 and 17
    df_merged = pd.DataFrame()
    groups_4 = list(chunks(df.columns[:68], 17))
    groups_5 = df.columns[68:]

    for i in range(len(groups_4)):
        col_name = 'chRNA_' + str(i + 1)
        df_merged[col_name] = df[groups_4[i]].sum(axis=1)
    col_name = 'chRNA_5'
    df_merged[col_name] = df[groups_5].sum(axis=1)


    df_merged.to_csv(output + ".counts.tab.gz", sep='\t', index=True, header=True)
        
#     else:
    
    df.to_csv(output + ".counts.all_samples.tab.gz", sep='\t', index=True, header=True)
        
        
        