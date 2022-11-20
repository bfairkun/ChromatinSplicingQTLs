import numpy as np
import pandas as pd
import os
import argparse

def get_expression_levels(df):
    sorted_states = df.groupby('state').counts.median().sort_values().index
    no_expression = sorted_states[0]
    low_expression = sorted_states[1]
    expression_1 = sorted_states[2]
    expression_2 = sorted_states[3]
    expression_3 = sorted_states[4]
    return no_expression, low_expression, expression_1, expression_2, expression_3

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--strand', type=str, required=True)

if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    strand = args.strand
    
    counts = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz'.format(
        chrom = chrom, strand=strand), index_col=0, sep='\t')
    if strand == 'plus':
        counts_rev = pd.read_csv('NonCodingRNA/tables/{chrom}.minus.counts.tab.gz'.format(
            chrom = chrom), index_col=0, sep='\t')
    else:
        counts_rev = pd.read_csv('NonCodingRNA/tables/{chrom}.plus.counts.tab.gz'.format(
            chrom = chrom), index_col=0, sep='\t')
    
    df = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.predicted_states.tab.gz'.format(
        chrom = chrom, strand=strand), sep='\t', index_col=0)
    
    df['counts'] = np.log1p(counts.median(axis=1))
    
    X = np.array(counts)
    Y = np.array(counts_rev)

    ratio = X/(X+Y+1e-20)
    median_ratio = pd.DataFrame(ratio, index=counts.index, columns=counts.columns).median(axis=1)
    
    no_expression, low_expression, expression_1, expression_2, expression_3 = get_expression_levels(df)
    
    expression = []
    for idx, row in df.iterrows():
        if (row.state == no_expression): #or (median_ratio.loc[idx] <= 0.03):
            expression.append(-1)
        elif (row.state == low_expression):
            expression.append(1)
        elif (row.state == expression_1):
            expression.append(2)
        elif (row.state == expression_2):
            expression.append(3)
        elif (row.state == expression_3):
            expression.append(4)
    
    df['expression'] = expression
    
    df = pd.merge(df, counts, left_index=True, right_index=True)
    samples = list(counts.columns)
    
    X = pd.DataFrame([x.split('_') for x in df.loc[df.expression == 1].index], 
                 columns=['chrom', 'start', 'end', 'strand']
                )

    X.index = df.loc[df.expression == 1].index
    if strand == 'plus':
        X.strand = ['+']*len(X)
    else:
        X.strand = ['-']*len(X)
    X['name'] = df.loc[df.expression == 1].index
    X['score'] = [1] * len(X)

    X = pd.merge(X, df.loc[X.index,samples], 
                 left_index=True, right_index=True
                )
    X[
        ['chrom', 'start', 'end', 'name', 'score', 'strand'] + samples
     ].to_csv('NonCodingRNA/tables/{chrom}.{strand}.low_expression.bed.gz'.format(
        chrom=chrom, strand=strand), 
              sep='\t', header=False, index=False)


    X = pd.DataFrame([x.split('_') for x in df.loc[df.expression >= 2].index], 
                     columns=['chrom', 'start', 'end', 'strand'])
    X.index = df.loc[df.expression >= 2].index
    if strand == 'plus':
        X.strand = ['+']*len(X)
    else:
        X.strand = ['-']*len(X)
    X['name'] = df.loc[df.expression >= 2].index
    X['score'] = df.loc[df.expression >= 2].expression
    X = pd.merge(X, df.loc[X.index,samples], 
                 left_index=True, right_index=True)
    X[
        ['chrom', 'start', 'end', 'name', 'score', 'strand'] + samples
    ].to_csv('NonCodingRNA/tables/{chrom}.{strand}.expression.bed.gz'.format(
        chrom=chrom, strand=strand), sep='\t', header=False, index=False)







    