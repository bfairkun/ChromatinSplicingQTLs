import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    input_file = args.input
    output = args.output

    df = pd.read_csv(input_file, sep=' ')
    samples = list(df.columns[1:])
    df.columns = ['cluster_id'] + samples

    df[['chrom', 'start', 'end', 'cluster']] = df['cluster_id'].str.split(':', 3, expand=True)

    df['strand'] = [x.split('_')[-1] for x in df.cluster]
    df['junction_id'] = df.chrom + ':' + df.start + ':' + df.end + ':' + df.strand

    for sample in samples:
        df[[sample, sample+'_total']] = df[sample].str.split('/', 1, expand=True)

    out = df[['chrom', 'start', 'end', 'junction_id', 'cluster_id', 'strand'] + samples]
    out.columns = ['#Chr', 'start', 'end', 'junction_id', 'cluster_id', 'strand'] + samples
    out.to_csv(output, sep='\t', index=False, header=True)