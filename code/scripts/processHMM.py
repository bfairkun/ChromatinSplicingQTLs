import numpy as np
import pandas as pd

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--nstates', type=int, required=True)


if __name__ == '__main__':

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    nstates = args.nstates

    chr_pred = pd.read_csv(input_file, sep='\t',
                    index_col=0)

    chr_bed =  pd.Series(chr_pred.index).str.split('_', expand=True)
    chr_bed.columns = ['chrom', 'start', 'end', 'strand']
    
    if nstates == 3:
        chr_bed['state'] = list((chr_pred.state*(-1)+3))
    else:
        chr_bed['state'] = list((chr_pred.state*(-1)+2))

    chr_bed.loc[chr_bed.state >= 1, ['chrom', 'start', 'end']].to_csv(output_file, sep='\t', index=False, header=False)
    