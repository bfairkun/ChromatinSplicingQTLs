import pandas as pd
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    
    qqnorm_df = pd.read_csv(input_file, sep='\t')
    samples = pd.read_csv('../data/igsr_samples.tsv.gz', sep='\t', index_col=0)
    YRI_samples = qqnorm_df.columns.intersection(samples.loc[samples['Superpopulation code'] == 'EUR'].index)
    YRI_columns = list(qqnorm_df.columns[:6]) + list(YRI_samples)
    qqnorm_df[YRI_columns].to_csv(output_file, sep='\t', index=False, header=True)

