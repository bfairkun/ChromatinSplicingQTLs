import numpy as np
import pandas as pd
import argparse

from sklearn.preprocessing import scale
from scipy.stats import rankdata
from scipy.stats import norm
from sklearn.impute import SimpleImputer

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))

def qqnorm_data(df_data):
    
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(df_data.T)).T
    
    imputed_df.columns = df_data.columns 
    imputed_df.index = df_data.index
    
    df_scale = pd.DataFrame()
    df_qqnorm = pd.DataFrame()
    
    for idx, row in imputed_df.iterrows():#, leave=True, position=0):
        
        scale_row = scale(row)
        df_scale[idx] = scale_row
        
    df_scale.index = imputed_df.columns
    
    df_scale = df_scale.T
    
    for sample in df_scale.columns:
        
        qqnorm_row = qqnorm(df_scale[sample])
        df_qqnorm[sample] = qqnorm_row
        
    df_qqnorm.index = imputed_df.index
    
    return df_qqnorm


parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--min_obs', type=int, required=False, default=10)

if __name__ == '__main__':
    
    args = parser.parse_args()
    input_file = args.input_file
    output = args.output
    min_obs = args.min_obs
    
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    samples = [x for x in df.columns if x != 'NA18855']
    
    df_in_order = pd.DataFrame()
    df_out_of_order = pd.DataFrame()

    for sample in df.columns:
        x = [int(x.split('/')[0]) for x in df[sample]]
        df_in_order[sample] = x

        x = [int(x.split('/')[1]) for x in df[sample]]
        df_out_of_order[sample] = x


    df_in_order.index = df.index
    df_out_of_order.index = df.index
    
    all_reads = df_in_order + df_out_of_order

    
    ratio = df_in_order/all_reads

    idx = ratio.loc[((~ratio[samples].isna()).sum(axis=1) >= min_obs) & (ratio[samples].var(axis=1) > 0), samples].index
    
    X = ratio.loc[idx, samples]
    
    qqnorm = qqnorm_data(X)
    
    
    df_bed = pd.DataFrame(pd.Series(idx).str.split(':', expand=True))
    df_bed.columns = ['#Chr', 'start', 'x', 'y', 'end', 'strand']

    df_bed['pid'] = pd.Series(idx)
    df_bed['gid'] = pd.Series(idx)

    df_bed = df_bed[['#Chr', 'start', 'end', 'pid', 'gid', 'strand']]
    df_bed.index = idx
    
    
    out = pd.concat([df_bed, qqnorm], axis=1)
    
    out.to_csv(output, sep='\t', header=True, index=False)
