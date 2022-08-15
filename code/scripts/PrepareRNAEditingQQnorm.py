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
    
    df_scale = pd.DataFrame(scale(imputed_df, axis=1))
    df_scale.columns = imputed_df.columns
    df_scale.index = imputed_df.index
    
    df_qqnorm = pd.DataFrame()
    
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
    samples = [x for x in df.columns if ((x[:2] == 'NA') and (x != 'NA18855'))]
    
    idx = df.loc[((~df[samples].isna()).sum(axis=1) >= min_obs) & (df[samples].var(axis=1) > 0), samples].index
    
    X = df.loc[idx, samples]
    
    qqnorm = qqnorm_data(X)
    
    
    df_bed = df.loc[idx, ['chr', 'BEDstart' 'BEDend' 'name' 'score', 'strand']]

    df_bed['pid'] = df_bed.name
    df_bed['gid'] = df_bed.name

    df_bed = df_bed[['chr', 'BEDstart' 'BEDend' 'pid' 'gid', 'strand']]
    
    
    df_bed.columns = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand']
    
    out = pd.concat([df_bed, qqnorm], axis=1)
    
    out.to_csv(output, sep='\t', header=True, index=False)
