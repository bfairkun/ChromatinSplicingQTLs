import numpy as np
import pandas as pd
import os
#from tqdm import tqdm
from scipy.stats import rankdata
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from sklearn.impute import SimpleImputer
import argparse 

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))


def combine_dataframe(samples, windowStyle, max_missing=0.4):
    df = pd.DataFrame()
    for sample in samples:#, leave=True, position=0):
        df_sample = pd.read_csv('IntronSlopes/slopes/{sample}.{windowStyle}.glm_nb.tab.gz'.format(
            sample=sample, windowStyle=windowStyle), sep='\t', index_col=0).dropna()

        fdr = multipletests(df_sample['Slope.p.value'], alpha=0.25, method='fdr_bh')
        df_sample = pd.DataFrame(df_sample.loc[fdr[0], 'Slope'])
        df_sample.columns = [sample]
        df = pd.concat([df, df_sample], axis=1)
        
    df = df.loc[df.isna().mean(axis=1) <= max_missing]
    return df


def qqnorm_slopes(df_slopes):
    
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(df_slopes).T).T
    
    imputed_df.columns = df_slopes.columns 
    imputed_df.index = df_slopes.index
    
    df_qqnorm = pd.DataFrame()
    
    for idx, row in imputed_df.iterrows():#, leave=True, position=0):
        
        qqnorm_row = qqnorm(row)
        
        df_qqnorm[idx] = qqnorm_row
        
    df_qqnorm.index = imputed_df.columns
    
    return df_qqnorm.T


def intron_list_to_bed(intron_list):
    df_list = []
    for intron in intron_list:
        df_list.append(intron.split('_'))
        
    df = pd.DataFrame(df_list)
    df.index = intron_list
    df.columns = ['gid', '#Chr', 'start', 'end', 'strand']
    df['pid'] = intron_list
    df = df[['#Chr', 'start', 'end', 'pid', 'gid', 'strand']]
    return df


parser = argparse.ArgumentParser()
parser.add_argument('--windowStyle', type=str, required=True)
parser.add_argument('--FDR', type=str, required=True)


if __name__ == '__main__':
    
    args = parser.parse_args()
    windowStyle = args.windowStyle
    FDR = np.float(args.FDR)
    
    samples = [x.split('.')[0] for x in os.listdir("IntronSlopes/slopes/") if "{windowStyle}.glm_nb.tab.gz".format(windowStyle=windowStyle) in x]
    
    df_slopes = combine_dataframe(samples, windowStyle)
    df_qqnorm = qqnorm_slopes(df_slopes)

    bed = intron_list_to_bed(df_qqnorm.index)
    
    output_bed = pd.concat([bed, df_qqnorm], axis=1).sort_values(['#Chr', 'start'])
    
    output_bed.to_csv('QTLs/QTLTools/chRNA.Slopes/OnlyFirstReps.qqnorm.bed.gz', sep='\t', index=False, header=True)
        
