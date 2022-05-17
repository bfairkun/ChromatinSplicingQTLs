import numpy as np
import pandas as pd
import os
#from tqdm import tqdm
from sklearn.preprocessing import scale
from scipy.stats import rankdata
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from sklearn.impute import SimpleImputer
import argparse 

def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))


def combine_dataframe(samples, windowStyle, skip_introns, max_missing=0.1, top_introns = 10000, parameter = 'Slope'):
    df = pd.DataFrame()
    df_error = pd.DataFrame()
    for sample in samples:#, leave=True, position=0):
            
        df_sample = pd.read_csv('IntronSlopes/AllSlopes/{sample}.{windowStyle}.tab.gz'.format(
#         df_sample = pd.read_csv('IntronSlopes/slopes/{sample}.{windowStyle}.glm_nb.tab.gz'.format(
            sample=sample, windowStyle=windowStyle), sep='\t', index_col=0).dropna()

#         fdr = multipletests(df_sample['Slope.p.value'], alpha=0.25, method='fdr_bh')
#         df_sample = pd.DataFrame(df_sample.loc[fdr[0], 'Slope'])
        df_sample_error = pd.DataFrame(df_sample['std.error'])
        df_sample_error.columns = [sample]
        df_error = pd.concat([df_error, df_sample_error], axis=1)
        df_sample = pd.DataFrame(df_sample[parameter])
        df_sample.columns = [sample]
        df = pd.concat([df, df_sample], axis=1)
        
    df = df.loc[df.isna().mean(axis=1) <= max_missing]
    df_error = df_error.loc[[x for x in df.index if x not in skip_introns]]
    
    best_samples = df_error.mean(axis=1).sort_values()[:top_introns].index
    
    return df.loc[best_samples]


def qqnorm_slopes(df_slopes):
    
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(df_slopes).T).T
    
    imputed_df.columns = df_slopes.columns 
    imputed_df.index = df_slopes.index
    
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
parser.add_argument('--skip_samples', type=str, required=True)
parser.add_argument('--skip_introns', type=str, required=True)
parser.add_argument('--max_missing', type=str, required=True)
parser.add_argument('--top_introns', type=str, required=True)
# parser.add_argument('--FDR', type=str, required=True)


if __name__ == '__main__':
    
    args = parser.parse_args()
    windowStyle = args.windowStyle
    skip_samples = args.skip_samples.split(':')
    skip_introns = args.skip_introns.split(':')
    max_missing = float(args.max_missing)
    top_introns = int(args.top_introns)
#     FDR = np.float(args.FDR)
    
    samples = [x.split('.')[0] for x in os.listdir("IntronSlopes/AllSlopes/") if "{windowStyle}.tab.gz".format(windowStyle=windowStyle) in x]
    
    samples = [x for x in samples if x not in  skip_samples]

    
    df_slopes = combine_dataframe(samples, windowStyle, skip_introns, max_missing=max_missing, 
                                  top_introns=top_introns, parameter = 'Slope')
    df_qqnorm = qqnorm_slopes(df_slopes)
    bed = intron_list_to_bed(df_qqnorm.index)
    output_bed = pd.concat([bed, df_qqnorm], axis=1).sort_values(['#Chr', 'start'])
    output_bed.to_csv('QTLs/QTLTools/chRNA.Slopes.All/OnlyFirstReps.qqnorm.bed.gz', sep='\t', index=False, header=True)
    
    slope_bed = pd.concat([bed, df_slopes], axis=1).sort_values(['#Chr', 'start'])
    slope_bed.to_csv('QTLs/QTLTools/chRNA.Slopes.All/OnlyFirstReps.Slopes.bed.gz', sep='\t', index=False, header=True)
    
    df_intercept = combine_dataframe(samples, windowStyle, skip_introns, max_missing=max_missing, 
                                  top_introns=top_introns, parameter = '(Intercept)')
    intercept_bed = pd.concat([bed, df_intercept], axis=1).sort_values(['#Chr', 'start'])
    intercept_bed.to_csv('QTLs/QTLTools/chRNA.Slopes.All/OnlyFirstReps.Intercept.bed.gz', sep='\t', index=False, header=True)
    
    
    
        
        
        
        
        
        
