import numpy as np
import pandas as pd
import os
import argparse
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import scale
from scipy.stats import rankdata
from scipy.stats import norm

def process_spliceq_IER(spliceq_file):
    
    sample = spliceq_file.split('/')[-1].split('.')[0]
    chIER = pd.read_csv(spliceq_file, sep='\t')
    chIER_ = chIER.drop(columns=['gene_ID', 'transcript_ID', 'intron_ID']).drop_duplicates()
    chIER_ = chIER_.loc[[x not in ['X', 'Y'] for x in chIER_.chr]]

    intronID = []
    for idx, row in chIER_.iterrows():
        intron = 'chr' + str(row.chr) + ':' + str(row.IStart+3) + ':' + str(row.IEnd-3) + ':' + row.strand
        intronID.append(intron)

    chIER_['IntronID'] = intronID
#     chIER_.groupby('IntronID').IER.mean()

    chIER_score = pd.DataFrame(chIER_.groupby('IntronID').IER.mean())
    chIER_score.columns = [sample]
        
    chIER_exon5_cov = pd.DataFrame(chIER_.groupby('IntronID').exon5_cov.max())
    chIER_exon5_cov.columns = [sample]
    
    chIER_sj5_cov_split = pd.DataFrame(chIER_.groupby('IntronID').sj5_cov_split.max())
    chIER_sj5_cov_split.columns = [sample]
    
    chIER_sj5_cov_nonsplit = pd.DataFrame(chIER_.groupby('IntronID').sj5_cov_nonsplit.max())
    chIER_sj5_cov_nonsplit.columns = [sample]
    
    chIER_intron_cov = pd.DataFrame(chIER_.groupby('IntronID').intron_cov.max())
    chIER_intron_cov.columns = [sample]
    
    
    
    chIER_exon3_cov = pd.DataFrame(chIER_.groupby('IntronID').exon3_cov.max())
    chIER_exon3_cov.columns = [sample]
    
    chIER_sj3_cov_split = pd.DataFrame(chIER_.groupby('IntronID').sj3_cov_split.max())
    chIER_sj3_cov_split.columns = [sample]
    
    chIER_sj3_cov_nonsplit = pd.DataFrame(chIER_.groupby('IntronID').sj3_cov_nonsplit.max())
    chIER_sj3_cov_nonsplit.columns = [sample]
    
    return chIER_score, chIER_exon5_cov, chIER_sj5_cov_split, chIER_sj5_cov_nonsplit, chIER_intron_cov, chIER_exon3_cov, chIER_sj3_cov_split, chIER_sj3_cov_nonsplit


def process_spliceq_IRjunctions(spliceq_file):
    
    sample = spliceq_file.split('/')[-1].split('.')[0]
    chIR_f3 = pd.read_csv(spliceq_file, sep='\t')
    chIR_f3 = chIR_f3.drop(columns=['transcript_ID', 'intron_ID']).drop_duplicates()
    
    chIR_f3_pos = chIR_f3.loc[chIR_f3.strand == '+']
    chIR_f3_pos.sj5end += 2
    chIR_f3_pos.sj3start += (-3)
    intronID_pos = []
    for idx, row in chIR_f3_pos.iterrows():
        intid = 'chr' + str(row.chr) + ':' + str(row.sj5end) + ':' + str(row.sj3start) + ':' + row.strand
        intronID_pos.append(intid)

    chIR_f3_pos.index = intronID_pos

    chIR_f3_neg = chIR_f3.loc[chIR_f3.strand == '-']
    chIR_f3_neg.sj3end += 2
    chIR_f3_neg.sj5start += (-3)
    intronID_neg = []
    for idx, row in chIR_f3_neg.iterrows():
        intid = 'chr' + str(row.chr) + ':' + str(row.sj3end) + ':' + str(row.sj5start) + ':' + row.strand
        intronID_neg.append(intid)

    chIR_f3_neg.index = intronID_neg

    chIR_f3_total = pd.concat([chIR_f3_pos, chIR_f3_neg])
    chIR_f3_total['IntronID'] = (intronID_pos + intronID_neg)
    chIR_f3_total = chIR_f3_total.loc[[x not in ['X', 'Y'] for x in chIR_f3_total.chr]]
    
    
    chIR_f3_score = pd.DataFrame(chIR_f3_total.groupby('IntronID').score.mean())
    chIR_f3_score.columns = [sample]
    
    return chIR_f3_score
    
def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))


def qqnormDF(dfIR, max_missing=0.4, skip_samples = []):
    
    samples = [x for x in dfIR.columns if (((x[:2] == 'NA') or (x[:2] == 'HG')) and (x not in skip_samples))]
    
    print(str(dfIR.shape[0]) + " total sites")    
    
    # var min = 0.005; same as in leafcutter
    dfIR = dfIR.loc[dfIR[samples].std(axis=1) >= 0.005, samples]
    print(str(dfIR.shape[0]) + " with + min variance")    
    
    # missing max = 0.4; same as in leafcutter
    dfIR = dfIR.loc[dfIR[samples].isna().mean(axis=1) <= max_missing]
    print(str(dfIR.shape[0]) + " with enough observations") 
    
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(dfIR).T).T
    
    imputed_df.columns = dfIR.columns 
    imputed_df.index = dfIR.index
    
    df_qqnorm = pd.DataFrame()
    
    df_scale = pd.DataFrame(scale(imputed_df[samples], axis=1))
    df_scale.index = imputed_df.index
    df_scale.columns = imputed_df.columns
    
    for sample in df_scale.columns:
        
        qqnorm_col = qqnorm(df_scale[sample])
        df_qqnorm[sample] = qqnorm_col
        
    df_qqnorm.index = imputed_df.index
    
    return df_qqnorm



def process_df_list(df_list):
    IRdf = pd.concat(df_list, axis=1)
    
    samples = list(IRdf.columns)
    
    IRdf['pid'] = list(IRdf.index)
    IRdf['gid'] = list(IRdf.index)
    
    IRdf[['#Chr', 'start', 'end', 'strand']] = IRdf['pid'].str.split(':', expand=True)
        
    column_order = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + samples
        
    IRdf = IRdf[column_order]
    
    return IRdf

parser = argparse.ArgumentParser()
parser.add_argument('--spliceq_dir', type=str, required=True)
parser.add_argument('--spliceq_mode', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--output_qqnorm', type=str, required=True)
parser.add_argument('--skip_samples', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    spliceq_dir = args.spliceq_dir
    spliceq_mode = args.spliceq_mode
    output = args.output
    output_qqnorm = args.output_qqnorm
    skip_samples = args.skip_samples
    
    df_score = []
    df_exon5 = []
    df_sj5 = []
    df_sj5_ncov = []
    df_intron = []
    df_exon3 = []
    df_sj3 = []
    df_sj3_ncov = []
    
    sample_files = os.listdir(spliceq_dir)
    
    for sample in sample_files:
        
        if spliceq_mode == 'IRjunctions':
            df = process_spliceq_IRjunctions(spliceq_dir + '/' + sample)
            df_score.append(df)
        else:
            
            df = process_spliceq_IER(spliceq_dir + '/' + sample)
            df_score.append(df[0])
            df_exon5.append(df[1])
            df_sj5.append(df[2])
            df_sj5_ncov.append(df[3])
            df_intron.append(df[4])
            df_exon3.append(df[5])
            df_sj3.append(df[6])
            df_sj3_ncov.append(df[7])
        
        
#     IRdf = pd.concat(df_list, axis=1)
    
#     samples = list(IRdf.columns)
    
#     IRdf['pid'] = list(IRdf.index)
#     IRdf['gid'] = list(IRdf.index)
    
#     IRdf[['#Chr', 'start', 'end', 'strand']] = IRdf['pid'].str.split(':', expand=True)
        
#     column_order = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + samples
        
#     IRdf = IRdf[column_order]
     
    IRdf = process_df_list(df_score)
    IRdf.to_csv(output, sep='\t', index=False)
    
    IRqqnorm = qqnormDF(IRdf, max_missing=0.1, skip_samples = [skip_samples])
    idx = IRqqnorm.index
    IRqqnorm_out = pd.concat([IRdf.loc[idx, IRdf.columns[:6]], IRqqnorm], axis=1)
    IRqqnorm_out.to_csv(output_qqnorm, sep='\t', index=False)
    
    
    if spliceq_mode == 'IER':
    
        IRdf = process_df_list(df_exon5)
        IRdf.to_csv(output[:-13] + 'exon5.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_sj5)
        IRdf.to_csv(output[:-13] + 'sj5.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_sj5_ncov)
        IRdf.to_csv(output[:-13] + 'sj5_ncov.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_intron)
        IRdf.to_csv(output[:-13] + 'intron.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_exon3)
        IRdf.to_csv(output[:-13] + 'exon3.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_sj3)
        IRdf.to_csv(output[:-13] + 'sj3.bed.gz', sep='\t', index=False)

        IRdf = process_df_list(df_sj3_ncov)
        IRdf.to_csv(output[:-13] + 'sj3_ncov.bed.gz', sep='\t', index=False)

        
    
    
    
    
    
    
    