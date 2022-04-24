import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from scipy.stats import rankdata
from scipy.stats import norm
from sklearn.impute import SimpleImputer
import argparse 


def SummarizeSS(df, sample):
    df_split = df[sample].str.split('/', expand=True).astype(int)
    df_split.columns = ['num', 'den']
    df_split['pid'] = df.pid
    
    num = df_split[['num', 'pid']].groupby('pid', sort=False).sum('num').num
    den = df_split[['den', 'pid']].groupby('pid', sort=False).max('den').den
    
    SS_usage = num/den
    
    num = list(num)
    den = list(den)
    
    # adding pseudocount of 0.5
    SS_usage_for_qqnorm = [(num[i]+0.5)/(den[i]+0.5) if den[i] > 0 else np.nan for i in range(len(num))]
    
    SS_usage = pd.DataFrame(SS_usage)
    SS_usage.columns = [sample]
    
    SS_usage_for_qqnorm = pd.DataFrame(SS_usage_for_qqnorm)
    SS_usage_for_qqnorm.columns = [sample]
    
#     print(SS_usage.shape)
#     print(SS_usage_for_qqnorm.shape)
    
    return SS_usage, SS_usage_for_qqnorm


def process_introns(counts, samples):
    
    PSI_bed, out_bed = SummarizeSS(counts, samples[0])
    
#     PSI_bed = pd.DataFrame()
#     out_bed = pd.DataFrame()
    
    for sample in samples[1:]:
        SS_usage, SS_usage_for_qqnorm = SummarizeSS(counts, sample)
        PSI_bed = pd.concat((PSI_bed, SS_usage), axis=1)
        out_bed = pd.concat((out_bed, SS_usage_for_qqnorm), axis=1)
#         PSI_bed[sample] = SS_usage
#         out_bed[sample] = SS_usage_for_qqnorm
        
    out_bed.index = PSI_bed.index

        
    pid = list(PSI_bed.index)
    
    SS_loc = list(out_bed.index.str.split(':'))
    chrom = [x[0] for x in SS_loc]
    gid = [x[0] + ':' + x[2] for x in SS_loc]
    
    bp_loc = [int(x[1]) for x in SS_loc]
    
    strand = [x[2].split('_')[-1] for x in SS_loc]
    
    chrom = pd.DataFrame(chrom, index = out_bed.index, columns = ['#Chr'])
    bp_loc_start = pd.DataFrame(bp_loc, index = out_bed.index, columns = ['start'])
    bp_loc_end = pd.DataFrame(bp_loc, index = out_bed.index, columns = ['end'])
    pid = pd.DataFrame(pid, index = out_bed.index, columns = ['pid'])
    gid = pd.DataFrame(gid, index = out_bed.index, columns = ['gid'])
    strand = pd.DataFrame(strand, index = out_bed.index, columns = ['strand'])
    
    out_bed = pd.concat((out_bed, chrom), axis=1)
    out_bed = pd.concat((out_bed, bp_loc_start), axis=1)
    out_bed = pd.concat((out_bed, bp_loc_end), axis=1)
    out_bed = pd.concat((out_bed, pid), axis=1)
    out_bed = pd.concat((out_bed, gid), axis=1)
    out_bed = pd.concat((out_bed, strand), axis=1)
#     out_bed['#Chr'] = chrom
#     out_bed['start'] = bp_loc
#     out_bed['end'] = bp_loc
#     out_bed['pid'] = pid
#     out_bed['gid'] = gid
    
#     out_bed['strand'] = strand
    
    out_bed = out_bed[['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + samples]
    
    return out_bed, PSI_bed
    


def ProcessLeafcutterBed(counts):
    
    samples = [x for x in counts.columns if ((x[:2] == 'NA') or (x[:2] == 'HG'))]
    
    bed = counts.chrom.str.split(':', expand=True)
    bed.columns = ['#Chr', 'start', 'end', 'cluster']
    bed['strand'] = list(bed.cluster.str.split('_', expand=True)[2])
    
    pid_start = [row['#Chr'] + ':' + row.start + ':' + row.cluster for idx, row in bed.iterrows()]
    pid_end = [row['#Chr'] + ':' + row.end + ':' + row.cluster for idx, row in bed.iterrows()]
    
    gid = [row['#Chr'] + ':' + row.cluster for idx, row in bed.iterrows()]
    
    bed['pid_start'] = pid_start
    bed['pid_end'] = pid_end
    bed['gid'] = gid
    
    bed = bed[['#Chr', 'start', 'end', 'pid_start', 'pid_end', 'gid', 'strand']]
    
    AltStart = bed.groupby('gid').pid_start.nunique()
    AltEnd = bed.groupby('gid').pid_end.nunique()
    
    StartClusters = AltStart.loc[(AltStart > 1)].index
    EndClusters = AltEnd.loc[(AltEnd > 1)].index
    
    StartIdx = bed.loc[[x in StartClusters for x in bed.gid]].index
    EndIdx = bed.loc[[x in EndClusters for x in bed.gid]].index
    
    StartCounts = counts.loc[StartIdx, samples]
    EndCounts = counts.loc[EndIdx, samples]
    
    StartCounts['pid'] = bed.pid_start
    EndCounts['pid'] = bed.pid_end
    
    print("Processing 5' positions")
    
    StartBed, StartPSI = process_introns(StartCounts, samples)
    
    print("Processing 3' positions")
    EndBed, EndPSI = process_introns(EndCounts, samples)
    
    StartBed.start = StartBed.start.astype(int)
    StartBed.end = StartBed.end.astype(int)

    EndBed.start = EndBed.start.astype(int)
    EndBed.end = EndBed.end.astype(int)
    
    StartBed.start += -3
    StartBed.end += 6

    EndBed.start += -7
    EndBed.end += 2
    
    FivePrimeSS = pd.concat([StartBed.loc[StartBed.strand == '+'], EndBed.loc[EndBed.strand == '-']], axis=0)
    ThreePrimeSS = pd.concat([StartBed.loc[StartBed.strand == '-'], EndBed.loc[EndBed.strand == '+']], axis=0)
    
    FivePrimePSI = pd.concat([StartPSI.loc[StartBed.strand == '+'], EndPSI.loc[EndBed.strand == '-']], axis=0)
    ThreePrimePSI = pd.concat([StartPSI.loc[StartBed.strand == '-'], EndPSI.loc[EndBed.strand == '+']], axis=0)
    
    return FivePrimeSS, ThreePrimeSS, FivePrimePSI, ThreePrimePSI


def qqnorm(x):
    n=len(x)
    a=3.0/8.0 if n<=10 else 0.5
    return(norm.ppf( (rankdata(x)-a)/(n+1.0-2.0*a) ))


def qqnormSpliceSite(dfSS, skip_samples = []):
    
    samples = [x for x in dfSS.columns if (((x[:2] == 'NA') or (x[:2] == 'HG')) and (x not in skip_samples))]
    
    print(str(dfSS.shape[0]) + " total sites")    
    
    # var min = 0.005; same as in leafcutter
    dfSS = dfSS.loc[dfSS[samples].std(axis=1) >= 0.005, samples]
    print(str(dfSS.shape[0]) + " with + min variance")    
    
    # missing max = 0.4; same as in leafcutter
    dfSS = dfSS.loc[dfSS[samples].isna().mean(axis=1) <= 0.4]
    print(str(dfSS.shape[0]) + " with < 40% missing values") 
    
    imputer = SimpleImputer()
    imputed_df = pd.DataFrame(imputer.fit_transform(dfSS).T).T
    
    imputed_df.columns = dfSS.columns 
    imputed_df.index = dfSS.index
    
    df_qqnorm = pd.DataFrame()
    
    df_scale = pd.DataFrame(scale(imputed_df[samples], axis=1))
    df_scale.index = imputed_df.index
    df_scale.columns = imputed_df.columns
    
    for sample in df_scale.columns:
        
        qqnorm_col = qqnorm(df_scale[sample])
        df_qqnorm[sample] = qqnorm_col
        
    df_qqnorm.index = imputed_df.index
    
    return df_qqnorm


parser = argparse.ArgumentParser()
parser.add_argument('--counts', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--skip_samples', type=str, default = '', required=False)

if __name__ == '__main__':
    
    args = parser.parse_args()
    counts = args.counts
    output = args.output
    skip_samples = args.skip_samples.split(':')
    
    print('Loading counts from ' + counts)
    
    counts = pd.read_csv(counts, sep=' ')
    
    
    normPsi5SS, normPsi3SS, Psi5SS, Psi3SS = ProcessLeafcutterBed(counts)
    
    print("Filtering 5' splice sites")
    
    qqnorm5SS = qqnormSpliceSite(normPsi5SS, skip_samples = skip_samples)
    
    print("Filtering 3' splice sites")
    
    qqnorm3SS = qqnormSpliceSite(normPsi3SS, skip_samples = skip_samples)
    
    bed5SS = normPsi5SS.loc[qqnorm5SS.index, normPsi5SS.columns[:6]]
    bed3SS = normPsi3SS.loc[qqnorm3SS.index, normPsi3SS.columns[:6]]
    
    qqnorm5SS = pd.concat([bed5SS, qqnorm5SS], axis=1)
    qqnorm3SS = pd.concat([bed3SS, qqnorm3SS], axis=1)
    
    output_5SS = output + '.5PrimeSS/'
    output_3SS = output + '.3PrimeSS/'
    
    qqnorm5SS.to_csv(output_5SS + 'OnlyFirstReps.qqnorm.bed.gz', sep = '\t', index = False, header = True)
    qqnorm3SS.to_csv(output_3SS + 'OnlyFirstReps.qqnorm.bed.gz', sep = '\t', index = False, header = True)
    
    Psi5SS.to_csv(output_5SS + 'OnlyFirstReps.PSI.bed.gz', sep = '\t', index = False, header = True)
    Psi3SS.to_csv(output_3SS + 'OnlyFirstReps.PSI.bed.gz', sep = '\t', index = False, header = True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    