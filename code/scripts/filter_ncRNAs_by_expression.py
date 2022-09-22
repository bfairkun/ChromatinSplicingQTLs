import numpy as np
import pandas as pd

if __name__ == '__main__':
    
    dir_name = "NonCodingRNA_merged"
    ncRNA_bed = pd.read_csv(dir_name + '/annotation/ncRNA.bed.gz', 
                             sep='\t', names = ['chrom', 'start', 'end', 'name', 'id', 'strand'])
    ncRNA_bed.index = ncRNA_bed.name

    RPKM = pd.read_csv(dir_name + "/Expression_HMM/OnlyFirstReps.RPKM.bed.gz", sep='\t')
    RPKM.index = RPKM.pid

    samples = [x for x in RPKM.columns if ((x[:2]=='NA') and (x!='NA18855'))]


    featureCounts = pd.read_csv('NonCodingRNA_merged/Expression_ncRNA_and_reverse/chRNA.Expression_featureCounts/Counts.txt',
                        sep='\t', comment='#', index_col=0)

    bed = featureCounts[featureCounts.columns[:5]]
    counts = featureCounts[[x for x in featureCounts.columns[5:] if x.split('/')[4] == '1']]
    counts.columns = [x.split('/')[3] for x in counts.columns]

    counts_reverse = counts.loc[[x for x in counts.index if x.split('_')[-1] == 'reverse']]
    counts = counts.loc[counts.index[:len(counts_reverse.index)]]
    counts_reverse.index = counts.index

    annotation = pd.read_csv('NonCodingRNA_merged/annotation/ncRNA.histone.tab.gz', sep='\t', index_col=0)

    incRNA = annotation.loc[['incRNA' in x for x in annotation.rna_type]].index
    uaRNA = annotation.loc[['uaRNA' in x for x in annotation.rna_type]].index
    srtRNA = annotation.loc[['srtRNA' in x for x in annotation.rna_type]].index
    rtRNA = annotation.loc[['rtRNA' in x for x in annotation.rna_type]].index.difference(srtRNA)
    coRNA = annotation.loc[['coRNA' in x for x in annotation.rna_type]].index.difference(uaRNA).difference(incRNA).difference(rtRNA).difference(srtRNA)


    other = annotation.index.difference(incRNA).difference(uaRNA).difference(srtRNA).difference(rtRNA).difference(coRNA)

    mean_ratio = (counts/(counts+counts_reverse)).mean(axis=1)
    median_ratio = (counts/(counts+counts_reverse)).median(axis=1)

    length = ncRNA_bed.loc[annotation.index].end - ncRNA_bed.loc[annotation.index].start

    idx_all = list(incRNA) + list(uaRNA) + list(srtRNA) + list(rtRNA) + list(coRNA) + list(other)
    idx_type = (['incRNA']*len(incRNA)) + (['uaRNA']*len(uaRNA)) + (['srtRNA']*len(srtRNA)) + (['rtRNA']*len(rtRNA))
    idx_type += (['coRNA']*len(coRNA)) + (['annotated']*len(other))
    RPKM_median = np.log10(RPKM[samples].loc[idx_all]+1e-3).median(axis=1)
    RPKM_q = np.log10(RPKM[samples].loc[idx_all]+1e-3).quantile(0.88, axis=1) # 10 samples
    ratio = median_ratio.loc[idx_all]

    df = pd.DataFrame()
    df['rna_type'] = idx_type
    df['RPKM'] = list(RPKM_median)
    df['RPKM_q'] = list(RPKM_q)
    df['strand/reverse ratio'] = list(ratio)
    df.index = idx_all
    df['pass_filter'] = (df.RPKM > -2) & (df['strand/reverse ratio'] >= 0.1)
    df['length'] = np.log10(np.array(list(bed.loc[idx_all].End - bed.loc[idx_all].Start)))
    df['pass_filter_strict'] = ((df.RPKM_q > np.log10(0.05)) & (df['strand/reverse ratio'] >= 0.25))
    # For srtRNAs we require a higher ratio and a RPKM > 0.05 in at least 10 samples
    df['selected'] = ((df.pass_filter) & (df.pass_filter_strict) & (df.rna_type=='srtRNA')) | ((df.pass_filter) & (df.rna_type!='srtRNA'))


    
    ncRNA_bed.loc[df.selected].to_csv(dir_name + '/annotation/NonCodingRNA.bed.gz', sep='\t', index=False, header=False)
    annotation.loc[df.selected].to_csv(dir_name + '/annotation/NonCodingRNA.annotation.tab.gz', sep='\t', index=True, header=True)
