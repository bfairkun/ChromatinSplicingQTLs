import numpy as np
import pandas as pd
import os

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--min_RPKM',  type=float, required=True)
parser.add_argument('--distance',  type=int, required=True)
parser.add_argument('--min_correlation',  type=float, required=True)
parser.add_argument('--RPKM_ratio',  type=float, required=True)
parser.add_argument('--merge', action='store_true', required=False)


from scipy.stats import pearsonr, spearmanr 
def get_ncRNA_distance(ncRNA_strand, ncRNA1, ncRNA2, strand):
    
    if strand == '+':
        FivePrimeEnd = int(ncRNA_strand.loc[ncRNA1].end)
        ThreePrimeStart = int(ncRNA_strand.loc[ncRNA2].start)
        distance = ThreePrimeStart - FivePrimeEnd
    else:
        FivePrimeEnd = int(ncRNA_strand.loc[ncRNA1].start)
        ThreePrimeStart = int(ncRNA_strand.loc[ncRNA2].end)
        distance = FivePrimeEnd - ThreePrimeStart
    return distance

def get_ncRNA_RPKM_ratio(ncRNA_RPKM_median, ncRNA1, ncRNA2):
    rpkm1 = np.array(np.exp(ncRNA_RPKM_median.loc[ncRNA1]))
    rpkm2 = np.array(np.exp(ncRNA_RPKM_median.loc[ncRNA2]))
    
    ratio = np.median(rpkm1/rpkm2)
    
    return ratio

def get_ncRNA_RPKM_correlation(qqnorm, ncRNA1, ncRNA2):
    r = pearsonr(qqnorm.loc[ncRNA1], qqnorm.loc[ncRNA2])[0]
    return r


def filter_trailing_ncRNA(ncRNA_bed, ncRNA_RPKM, ncRNA_qqnorm, distance_thre, cor_thre, RPKM_ratio_thre):
    ncRNA_plus = ncRNA_bed.loc[(ncRNA_bed.strand == '+') & (ncRNA_bed.chrom != 'chrY')]
    ncRNA_minus = ncRNA_bed.loc[(ncRNA_bed.strand == '-') & (ncRNA_bed.chrom != 'chrY')]
#     ncRNA_corr = ncRNA_qqnorm.T.corr()
    ncRNA_RPKM_median = ncRNA_RPKM.median(axis=1)
    
    ncRNA_plus_list = ncRNA_plus.index
    ncRNA_minus_list = ncRNA_minus.index[::-1]
    
    plus_trail = []
    
    anchor = ncRNA_plus_list[0]
    current_trail = ncRNA_plus_list[0]
    
    for ncRNA in ncRNA_plus_list[1:]:
        
        if ncRNA[:3] == 'MIR':
            continue
            
        
        
        anchor_chrom = ncRNA_plus.loc[anchor, 'chrom']
        ncRNA_chrom = ncRNA_plus.loc[ncRNA, 'chrom']
        
        
        #if (anchor_chrom != ncRNA_chrom) or (ncRNA[:5] != 'ncRNA'):
        #    anchor = ncRNA
        #    current_trail = ncRNA
        #    continue
            
        ##################################################################
            
        if (anchor_chrom != ncRNA_chrom):# or (ncRNA[:5] != 'ncRNA'):
            anchor = ncRNA
            current_trail = ncRNA
            continue
            
        if int(ncRNA_plus.loc[anchor, 'end']) > int(ncRNA_plus.loc[ncRNA, 'end']):
            continue
            
        if ncRNA[:5] != 'ncRNA':
            anchor = ncRNA
            current_trail = ncRNA
            continue
        
        ###############################################################3
            
        distance =  get_ncRNA_distance(ncRNA_plus, current_trail, ncRNA, '+')
        corr = get_ncRNA_RPKM_correlation(ncRNA_qqnorm, anchor, ncRNA)
        corr_trail = get_ncRNA_RPKM_correlation(ncRNA_qqnorm, current_trail, ncRNA)
        ratio = get_ncRNA_RPKM_ratio(ncRNA_RPKM, anchor, ncRNA)
        
        
        is_trail = (distance <= distance_thre) and ((corr >= cor_thre) or (corr_trail >= cor_thre)) and (ratio >= RPKM_ratio_thre)
        is_trail = is_trail or ((distance < 1000) and ((corr >= 0.2) or (corr_trail >= 0.2)) and (ratio >= 10))
        
        if is_trail:
            plus_trail.append(ncRNA)
            current_trail = ncRNA
        else:
            anchor = ncRNA
            current_trail = ncRNA
            
            
    #######
    
    minus_trail = []
    
    anchor = ncRNA_minus_list[0]
    current_trail = ncRNA_minus_list[0]
    
    for ncRNA in ncRNA_minus_list[1:]:
        
        if ncRNA[:3] == 'MIR':
            continue
        
        anchor_chrom = ncRNA_minus.loc[anchor, 'chrom']
        ncRNA_chrom = ncRNA_minus.loc[ncRNA, 'chrom']
        
#         if len(ncRNA_chrom) > 1:
#             continue
        
        #if (anchor_chrom != ncRNA_chrom) or (ncRNA[:5] != 'ncRNA'):
        #    anchor = ncRNA
        #    current_trail = ncRNA
        #    continue
        
        #if (anchor_chrom != ncRNA_chrom) or (ncRNA[:5] != 'ncRNA'):
        #    anchor = ncRNA
        #    current_trail = ncRNA
        #    continue
            
        ##################################################################
            
        if (anchor_chrom != ncRNA_chrom):# or (ncRNA[:5] != 'ncRNA'):
            anchor = ncRNA
            current_trail = ncRNA
            continue
            
        if int(ncRNA_minus.loc[anchor, 'start']) < int(ncRNA_minus.loc[ncRNA, 'start']):
            continue
            
        if ncRNA[:5] != 'ncRNA':
            anchor = ncRNA
            current_trail = ncRNA
            continue
        
        ###############################################################3
            
        distance =  get_ncRNA_distance(ncRNA_minus, current_trail, ncRNA, '-')
        corr = get_ncRNA_RPKM_correlation(ncRNA_qqnorm, anchor, ncRNA)
        corr_trail = get_ncRNA_RPKM_correlation(ncRNA_qqnorm, current_trail, ncRNA)
        ratio = get_ncRNA_RPKM_ratio(ncRNA_RPKM, anchor, ncRNA)
        
        
        is_trail = (distance <= distance_thre) and ((corr >= cor_thre) or (corr_trail >= cor_thre)) and (ratio >= RPKM_ratio_thre)
        is_trail = is_trail or ((distance < 1000) and ((corr >= 0.2) or (corr_trail >= 0.2)) and (ratio >= 10))
        
        if is_trail:
            minus_trail.append(ncRNA)
            current_trail = ncRNA
        else:
            anchor = ncRNA
            current_trail = ncRNA
            
    ncRNA_trail = plus_trail + minus_trail
            
    return ncRNA_trail #plus_trail, minus_trail
        

if __name__ == '__main__':
    args = parser.parse_args()
    min_RPKM = args.min_RPKM
    distance = args.distance
    min_correlation = args.min_correlation
    RPKM_ratio = args.RPKM_ratio
    merge = args.merge
    
    if merge:
        dir_name = "NonCodingRNA_merged"
    else:
        dir_name = "NonCodingRNA_annotation"
    
    ncRNA_bed = pd.read_csv(dir_name + '/annotation/allHMM.merged.bed.gz', 
                         sep='\t', names = ['chrom', 'start', 'end', 'name', 'id', 'strand'])
    ncRNA_bed.index = ncRNA_bed.name

    RPKM = pd.read_csv(dir_name + "/Expression_HMM/OnlyFirstReps.RPKM.bed.gz", sep='\t')
    RPKM.index = RPKM.pid

    samples = [x for x in RPKM.columns if ((x[:2]=='NA') and (x!='NA18855'))]
    
    
    
    qqnorm = pd.read_csv(dir_name + "/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz", sep='\t')
    qqnorm.index = qqnorm.pid
    
    
    ncRNA_trail = filter_trailing_ncRNA(ncRNA_bed, RPKM[samples], qqnorm[samples], distance_thre = distance, 
                                   cor_thre = min_correlation, RPKM_ratio_thre = RPKM_ratio)
    
    print(ncRNA_trail)
    filter1 = pd.Index([x for x in RPKM.index if x not in ncRNA_trail])
    filter2 = filter1[np.exp(RPKM.loc[filter1, samples]).quantile(0.9, axis=1) >= min_RPKM]
    
    ncRNAs = [x for x in filter2 if x[:5] == 'ncRNA']
    
    ncRNA_filtered_bed = ncRNA_bed.loc[ncRNAs]
    ncRNA_filtered_RPKM  = RPKM.loc[filter2]
    
    ncRNA_filtered_bed.to_csv(dir_name + "/annotation/ncRNA.bed.gz", sep='\t', index=False, header=False)
    #ncRNA_filtered_RPKM.to_csv("QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.RPKM.filtered.bed.gz", 
    #                           sep='\t', index=False, header=True)
    
    
    
    
