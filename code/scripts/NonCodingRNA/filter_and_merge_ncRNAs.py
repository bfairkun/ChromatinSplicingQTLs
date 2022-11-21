import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import argparse

import warnings
warnings.filterwarnings("error")

def get_distance(ncRNA_strand, ncRNA1, ncRNA2, strand):
    
    if strand == 'plus':
        FivePrimeEnd = int(ncRNA_strand.loc[ncRNA1].end)
        ThreePrimeStart = int(ncRNA_strand.loc[ncRNA2].start)
        distance = ThreePrimeStart - FivePrimeEnd
    elif strand == 'minus':
        FivePrimeEnd = int(ncRNA_strand.loc[ncRNA1].start)
        ThreePrimeStart = int(ncRNA_strand.loc[ncRNA2].end)
        distance = FivePrimeEnd - ThreePrimeStart
    else:
        raise Exception('invalid strand')
    return distance

def get_RPKM_ratio(RPKM, ncRNA1, ncRNA2):
    rpkm1 = np.array(RPKM.loc[ncRNA1])
    rpkm2 = np.array(RPKM.loc[ncRNA2])
    
    ratio = np.median(rpkm1/(rpkm2+1e-20))
    
    return ratio

def get_correlation(RPKM, ncRNA1, ncRNA2):
    try:
        r = pearsonr(RPKM.loc[ncRNA1], RPKM.loc[ncRNA2])[0]
    except:
        
        if RPKM.loc[ncRNA2].max() == 0:
            raise Exception (ncRNA2 + ' should have been removed')
        else:
            print('RPKM for ' + ncRNA1)
            print(RPKM.loc[ncRNA1].max())
            r = -1
    return r

def get_hmm_rev_ratio(hmm, hmm_rev, annot, samples):
    if annot[:3] == 'pc_':
        x = 1
    else:
        x = (hmm.loc[annot, samples]/hmm_rev.loc[annot, samples]).median()
    return x


def filter_trailing_ncRNA(hmm, hmm_rev, ncRNA_RPKM, distance_thre, cor_thre, RPKM_ratio_thre, strand, verbose=True):
    
    
    
    if strand == 'minus':
        hmm = hmm.sort_values('end', ascending=False)
    else:
        hmm = hmm.sort_values('start', ascending=True)
        
    hmm_list = hmm.index
    
    trail_list = []
    
    anchor = hmm_list[0]
    current_trail = hmm_list[0]
    current_tail = [hmm_list[0]]
    merge_list = []
    current_merge = []
    
    for ncRNA in hmm_list[1:]:
        
        if verbose:
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('##############################')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('filtering: ' + ncRNA)
            print('anchor: ' + anchor)
            print('current trail: ' + current_trail)
        
        anchor_chrom = hmm.loc[anchor, 'chrom']
        ncRNA_chrom = hmm.loc[ncRNA, 'chrom']
        
        if (anchor_chrom != ncRNA_chrom):# or (ncRNA[:5] != 'ncRNA'):
            anchor = ncRNA
            current_trail = ncRNA
            continue
            
        if RPKM.loc[ncRNA].max() == 0:
            if verbose:
                print('removing ' + ncRNA + ' for 0 expression')
            continue
            
        if strand == 'plus':
            
            if int(hmm.loc[anchor, 'end']) > int(hmm.loc[ncRNA, 'end']):
                if verbose:
                    print('removing ' + ncRNA + ' for overlap')
                trail_list.append(ncRNA)
                continue

            if int(hmm.loc[current_trail, 'end']) > int(hmm.loc[ncRNA, 'end']):
                if verbose:
                    print('removing ' + ncRNA + ' for overlap')
                trail_list.append(ncRNA)
                continue
                
            
        elif strand == 'minus':
            if int(hmm.loc[anchor, 'start']) < int(hmm.loc[ncRNA, 'start']):
                if verbose:
                    print('removing ' + ncRNA + ' for overlap')
                trail_list.append(ncRNA)
                continue

            if int(hmm.loc[current_trail, 'start']) < int(hmm.loc[ncRNA, 'start']):
                if verbose:
                    print('removing ' + ncRNA + ' for overlap')
                trail_list.append(ncRNA)
                continue
        else:
            raise Exception('invalid strand')
            
        distance =  get_distance(hmm, current_trail, ncRNA, strand)
        distance_from_anchor =  get_distance(hmm, anchor, ncRNA, strand)
        corr = get_correlation(ncRNA_RPKM, anchor, ncRNA)
        corr_trail = get_correlation(ncRNA_RPKM, current_trail, ncRNA)
        ratio = get_RPKM_ratio(ncRNA_RPKM, anchor, ncRNA)
        ratio_trail = get_RPKM_ratio(ncRNA_RPKM, current_trail, ncRNA)
        RPKM_current_ncRNA = ncRNA_RPKM.loc[ncRNA].median()
        RPKM_anchor = ncRNA_RPKM.loc[anchor].median()
        RPKM_current_trail = ncRNA_RPKM.loc[current_trail].median()
        
        corr_tail = [get_correlation(ncRNA_RPKM, x, ncRNA) for x in current_tail]
        
        samples = ['chRNA_'+str(i+1) for i in range(86)]
        
        current_trail_rev_ratio = get_hmm_rev_ratio(hmm, hmm_rev, current_trail, samples)
        ncRNA_rev_ratio = get_hmm_rev_ratio(hmm, hmm_rev, ncRNA, samples)
        
        token = True
        is_trail = False
        
        #remove_distance = distance_thre
            
#         max_merge_ratio = 10
#         min_merge_ratio = 0.25
        
        
        if len(current_merge) >= 2:
            current_merge_anchor = current_merge[0]
        else:
            current_merge_anchor = current_trail
        if strand == 'plus':
            len_trail = hmm.loc[current_trail].end - hmm.loc[current_merge_anchor].start
        elif strand == 'minus':
            len_trail = hmm.loc[current_merge_anchor].end - hmm.loc[current_trail].start
            
        
        if ratio > 10:
            remove_distance = 10000
        else:
            remove_distance = 3000
            
        remove_distance = np.max([remove_distance, len_trail])
        if corr >= 0.75:
            remove_distance = np.min([remove_distance, 250000])
        else:
            remove_distance = np.min([remove_distance, 50000])
          
        low_reverse =  (ncRNA_rev_ratio < 0.1) and (current_trail_rev_ratio < 0.1)
        low_expression = (RPKM_current_ncRNA <= 0.5) and (RPKM_current_trail <= 0.5)
        long_trail = (len_trail > 100000)
        
        high_corr = ((corr >= 0.5) and (corr_trail >= 0.5)) or (corr >= 0.75)
        high_ratio = (ratio >= 50)
        
        
        
        high_corr_and_ratio = high_corr and high_ratio
        
        anchor_is_pc = anchor[:3]=='pc_'
        
        corr_min = 0.2
        
#         if anchor_is_pc or (RPKM_current_ncRNA < 5):
#             RPKM_ratio_thre = 1
#         else:
        RPKM_ratio_thre = 2
    
        high_corr_long_trail = high_corr and (len_trail > 10000)
            
             
        if high_ratio or long_trail:
            merge_distance = 20000
        elif low_reverse or low_expression or high_corr_long_trail:
            merge_distance = 5000
        else:
            merge_distance = 1000
#             remove_distance = 100000
#             corr_min = 0
        
        if low_reverse or low_expression:

            RPKM_ratio_thre = 1
            max_merge_ratio = 1e100
            min_merge_ratio = -1
            merge_corr_min = 0.1
        
        else:
            merge_corr_min = 0.1
            max_merge_ratio = 1e2
            min_merge_ratio = 0.1
            
            
        
        
        if verbose:

            print('reverse ratio: ' + str(ncRNA_rev_ratio))
            print('distance: ' + str(distance))
            print('distance from anchor: ' + str(distance_from_anchor))
            print('correlation: ' + str(corr))
            print('correlation tail: ' + str(np.mean(corr_tail)))
            print('ratio anchor: ' + str(ratio))
            print('ratio trail: ' + str(ratio_trail))
            print('correlation trail: ' + str(corr_tail))
            print('length trail: ' + str(len_trail))
            print('remove distance: ' + str(remove_distance))
            print('RPKM: ' + str(RPKM_current_ncRNA))
            print('RPKM anchor: ' + str(RPKM_anchor))
            print('merge distance: ' + str(merge_distance))
            print('min_merge_ratio: ' + str(min_merge_ratio))
            
            print('max_merge_ratio: ' + str(max_merge_ratio))
            print('correlation min: ' + str(merge_corr_min))
            print('')
            
        rev_ratio_too_low = ncRNA_rev_ratio < 0.03
        print("criteria 0: reverse ratio is too low " + str(rev_ratio_too_low))
                
        too_close_to_pc = ( (distance < 1) and ((anchor[:5] != 'ncRNA') or (current_trail[:5] != 'ncRNA')))
        print("criteria 1: It's extension of pc gene " + str(too_close_to_pc))
        
        standard_criteria =  (distance < remove_distance and ( (corr >= corr_min) or (np.mean(corr_tail) >= corr_min)) and (ratio >= RPKM_ratio_thre) and (ratio_trail >= 0.1)) #any([x >= corr_min for x in corr_tail])       ) 
        
        print('distance:')
        print(distance)
        print(remove_distance)
        print('correlation:')
        print(corr)
        print(corr_tail)
        print(corr_min)
        print('ratio:')
        print(ratio)
        print(ratio_trail)
        print(RPKM_ratio_thre)
        print((ratio >= RPKM_ratio_thre))
        print("criteria 2: thresholds " + str(standard_criteria))
                              
        
        
#         is_trail = is_trail or ((distance < distance_thre) and any([x >= cor_thre for x in corr_tail]) and (ratio >= RPKM_ratio_thre) and (ratio_trail >= 0.5))
#         print('criteria 2: ' + str(is_trail))
#         is_trail = is_trail or ((distance < distance_thre) and any([x >= cor_thre for x in corr_tail]) and (ratio >= 1) and (ratio_trail >= 0.1) and (anchor[:3]=='pc_'))
#         print('criteria 3: ' + str(is_trail))
#         is_trail = is_trail or ((distance < remove_distance) and ((corr > 0) or (corr_trail > 0)) and (((ratio >= 50) and (ratio_trail >= 0.1)) or ((len_trail > 100000) and (ratio > 1))))
#         print('criteria 4: ' + str(is_trail))
#         is_trail = is_trail or ((distance < remove_distance) and ((corr >= 0.5) or (corr_trail >= 0.5)) and (ratio >= 2) and (ratio_trail >= 0.1))
#         print('criteria 5: ' + str(is_trail))
#         is_trail = is_trail or ((distance < remove_distance) and ((corr >= 0.5) or (corr_trail >= 0.5)) and ((ratio >= 1) and (anchor[:3]=='pc_')) and (ratio_trail >= 0.1))
#         print('criteria 6: ' + str(is_trail))
       
        too_close_high_ratio = ((distance < 1000) and high_ratio)
        print('criteria 3: too close and ratio is too high: ' + str(too_close_high_ratio))
        
        is_trail = rev_ratio_too_low or too_close_to_pc or standard_criteria or too_close_high_ratio
        
        print('is trail: ' + str(is_trail))
        print('')
        
        m = merge_pair(hmm, current_merge, current_trail, ncRNA, ncRNA_RPKM, strand = strand, 
                   distance_max=merge_distance, corr_min = merge_corr_min, RPKM_max = 1e100,
                   max_ratio = max_merge_ratio, min_ratio = min_merge_ratio)
                              
        
        if m and (current_trail not in trail_list) and token:
            print('Merge: ' + str(m))    
            if current_trail not in current_merge:
                current_merge.append(current_trail)
                
            current_merge.append(ncRNA)
            if verbose:
                print('Merging pair')
                print(current_trail, ncRNA)
                print('')
            
            current_trail = ncRNA
            current_tail.append(ncRNA)
            #anchor = ncRNA
            #current_tail[ncRNA]
            
        elif is_trail:
            print("Removing because it's trail")
                              
            trail_list.append(ncRNA)
            current_trail = ncRNA
            current_tail.append(ncRNA)
            merge_list.append(current_merge)
            if len(current_merge)> 0:
                print('Adding current merge:')
                print(current_merge)
            current_merge = []
        else:
            anchor = ncRNA
            current_trail = ncRNA
            current_tail = [ncRNA]
            merge_list.append(current_merge)
            current_merge = []
            
    merge_list.append(current_merge)
            
    merge_list = [x for x in merge_list if len(x) > 1]
    
    return trail_list, merge_list


def merge_pair(hmm, current_merge, ncRNA_1, ncRNA_2, RPKM, strand, 
               distance_max=10000, corr_min = 0.2, RPKM_max = 1e100, verbose=True,
               max_ratio = 10, min_ratio = 0.5):
                              
    
    if (ncRNA_1[:5] != 'ncRNA') or (ncRNA_2[:5] != 'ncRNA'):
        if verbose:
            print('One of the genes is not an ncRNA.')
                              
            
        return False
    
    length_1 = hmm.loc[ncRNA_1].length
    length_2 = hmm.loc[ncRNA_2].length
    
    RPKM_med1 = RPKM.loc[ncRNA_1].median()
    RPKM_med2 = RPKM.loc[ncRNA_2].median()
    
    RPKM_ratio = get_RPKM_ratio(RPKM, ncRNA_1, ncRNA_2)
    
    distance = get_distance(hmm, ncRNA_1, ncRNA_2, strand)
    
    corr_list = [get_correlation(RPKM, x, ncRNA_2) for x in current_merge]
    corr_list.append(get_correlation(RPKM, ncRNA_1, ncRNA_2))
    corr = np.mean(corr_list)
#     corr = get_correlation(RPKM, ncRNA_1, ncRNA_2)
    
    merge = (distance <= distance_max)
    merge = merge and (RPKM_med1 < RPKM_max) and (RPKM_med2 < RPKM_max)
#     merge = merge and (RPKM_med1 > 0) and (RPKM_med2 > 0)
    merge = merge and (corr >= corr_min)
    
    merge = merge and (RPKM_ratio < max_ratio) and (RPKM_ratio > min_ratio)
                              
    
    if verbose:
        print('Check if merging {nc1} and {nc2}:'.format(nc1=ncRNA_1, nc2=ncRNA_2))
        print('distance = ' + str(distance) + ', max = ' + str(distance_max))
        print(distance <= distance_max)
        print('correlation = ' + str(corr) + ', min = ' + str(corr_min))
        print(corr >= corr_min)
        print('RPKM ' + str(ncRNA_1) + ' = ' +str(RPKM_med1) + ', max = ' + str(RPKM_max))
        print(RPKM_med1 < RPKM_max)
        print(RPKM_med1 > 0)
        print('RPKM ' + str(ncRNA_2) + ' = ' +str(RPKM_med2) + ', max = ' + str(RPKM_max))
        print(RPKM_med2 < RPKM_max)
        print(RPKM_med2 > 0)
                              
        print('RPKM ratio = ' +str(RPKM_ratio) + ', max = {m1}, min = {m2}'.format(m1=max_ratio,
                                                                                   m2=min_ratio))
                              
        print(RPKM_ratio < max_ratio)
        print(RPKM_ratio > min_ratio)
        print('merge: ' + str(merge))
                              
              
                              
    
    return merge
                              
def merge_ncrnas(hmm, ncrnas, RPKM, strand, verbose=True):
    
    names = []
    chrom = []
    start = []
    end = []
    score = []
    strand_ = []
    
    list_of_merges = [item for sublist in ncrnas for item in sublist]
    
    no_merge = [x for x in hmm.index if x not in list_of_merges]
    
    for nc in ncrnas:
        print('')
        print('merging list')
        print(nc[0])
        print(nc)
        if strand == 'minus':
            start.append(hmm.loc[nc[-1]].start)
            end.append(hmm.loc[nc[0]].end)
            
        elif strand == 'plus':
            start.append(hmm.loc[nc[0]].start)
            end.append(hmm.loc[nc[-1]].end)
            
            
            
        names.append(nc[0])
        chrom.append(hmm.loc[nc[0]].chrom)
        strand_.append(hmm.loc[nc[0]].strand)
        RPKM_med = RPKM.loc[nc].median(axis=1).mean()
        score.append(RPKM_med)
        
        
    for nc in no_merge:
        names.append(nc)
        chrom.append(hmm.loc[nc].chrom)
        strand_.append(hmm.loc[nc].strand)
        score.append(RPKM.loc[nc].median())
        start.append(hmm.loc[nc].start)
        end.append(hmm.loc[nc].end)
    
    out_hmm = pd.DataFrame()
    out_hmm['chrom'] = chrom
    out_hmm['start'] = start
    out_hmm['end'] = end
    out_hmm['name'] = names
    out_hmm['score'] = score
    out_hmm['strand'] = strand_
    
    return out_hmm

"""
def merge_ncrnas(hmm, ncrnas, RPKM, strand, verbose=True):
    
    names = []
    chrom = []
    start = []
    end = []
    score = []
    strand_ = []
    
    current_merge_ncrna = ncrnas[0][0]
    current_merge_list = ncrnas[0]
    
    i = 0
    token = False
    for idx, row in hmm.iterrows():
        if idx == current_merge_ncrna:
            if verbose:
                print('')
                print('merging list')
                print(idx)
                print(current_merge_list)
            names.append(idx)
            chrom.append(row.chrom)
            if strand == 'minus':
                start.append(hmm.loc[current_merge_list[-1]].start)
                end.append(row.end)
            elif strand == 'plus':
                start.append(row.start)
                end.append(hmm.loc[current_merge_list[-1]].end)
            else:
                raise Exception('strand bug')
            RPKM_med = RPKM.loc[current_merge_list].median(axis=1).mean()
            score.append(RPKM_med)
            strand_.append(row.strand)
                
        elif idx in current_merge_list:
            token = True
            continue
        else:
            if token and (len(ncrnas) >= (i+2)):
                i += 1
                current_merge_ncrna = ncrnas[i][0]
                current_merge_list = ncrnas[i]
                token = False

            names.append(idx)
            chrom.append(row.chrom)
            start.append(row.start)
            end.append(row.end)
            score.append(RPKM.loc[idx].median())
            strand_.append(row.strand)
    
    out_hmm = pd.DataFrame()
    out_hmm['chrom'] = chrom
    out_hmm['start'] = start
    out_hmm['end'] = end
    out_hmm['name'] = names
    out_hmm['score'] = score
    out_hmm['strand'] = strand_
    return out_hmm
"""
    
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--strand', type=str, required=True)

if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    strand = args.strand
    
    samples = ['chRNA_'+str(i+1) for i in range(86)]
    
    hmm_nc = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm_w_counts.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', 
                     names = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + samples
                     )
    
    hmm_rev = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.reverse_expression.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', 
                     names = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + samples
                     )
    
    hmm_rev.index = hmm_rev.name
    
#     ncRNAs = pd.Index([x for x in hmm.index if x[:5]=='ncRNA'])
    
    hmm_pc = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', 
                     names = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + samples
                     )
    
    hmm_pc = hmm_pc.loc[[x[:3] == 'pc_' for x in hmm_pc.name]]
    
    hmm = pd.concat([hmm_nc, hmm_pc], axis=0, ignore_index=True)
    hmm.index = hmm.name
    
    if strand == 'plus':
        hmm = hmm.sort_values('start')
    elif strand == 'minus':
        hmm = hmm.sort_values('end', ascending=False)
    else:
        raise Exception('strand bug')


    hmm['length'] = hmm.end - hmm.start
    
    RPKM = ((hmm[samples]/np.array(hmm[samples].sum(axis=0))).T/np.array(hmm.length)).T*1e9
    RPKM.index = hmm.name
    
    print(RPKM.shape)
    
    samples = ['chRNA_'+str(i+1) for i in range(86)]
    
    RPKM = ((hmm[samples]/np.array(hmm[samples].sum(axis=0))).T/np.array(hmm.length)).T*1e9
    RPKM.index = hmm.name
    
    trail_list, merge_list = filter_trailing_ncRNA(hmm, hmm_rev, RPKM, 10000, 0.3, 2, strand=strand, verbose = True)
    
    print('')
    print('Merge list:')
    print(merge_list)
    
    print('')
    print('Trail list:')
    print(trail_list)
    
    ncRNAs = pd.Index([x for x in hmm.index if x[:5]=='ncRNA'])
    
    if strand == 'plus':
        kept = hmm.loc[
            hmm.index.difference(pd.Index(trail_list)).intersection(ncRNAs)
        ].sort_values('start').index
    else:
        kept = hmm.loc[
            hmm.index.difference(pd.Index(trail_list)).intersection(ncRNAs)
        ].sort_values('end', ascending=False).index
    
    hmm_out = merge_ncrnas(hmm.loc[kept], merge_list, RPKM, strand=strand)
    hmm_out = hmm_out.sort_values('start')
    hmm_out = hmm_out.loc[hmm_out.score >= 0.1]
    hmm_out = hmm_out.loc[(hmm_out.end-hmm_out.start)>=500]
    hmm_out = hmm_out.loc[((hmm_out.end - hmm_out.start) >= 1000) | (hmm_out.score >= 1)]
#     (((hmm_.end - hmm_.start) < 1000) & (hmm_.score < 1)) | (((hmm_.end-hmm_.start)<500)
    
    
    hmm_out.to_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm.filtered.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t',
                                        index = False, header=False)
    
    

    

    