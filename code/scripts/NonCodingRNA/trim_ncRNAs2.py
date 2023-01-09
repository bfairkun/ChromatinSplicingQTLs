import numpy as np
import pandas as pd
import argparse

def process_hmm(chrom, strand):
    

    hmm = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', 
                 names = ['chrom', 'start', 'end', 'names', 'score', 'strand'] + [
                     'chRNA_'+str(i+1) for i in range(86)
                 ]
                     ).sort_values('start')

    hmm.index = hmm.names

    hmm['length'] = hmm.end - hmm.start
    
    return hmm

def get_largest_string_non_zero(x):
    cum_sum = []
    suma= 0
    for i in x:
        if (i > 1.1):
            suma+= 1
        else:
            suma = 0
        cum_sum.append(suma)
    
    max_cum_sum = np.max(cum_sum)
    largest_string_non_zero = max_cum_sum/len(x)
    return largest_string_non_zero

def is_fractured(counts):
    largest_string_non_zero = get_largest_string_non_zero(counts)
    mean_non_zero = np.mean(np.array(counts)>1)
    
    if (mean_non_zero > 0.4) and (largest_string_non_zero < 0.2):
        return True
    else:
        return False
    
def split_and_int(counts):
    x = np.array([int(y) for y in counts.split(',')])
    return x

def is_high_expressed(transcript_counts):
    expression_max = np.max(transcript_counts) >= 4
    expression_median = np.median(transcript_counts) >= 3
    if expression_max:# or expression_median:
        return True
    
def is_very_high(counts):
    expression_median = np.max(counts) >= 4
    if expression_median:
        return True
    
def is_bimodal(counts):
    len_counts = len(counts)
    len_split = int(len_counts/5)
    if (np.median(counts[:len_split]) > 2) and (np.median(counts[-len_split:]) > 2):
        for i in range(len_split, len_counts-len_split):
            if np.median(counts[i:i+len_split]) <= 2:
                return True
    else:
        return False
#     len_split = int(len(counts)/3)
    
#     first = np.median(counts[:len_split]) >= 3
#     mid = np.median(counts[len_split:2*(len_split)]) <= 2
#     last = np.median(counts[2*(len_split):]) >= 3
    
#     if first and last and mid:
#         return True
    
def split_transcript(counts, min_len, zero = 1, grace = 0.1):
    
    grace_len = np.min([int(len(counts)*grace), 50]) 
    print(grace_len)#int(len(counts)*grace)#np.min([int(len(counts)*grace), 20])
    transcript_list = []
    transcript = False
    high_expression = False
    
    grace_counts = 0
    
    consider_grace = False
    
    for i in range(len(counts)):
            
        if counts[i] <= zero:
            if transcript:
                grace_counts += 1
                if (grace_counts > grace_len) or (not consider_grace):
                    end = i #- 1

                    if zero == 1:
                        is_high = is_high_expressed(counts[start:end])
                    else:
                        is_high = is_very_high(counts[start:end])

                    if ((end - start) > min_len) or is_high:
                        transcript_list.append((start, end))
                        high_expression = False
                    transcript = False
                    grace_counts = 0
        else:
            if counts[i] >= 3:
                consider_grace = True
#             grace_counts = 0
            if not transcript:
                start = i
            transcript = True
    if transcript:
        end = len(counts) #- 1
        if zero == 1:
            is_high = is_high_expressed(counts[start:end])
        else:
            is_high = is_very_high(counts[start:end])
        if ((end - start) > min_len) or is_high:
            transcript_list.append((start, end))
            
    if len(transcript_list) == 0:
        transcript_list = [(0, len(counts))]
    return transcript_list
    
    
def easy_trim(counts, min_len_percent = 0.2):
    len_counts = len(counts)
    
    
    if np.max(counts) == 2:
        if is_fractured(counts):
            return [(0, len_counts)]
        else:
            transcript_trimmed = trim_transcript(counts, (0, len_counts))
            return [transcript_trimmed]
    
    min_len = np.min([len_counts*min_len_percent, 500])
    
    transcript_list = split_transcript(counts, min_len)
    
    
    
    transcript_trimmed_list = []
    for transcript in transcript_list:
        transcript_trimmed = trim_transcript(counts, transcript)
        
        high_transcript = is_high_expressed(counts[transcript_trimmed[0]:transcript_trimmed[1]])
        bimodal_transcript = is_bimodal(counts[transcript_trimmed[0]:transcript_trimmed[1]])
        
        if (transcript_trimmed[1] - transcript_trimmed[0] >= min_len) and (not high_transcript): 
            transcript_trimmed_list.append(transcript_trimmed)
        elif high_transcript and (not bimodal_transcript):
            transcript_trimmed_list.append(transcript_trimmed)
        elif high_transcript and bimodal_transcript:
            
            bcounts = counts[transcript_trimmed[0]:transcript_trimmed[1]]
            
            add_start = transcript_trimmed[0]
            
            bimodal_transcript_list = split_transcript(bcounts, min_len, zero = 2)
            
            
            for btranscript in bimodal_transcript_list:
                btranscript = (btranscript[0] + add_start, btranscript[1] + add_start)
                btranscript_trimmed = trim_transcript(counts, btranscript)
                transcript_trimmed_list.append(btranscript_trimmed)
            
            
    if len(transcript_trimmed_list) == 0:
        transcript_trimmed_list = [(0, len_counts)]
    
    return transcript_trimmed_list


def trim_transcript(counts, transcript):
    start, end = transcript
    transcript_counts = counts[start:end]
    if np.median(transcript_counts) == 4:
        zero = 3
    elif is_high_expressed(transcript_counts):
        zero = 2
    else:
        zero = 1
        
    transcript_range = range(start, end)
        
    for i in transcript_range:
        if counts[i] <= zero:
            start += 1
        else:
            break
            
    for i in transcript_range[::-1]:
        if counts[i] <= zero:
            end += -1
        else:
            break
            
    return start, end

def main(chrom, strand):
    hmm = process_hmm(chrom, strand)
    
    samples = ['chRNA_'+str(i+1) for i in range(86)]
    RPKM = ((hmm[samples]/np.array(hmm[samples].sum(axis=0))).T/np.array(hmm.length)).T*1e9
    RPKM.index = hmm.names

    chrom_list = []
    start_list = []
    end_list = []
    names_list = []
    score_list = []
    hmm_strand_list = []

    ncRNAs = pd.Index([x for x in hmm.index if x[:5] == 'ncRNA'])
    ncRNAs = ncRNAs.intersection(RPKM.loc[RPKM.median(axis=1)>0.1].index)

    for ncRNA in ncRNAs:
        start = hmm.loc[ncRNA].start
        end = hmm.loc[ncRNA].end
        hmm_strand = hmm.loc[ncRNA].strand

        counts_string = hmm.loc[ncRNA].score

        print(ncRNA)
        counts = split_and_int(counts_string)

        if hmm_strand == '-':
            counts = counts[::-1]

        trimmed_transcripts = easy_trim(counts)

        print(trimmed_transcripts)

        i = 0

        for transcript in trimmed_transcripts:
            i += 1
            tt_position_start, tt_position_end = transcript



            if hmm_strand == '+':
                start_add = 50*tt_position_start
                end_add = 50*tt_position_end

            else:
                len_counts = len(counts)
                start_add = 50*(len_counts - tt_position_end)
                end_add = 50*(len_counts - tt_position_start)

            tt_start = start + start_add
            tt_end = start + end_add

            chrom_list.append(chrom)
            hmm_strand_list.append(hmm_strand)
            start_list.append(tt_start)
            end_list.append(tt_end)
            names_list.append(ncRNA + '_' + str(i+1))
            score_list.append('.')


    df = pd.DataFrame()
    df['chrom'] = chrom_list
    df['start'] = start_list
    df['end'] = end_list
    df['names'] = names_list
    df['score'] = score_list
    df['strand'] = hmm_strand_list
    
    df.to_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm_trimmed.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', index=False, header=False)



parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--strand', type=str, required=True)

if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    strand = args.strand
    
    main(chrom, strand)