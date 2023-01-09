import numpy as np
import pandas as pd
import argparse

def get_tail(hmm, ncRNA, strict=False):
    expression = [int(x) for x in hmm.loc[ncRNA].score.split(',')]
    if (np.mean(expression) > 1.5) and (len(expression) > 1000):
        trim_length = trim_expression(expression, strict)
        return trim_length
    else:
        return False
    
def trim_expression(expression_list, strict):
    trim_length = 0
    while True:
        if strict:
            if np.mean(np.array(expression_list[-100:])) < 1.1:
                expression_list = expression_list[:-100]
                trim_length += 5000
            else:
                break
        else:
            
            if all(np.array(expression_list[-100:]) == 1):
                expression_list = expression_list[:-100]
                trim_length += 5000
            else:
                break
    return trim_length

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

def is_fractured(transcripto):
    x = [int(x) for x in transcripto.split(',')]
    if np.max(x) >= 3.9:
        return False
    largest_string_non_zero = get_largest_string_non_zero(x)
    mean_non_zero = np.mean(np.array(x)>1)
    
    if (mean_non_zero > 0.4) and (largest_string_non_zero < 0.2):
        return True
    else:
        return False

def trim_and_split(hmm, ncRNA, split_len = 50, min_len=200, min_portion = 0.2, strict=True,
                  strand = 'plus'):
    
    current_pos = 0
    start_list = []
    end_list = []
    
    transcript_length = hmm.loc[ncRNA].end - hmm.loc[ncRNA].start
    
    transcripto = hmm.loc[ncRNA].score
    
#     if strand == 'minus':
#         transcripto = transcripto[::-1]
    
    if strict:
        transcripto_split = transcripto.split(',')
        counts = [int(x) for x in transcripto_split]
        
        
        largest_non_zero = get_largest_string_non_zero(counts)
        percent_zero = np.mean(np.array(counts) == 1)
        
        
#         if is_fractured(counts):#np.mean(counts) > 1.5:
#             split_len = 1e100
        
        #if (np.max(counts)==4):
        if (np.max(counts)==4) or (np.quantile(counts, 0.5)>=2.9):
            if np.quantile(counts, 0.75) <= 3.1:
                transcripto = transcripto.replace('2,'*100, '3,'*100)
                transcripto_split = transcripto.split(',')
            transcripto_split = ['1' if x == '2' else x for x in transcripto_split]
            if np.quantile(counts, 0.75) >= 3.9:
                print('enormous!')
                transcripto_split = ['1' if x == '3' else x for x in transcripto_split]
        elif (largest_non_zero > 0.25) and (percent_zero > 0.5):
            split_len = 5
            
        transcripto = ','.join(transcripto_split)
        
    if is_fractured(transcripto):#np.mean(counts) > 1.5:
        print('fractured!')
        end = len(transcripto.split(','))*50
        return [0], [end]
    split_transcripts = transcripto.split(','.join(['1']*split_len))
    
    longest = max([len(x) for x in split_transcripts])
    
    high_scores = []
    for transcript in split_transcripts:
        try:
            sc = np.mean([int(x) for x in transcript.strip(',').split(',')])
            high_scores.append(sc)
        except:
            continue
    highest = np.max(high_scores)
    inherited = False
    print(counts)
    print(split_transcripts)
    len_to_add = 0
    for transcript in split_transcripts:
        current_longest = ((len(transcript) == longest) and (len(transcript) > 3))
        
        if transcript == ',':
            current_pos += (50*split_len)
        elif transcript == '':
            current_pos += (50*split_len)
        else:
            
            transcript_counts = transcript.strip(',').split(',')
            if (transcript[0] == ',') and (strand == 'plus'):
                current_pos += 50
#             if (strand == 'minus') and inherited:
#                 current_pos += (50*split_len)

                
            transcript_counts = [int(x) for x in transcript_counts]
        
            print(transcript_counts)
            
            current_highest = np.abs(np.mean(transcript_counts) - highest) <= 1e-100
            
            is_max_peak = transcript_counts.count(4) >= 5
        
            transcript_portion = (len(transcript_counts)*50)/transcript_length
            if ((len(transcript_counts) > min_len) and (transcript_portion >= min_portion)) or is_max_peak or current_highest or current_longest:
                
                start_pos = current_pos

                for i in range(len(transcript_counts)):
                    
                    if transcript_counts[i] == 1:
                        print('trimmin')
                        start_pos += 50
                        
                    else:
                        break

                if np.mean(transcript_counts) > 1.1:
                    start_list.append(start_pos)
                    
                trim = trim_expression(transcript_counts, True)
                
                if trim:
                    if np.mean(transcript_counts) > 1.1:
                        end_list.append(current_pos + ((len(transcript_counts)*50)-trim))
                        
                else:
                    if np.mean(transcript_counts) > 1.1:
                        trimmer = 50
#                         if strand == 'plus':
#                             trimmer = 50
#                         elif strand == 'minus':
#                             trimmer = 0
#                         else:
#                             raise Exception('strand error')
                        end_list.append(current_pos + (len(transcript_counts)*50) - trimmer)

            print('adding')
            if transcript[-1] == ',':
                current_pos += (50*split_len)
            current_pos += (len(transcript_counts)*50)
            
            if transcript[-1] == ',':
                inherited = True
            
#     if strand == 'minus':
#         return end_list, start_list
        
    return start_list, end_list

# def trim_and_split(hmm, ncRNA, split_len = 50, min_len=200, min_portion = 0.2, strict=True,
#                   strand = 'plus'):
    
#     current_pos = 0
#     start_list = []
#     end_list = []
    
#     transcript_length = hmm.loc[ncRNA].end - hmm.loc[ncRNA].start
    
#     transcripto = hmm.loc[ncRNA].score
    
#     if strand == 'minus':
#         transcripto = transcripto[::-1]
    
#     if strict:
#         transcripto_split = transcripto.split(',')
#         counts = [int(x) for x in transcripto_split]
#         #if (np.max(counts)==4):
#         if (np.max(counts)==4) or (np.quantile(counts, 0.75)>=3):
#             transcripto_split = ['1' if x == '2' else x for x in transcripto_split]
#         transcripto = ','.join(transcripto_split)
        
#     split_transcripts = transcripto.split(','.join(['1']*split_len))
    
#     longest = max([len(x) for x in split_transcripts])
    
#     high_scores = []
#     for transcript in split_transcripts:
#         try:
#             sc = np.mean([int(x) for x in transcript.strip(',').split(',')])
#             high_scores.append(sc)
#         except:
#             continue
#     highest = np.max(high_scores)
#     inherited = False
#     print(counts)
#     print(split_transcripts)
#     len_to_add = 0
#     for transcript in split_transcripts:
#         current_longest = ((len(transcript) == longest) and (len(transcript) > 3))
        
#         if transcript == ',':
#             current_pos += (50*split_len)
#         elif transcript == '':
#             current_pos += (50*split_len)
#         else:
            
#             transcript_counts = transcript.strip(',').split(',')
#             if (transcript[0] == ',') and (strand == 'plus'):
#                 current_pos += 50
#             if (strand == 'minus') and inherited:
#                 current_pos += (50*split_len)
# #             if (transcript[-1] == ',') and (strand == '-'):
# #                 current_pos += (50*split_len)
                
#             transcript_counts = [int(x) for x in transcript_counts]
        
#             print(transcript_counts)
            
#             current_highest = np.abs(np.mean(transcript_counts) - highest) <= 1e-100
        
#             transcript_portion = (len(transcript_counts)*50)/transcript_length
#             if ((len(transcript_counts) > min_len) and (transcript_portion >= min_portion)) or current_highest or current_longest:
                
#                 start_pos = current_pos

#                 for i in range(len(transcript_counts)):
                    
#                     if transcript_counts[i] == 1:
#                         print('trimmin')
#                         start_pos += 50
                        
#                     else:
#                         break

#                 if np.mean(transcript_counts) > 1.1:
#                     start_list.append(start_pos)
                    
#                 trim = trim_expression(transcript_counts, True)
                
#                 if trim:
#                     if np.mean(transcript_counts) > 1.1:
#                         end_list.append(current_pos + ((len(transcript_counts)*50)-trim))
                        
#                 else:
#                     if np.mean(transcript_counts) > 1.1:
#                         if strand == 'plus':
#                             trimmer = 50
#                         elif strand == 'minus':
#                             trimmer = 0
#                         else:
#                             raise Exception('strand error')
#                         end_list.append(current_pos + (len(transcript_counts)*50) - trimmer)

#             print('adding')
#             if transcript[-1] == ',':
#                 current_pos += (50*split_len)
#             current_pos += (len(transcript_counts)*50)
            
#             if transcript[-1] == ',':
#                 inherited = True
            
    
        
#     return start_list, end_list

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True)
parser.add_argument('--strand', type=str, required=True)


if __name__ == '__main__':

    args = parser.parse_args()
    chrom = args.chrom
    strand = args.strand
    
    hmm = pd.read_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', 
                 names = ['chrom', 'start', 'end', 'names', 'score', 'strand'] + [
                     'chRNA_'+str(i+1) for i in range(86)
                 ]
                     ).sort_values('start')

    hmm.index = hmm.names

    hmm['length'] = hmm.end - hmm.start
    samples = ['chRNA_'+str(i+1) for i in range(86)]
    RPKM = ((hmm[samples]/np.array(hmm[samples].sum(axis=0))).T/np.array(hmm.length)).T*1e9
    RPKM.index = hmm.names
    
    chrom_list = []
    start_list = []
    end_list = []
    names_list = []
    score_list = []
    strand_list = []

    ncRNAs = pd.Index([x for x in hmm.index if x[:5] == 'ncRNA'])
    ncRNAs = ncRNAs.intersection(RPKM.loc[RPKM.median(axis=1)>0.1].index)
    
    
    
    for n in ncRNAs:
        start = hmm.loc[n].start
        transcript_end = hmm.loc[n].end
        scores = [int(x) for x in hmm.loc[n].score.split(',')]
        transcript_ln = transcript_end - start
        ln = len(scores)
        if (RPKM.loc[n].median() >= 1) or (ln > 100):
#             if (ln < 750):
            if (ln < 650):
                srt, end = trim_and_split(hmm, n, min_len=50, strict=True, split_len=10, strand=strand)
#             elif (ln >= 750):
            else:
                srt, end = trim_and_split(hmm, n, min_len=50, strict=True, split_len=100, strand=strand)
        else:
            srt, end = [], []
        
        print('trimming: ' + n)
        print(srt, end)
        print('')
        if len(srt) > 0:
            
            for i in range(len(srt)):
            
#                 if strand == 'plus':

                start_ = start + srt[i]
                end_ = start + end[i]
                if strand == 'plus':
                    strand_list.append('+')
                else:
#                     start_ = start + (transcript_ln - end[i])
#                     print(end[i])
#                     end_ = start + (transcript_ln - srt[i])
#                     print(srt[i])
                    strand_list.append('-')
                    
                chrom_list.append(chrom)
                start_list.append(start_)
                end_list.append(end_)
                names_list.append(n + '_' + str(i+1))
                score_list.append('.')
        else:
            chrom_list.append(chrom)
            start_list.append(start)
            end_list.append(transcript_end)
            names_list.append(n)
            score_list.append('.')
            if strand == 'plus':
                strand_list.append('+')
            else:
                strand_list.append('-')
            
#     pc_names = pd.Index([x for x in hmm.index if x[:5] != 'ncRNA'])
#     chrom_list += list(hmm.loc[pc_names].chrom)
#     start_list += list(hmm.loc[pc_names].start)
#     end_list += list(hmm.loc[pc_names].end)
#     names_list += list(hmm.loc[pc_names].names)
#     score_list += list(hmm.loc[pc_names].score)
#     strand_list += list(hmm.loc[pc_names].strand)
    
    df = pd.DataFrame()
    df['chrom'] = chrom_list
    df['start'] = start_list
    df['end'] = end_list
    df['names'] = names_list
    df['score'] = score_list
    df['strand'] = strand_list
    
    df.to_csv('NonCodingRNA/tables/{chrom}.{strand}.hmm_trimmed.bed.gz'.format(
        chrom = chrom, strand = strand), sep='\t', index=False, header=False)



