import numpy as np
import pandas as pd

def read_junc(junc_file, nmd_list):
    junc = pd.read_csv(junc_file, sep='\t', 
                  names = ['chrom', 'start', 'end', 'junction_id', 'read', 'strand', 'annot', 'gene'])
    
    junc['isNMD'] = np.array(junc.annot.isin(nmd_list))
    junc['isProteinCoding'] = np.array(junc.annot == 'protein_coding.gencode')
    junc['count_exons'] = [1]*len(junc)
    
    return junc

def get_counts_and_NMD_reads(junc, juncs_to_remove, max_junc):
    
    remove_junc = np.array(junc.junction_id.isin(juncs_to_remove))
    junc['remove_junc'] = remove_junc
    reads_to_remove = junc.groupby('read').remove_junc.any().loc[
        junc.groupby('read').remove_junc.any()
    ].index
    
    pc_reads = junc.groupby('read').isProteinCoding.any().loc[
        junc.groupby('read').isProteinCoding.any()
    ].index
    
    junctions_per_read = junc.loc[
        junc.read.isin(pc_reads) & (~junc.read.isin(reads_to_remove))
    ].groupby('read').count_exons.sum()
    
    NMD_junctions_in_read = junc.loc[
        junc.read.isin(pc_reads) & (~junc.read.isin(reads_to_remove))
    ].groupby('read').isNMD.any()

    NMD_reads = [NMD_junctions_in_read.loc[
        junctions_per_read.loc[junctions_per_read == i].index
    ].mean() for i in range(1, max_junc)]
    
    
    
    pc_reads_id = junc.groupby('read').isProteinCoding.any().loc[
        junc.groupby('read').isProteinCoding.any()
    ].index


    NMD_reads_id = junc.groupby('read').isNMD.any().loc[
        junc.groupby('read').isNMD.any()
    ].index
    
    
    junc['isProteinCoding_read'] = junc.read.isin(pc_reads_id)
    junc['isNMD_read'] = junc.read.isin(pc_reads_id.intersection(NMD_reads_id))
    
    pc_reads_out = pc_reads.difference(reads_to_remove)
    
    return junc, NMD_reads, junctions_per_read, pc_reads_out

def wrap_junc(junc_file, nmd_list, juncs_to_remove = ['chr17:81510329-81510479:-'], max_junc = 15):
    junc = read_junc(junc_file, nmd_list)
    junc, NMD_reads, x, y = get_counts_and_NMD_reads(junc, juncs_to_remove, max_junc)
    return junc, NMD_reads


def select_random_read(GM, GM_junc_per_read, n):
    GM_df = GM.loc[GM.read == np.random.choice(GM_junc_per_read.loc[GM_junc_per_read == n].index)]
    return GM_df

def subsample_junctions(GM_df, k):
    junc_subset = np.random.choice(np.array(GM_df.annot), k, replace=False)
    return junc_subset

def subsample_is_nmd(GM, GM_junc_per_read, n, k):
    nmd_list = [x for x in GM.annot.unique() if 'nonsense' in x]
    GM_df = select_random_read(GM, GM_junc_per_read, n)
    junc_subset = subsample_junctions(GM_df, k)
    isNMD = np.any([x in nmd_list for x in junc_subset])
    return isNMD

def select_reads_by_number_of_junctions(GM, GM_junctions_per_read, GM_pc_reads, n):
    
    GM_junctions_per_read = GM_junctions_per_read.loc[GM_pc_reads]
    GM_subsample = GM.loc[GM.read.isin(GM_junctions_per_read.loc[
        (GM_junctions_per_read >= n) & (GM_junctions_per_read <= 15)].index)]
    
    is_nmd = []
    
    for x, y in GM_subsample.groupby('read'):
        is_nmd.append(any(np.random.choice(y.isNMD, n, replace=False)))
        
    nmd_mean = np.mean(is_nmd)
    return nmd_mean
        

def get_NMD_avg(junc, juncs_to_remove, max_junc, i):
    
    remove_junc = np.array(junc.junction_id.isin(juncs_to_remove))
    junc['remove_junc'] = remove_junc
    reads_to_remove = junc.groupby('read').remove_junc.any().loc[
        junc.groupby('read').remove_junc.any()
    ].index
    
    pc_reads = junc.groupby('read').isProteinCoding.any().loc[
        junc.groupby('read').isProteinCoding.any()
    ].index
    
    junctions_per_read = junc.loc[
        junc.read.isin(pc_reads) & (~junc.read.isin(reads_to_remove))
    ].groupby('read').count_exons.sum()
    
    NMD_junctions_in_read = junc.loc[
        junc.read.isin(pc_reads) & (~junc.read.isin(reads_to_remove))
    ].groupby('read').isNMD.any()

    NMD_reads = [NMD_junctions_in_read.loc[
        junctions_per_read.loc[junctions_per_read == i].index
    ].mean() for i in range(1, max_junc)]
    
    pc_reads_id = junc.groupby('read').isProteinCoding.any().loc[
        junc.groupby('read').isProteinCoding.any()
    ].index

    NMD_reads_id = junc.groupby('read').isNMD.any().loc[
        junc.groupby('read').isNMD.any()
    ].index
    
    junc['isProteinCoding_read'] = junc.read.isin(pc_reads_id)
    junc['isNMD_read'] = junc.read.isin(pc_reads_id.intersection(NMD_reads_id))
    
    y = select_reads_by_number_of_junctions(junc, junctions_per_read, 
                                    pc_reads.difference(reads_to_remove), i)
    
    return y
    

def subsample_juncs(junc_file, nmd_list, juncs_to_remove):
    
    junc = read_junc(junc_file, nmd_list)
    
    
    
    gm_nmd_avg = []
    for i in range(1, 16):
        y = get_NMD_avg(junc, juncs_to_remove, 16, i)
        gm_nmd_avg.append(y)
    return gm_nmd_avg
    
    
if __name__ == '__main__':
    juncs_to_remove = ['chr17:81510329-81510479:-']
    
    nmd_list = ['nonsense_mediated_decay.pstopcodon',
                'nonsense_mediated_decay.gencode',
                'nonsense_mediated_decay.YN',
                'nonsense_mediated_decay.far5p',
                'nonsense_mediated_decay.far3p',
                'nonsense_mediated_decay.reason2',
                'nonsense_mediated_decay.novel_junctions',
                'nonsense_mediated_decay.reason1',
                'retained_intron.gencode',
                'retained_intron.novel_junctions',
                'processed_transcript.gencode',
                'processed_transcript.novel_junctions']
    
    stable_list = ['stable.NY',
                   'stable.UTR_junction',
                   'stable.YY']

    nmd_avgs = []
    stable_avgs = []

    for i in range(1, 11):
        
        print('Processing G' + str(i))

        x = subsample_juncs('LongReads/Junctions/GM{i}.annotated.junc.gz'.format(i=str(i)), 
                            nmd_list, juncs_to_remove)
        nmd_avgs.append(x)
        
        y = subsample_juncs('LongReads/Junctions/GM{i}.annotated.junc.gz'.format(i=str(i)), 
                        stable_list, juncs_to_remove)
        stable_avgs.append(y)
        
    nmd_df = pd.DataFrame(nmd_avgs)
    nmd_df.columns = ['junc_' + str(i) for i in range(1, 16)]
    nmd_df.index = ['GM' + str(i) for i in range(1, 11)]
    
    stable_df = pd.DataFrame(stable_avgs)
    
    stable_df.columns = ['junc_' + str(i) for i in range(1, 16)]
    stable_df.index = ['GM' + str(i) for i in range(1, 11)]
    
    nmd_df.to_csv('LongReads/Analysis/nmd_reads.tab.gz', sep='\t', index=True, header=True)
    stable_df.to_csv('LongReads/Analysis/stable_reads.tab.gz', sep='\t', index=True, header=True)

    
    
