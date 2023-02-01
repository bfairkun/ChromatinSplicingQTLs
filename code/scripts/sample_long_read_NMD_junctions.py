import numpy as np
import pandas as pd

def read_junc(junc_file, nmd_list):
    junc = pd.read_csv(junc_file, sep='\t', 
                  names = ['chrom', 'start', 'end', 'junction_id', 'read', 'strand', 
                           'annot', 'ensembl', 'gene_id'])
    
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

def wrap_junc(junc_file, nmd_list, juncs_to_remove = ['chr17:81510329-81510479:-'], max_junc = 16):
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

def process_junctions_file(junctions_file, juncs_to_remove):
    
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
    
    ##### No Sample #####
    junc, nmd, = wrap_junc(junctions_file, 
                                    nmd_list, juncs_to_remove = juncs_to_remove)
        
    nmd.append(nmd)

    junc, stable, = wrap_junc(junctions_file, 
                                stable_list, juncs_to_remove = juncs_to_remove)

    stable.append(stable)

    ##### Sample #####

    nmd_avg = subsample_juncs(junctions_file, nmd_list, juncs_to_remove)

    stable_avg = subsample_juncs(junctions_file, stable_list, juncs_to_remove)
    
    return nmd, stable, nmd_avg, stable_avg
    
    
def process_sample_list(sample_list, juncs_to_remove):
    
    nmd_list = []
    stable_list = []
    nmd_avg_list = []
    stable_avg_list = []
    
    for sample in sample_list:
        junc_file = "LongReads/Junctions/{sample}.annotated.junc.gz".format(sample=sample)
        print("Processing " + junc_file)
        nmd, stable, nmd_avg, stable_avg = process_junctions_file(junc_file, juncs_to_remove)
        
        nmd_list.append(nmd)
        stable_list.append(stable)
        nmd_avg_list.append(nmd_avg)
        stable_avg_list.append(stable_avg)
        
    return nmd_list, stable_list, nmd_avg_list, stable_avg_list

def make_df_from_list(percent_list, sample_list):
    df = pd.DataFrame(percent_list)
    df.columns = ['junc_' + str(i+1) for i in range(df.shape[1])]
    df.index = sample_list
    
    return df

def make_dataset_df(sample_list, juncs_to_remove, output_prefix):
    
    nmd_list, stable_list, nmd_avg_list, stable_avg_list = process_sample_list(sample_list, juncs_to_remove)
    
    nmd_df = make_df_from_list(nmd_list, sample_list)
    stable_df = make_df_from_list(stable_list, sample_list)
    nmd_avg_df = make_df_from_list(nmd_avg_list, sample_list)
    stable_avg_df = make_df_from_list(stable_avg_list, sample_list)
    
    nmd_df.to_csv(
        'LongReads/Analysis/{output_prefix}.nmd.tab.gz'.format(output_prefix=output_prefix), 
        sep='\t', index=True, header=True
    )
    
    stable_df.to_csv(
        'LongReads/Analysis/{output_prefix}.stable.tab.gz'.format(output_prefix=output_prefix), 
        sep='\t', index=True, header=True
    )
    
    nmd_avg_df.to_csv(
        'LongReads/Analysis/{output_prefix}.nmd_avg.tab.gz'.format(output_prefix=output_prefix), 
        sep='\t', index=True, header=True
    )
    
    stable_avg_df.to_csv(
        'LongReads/Analysis/{output_prefix}.stable_avg.tab.gz'.format(output_prefix=output_prefix), 
        sep='\t', index=True, header=True
    )
    
    
if __name__ == '__main__':
    
    print('At least this should run')
    
    juncs_to_remove_isoseq = ['chr17:81510329-81510479:-']
    iso_seq_samples = ['GM' + str(i) for i in range(1, 11)]
    print("Processing Iso-seq samples")
    make_dataset_df(iso_seq_samples, juncs_to_remove_isoseq, 'IsoSeq')
    
    juncs_to_remove_NMD_KD = ['chr17:81510329-81510479:-',
                              'chr5:181241639-181242170:+', 
                              'chr5:181241639-181242170:-']
    
    NMD_KD_samples = ['CTRL1_shRNA.SAMEA8691110',
                      'CTRL2_shRNA.SAMEA8691111',
                      'SMG6_shRNA.SAMEA8691112',
                      'SMG6_SMG7_shRNA.SAMEA8691113',
                      'SMG7_shRNA.SAMEA8691114',
                      'UPF1_shRNA.SAMEA8691115']
    
    print("Processing NMD KD samples")
    make_dataset_df(NMD_KD_samples, juncs_to_remove_NMD_KD, 'NMD_KD')
    
    
    
    juncs_to_remove_Churchman = ['chr17:81510329-81510479:-',
                                 'chr5:181241639-181242170:+', 
                                 'chr5:181241639-181242170:-']
    
    Churchman_samples = ['K562_ONT_DMSO_1.SAMN11467430',
                         'K562_ONT_DMSO_2.SAMN11467429',
                         'K562_4sUchr_ONT_1.SAMN10505969',
                         'K562_4sUchr_ONT_2.SAMN10505968',
                         'K562_4sUchr_ONT_3.SAMN10505967',
                         'K562_4sUchr_ONT_4.SAMN12726878',
                         'K562_4sUchr_ONT_5a.SAMN12726877',
                         'K562_4sUchr_ONT_5b.SAMN12726876',
                         'K562_4sUchr_ONT_noA.SAMN10505966']
    
    print("Processing Churchman samples")
    make_dataset_df(Churchman_samples, juncs_to_remove_Churchman, 'Churchman')


    
    
    
    
    
    