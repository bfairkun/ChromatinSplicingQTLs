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

def clean_junc_file(gm, juncs_to_remove):
    PC_list = gm.groupby('read').isProteinCoding.any()
    PC_list = PC_list.loc[PC_list].index
    
    gm['remove'] = gm.junction_id.isin(juncs_to_remove)
    reads_to_remove = gm.groupby('read').remove.any()
    reads_to_remove = reads_to_remove.loc[reads_to_remove].index
    
    reads_to_keep = PC_list.difference(reads_to_remove)
    
    junc_out = gm.loc[gm.read.isin(reads_to_keep)]
    
    return junc_out

def get_expression_quartiles(gm, RPKM):
    
    shared_genes = RPKM.index.intersection(pd.Index(gm.ensembl))
    genes_exp = RPKM[RPKM.columns[5:]].median(axis=1).sort_values()
    Q25 = genes_exp.loc[shared_genes].sort_values().quantile(0.25)
    Q50 = genes_exp.loc[shared_genes].sort_values().quantile(0.5)
    Q75 = genes_exp.loc[shared_genes].sort_values().quantile(0.75)
    
    genes_Q25 = genes_exp.loc[genes_exp <= Q25].index
    genes_Q50 = genes_exp.loc[(genes_exp > Q25) & (genes_exp <= Q50)].index
    genes_Q75 = genes_exp.loc[(genes_exp > Q50) & (genes_exp <= Q75)].index
    genes_Q100 = genes_exp.loc[(genes_exp > Q75)].index
    
    return genes_Q25, genes_Q50, genes_Q75, genes_Q100

def get_quartile_dict(gm, genes_Q, n_samples=20):
    n_list = {}

    for s in range(n_samples):

        n_list_sample = {}

        for read_id, read in gm.loc[gm.ensembl.isin(genes_Q)].groupby('read'):
            juncs_in_read = len(read.count_exons.values)

            read_NMD = read.isNMD.values


            if (juncs_in_read >=2) and (juncs_in_read <= 20):
                for n in range(2, juncs_in_read+1):
        #             for j in range(20):
                    sample_nmd = any(np.random.choice(read_NMD, n, replace=False))
                    if n not in n_list_sample.keys():
                        n_list_sample.update({n:[sample_nmd]})
                    else:
                        n_list_sample[n].append(sample_nmd)
            
        nmd_by_junctions = [np.mean(n_list_sample[k]) for k in sorted(n_list_sample.keys())]
        
        n_list.update({s:nmd_by_junctions})
        
    return n_list

def get_junc_file_subsample(junc_file, nmd_list, RPKM, juncs_to_remove):
    gm = read_junc(junc_file, nmd_list)
    gm = clean_junc_file(gm, juncs_to_remove)
    genes_Q25, genes_Q50, genes_Q75, genes_Q100 = get_expression_quartiles(gm, RPKM)
    
    n_list_Q25 = get_quartile_dict(gm, genes_Q25, n_samples=20)
    n_list_Q50 = get_quartile_dict(gm, genes_Q50, n_samples=20)
    n_list_Q75 = get_quartile_dict(gm, genes_Q75, n_samples=20)
    n_list_Q100 = get_quartile_dict(gm, genes_Q100, n_samples=20)
    
    return n_list_Q25, n_list_Q50, n_list_Q75, n_list_Q100


def make_df_from_samples(Q_samples, name, Q):
    Q_ = pd.DataFrame(Q_samples).T
    junc_list = ['junc_' + str(i+2) for i in range(len(Q_.columns))]
    Q_.columns = junc_list
    Q_['quartile'] = [Q] * 20
    Q_['name'] = [name] * 20
    
    return Q_

def make_df_from_quartile_list(Q_samples_list, name):
    Q_25 = make_df_from_samples(Q_samples_list[0], name, '25')
    Q_50 = make_df_from_samples(Q_samples_list[1], name, '50')
    Q_75 = make_df_from_samples(Q_samples_list[2], name, '75')
    Q_100 = make_df_from_samples(Q_samples_list[3], name, '100')
    
    X = pd.concat([Q_25, Q_50, Q_75, Q_100], axis=0)
    return X

def make_dataset_df(sample_list, nmd_list, juncs_to_remove, RPKM):
    GM_df = pd.DataFrame()
    for GM in sample_list:
        X = get_junc_file_subsample('LongReads/Junctions/{GM}.annotated.junc.gz'.format(GM=GM), 
                                    nmd_list, RPKM, juncs_to_remove)
        df = make_df_from_quartile_list(X, GM)
        GM_df = pd.concat([GM_df, df], axis=0)
        
    return GM_df
        


if __name__ == '__main__':
    
    print('Reading RPKM table')
    
    RPKM = pd.read_csv(
        'QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstRepsUnstandardized.qqnorm.bed.gz', 
        sep='\t', index_col=3
    )
    
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
    
    print('Processing Iso-seq data...')
    juncs_to_remove = ['chr17:81510329-81510479:-']
    GM_samples = ['GM' + str(i) for i in range(1, 11)]
    GM_df = make_dataset_df(GM_samples, nmd_list, juncs_to_remove, RPKM)
    print('Done. Saving.')
    GM_df.to_csv('LongReads/Analysis/IsoSeq.ByQuartile.tab.gz', sep='\t', index=False, header=True)
    
    juncs_to_remove = ['chr17:81510329-81510479:-',
                       'chr5:181241639-181242170:+', 
                       'chr5:181241639-181242170:-']
    
    GM_samples = ['CTRL1_shRNA.SAMEA8691110',
                  'CTRL2_shRNA.SAMEA8691111',
                  'SMG6_shRNA.SAMEA8691112',
                  'SMG6_SMG7_shRNA.SAMEA8691113',
                  'SMG7_shRNA.SAMEA8691114',
                  'UPF1_shRNA.SAMEA8691115']
    
    print('Processing NMD KD data...')
    
    GM_df = make_dataset_df(GM_samples, nmd_list, juncs_to_remove, RPKM)
    GM_df.to_csv('LongReads/Analysis/NMD_KD.ByQuartile.tab.gz', sep='\t', index=False, header=True)
    
    GM_samples = ['K562_ONT_DMSO_1.SAMN11467430',
                  'K562_ONT_DMSO_2.SAMN11467429',
                  'K562_4sUchr_ONT_1.SAMN10505969',
                  'K562_4sUchr_ONT_2.SAMN10505968',
                  'K562_4sUchr_ONT_3.SAMN10505967',
                  'K562_4sUchr_ONT_4.SAMN12726878',
                  'K562_4sUchr_ONT_5a.SAMN12726877',
                  'K562_4sUchr_ONT_5b.SAMN12726876',
                  'K562_4sUchr_ONT_noA.SAMN10505966']
    
    print('Processing Churchman data...')
    
    GM_df = make_dataset_df(GM_samples, nmd_list, juncs_to_remove, RPKM)
    GM_df.to_csv('LongReads/Analysis/Churchman.ByQuartile.tab.gz', sep='\t', index=False, header=True)
    
    