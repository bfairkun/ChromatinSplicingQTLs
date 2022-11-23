import numpy as np
import pandas as pd
import os

def read_bed(bed_file, flag='-wa'):
    if flag == '-wa':
        bed = pd.read_csv(bed_file, sep='\t',
                    names = ['chrom', 'start', 'end', 'gene_name', 'id', 'strand'])
    elif flag == '-wo':
        bed = pd.read_csv(bed_file, sep='\t', 
                         names = ['chrom', 'start', 'end', 'gene_name', 'id', 'strand',
                                  'chrom_', 'start_', 'end_', 'gene_name_', 'id_', 'strand_','overlap'])
        
    return bed


def check_true_uaRNA(bed, tss_bed, uaRNA, gene, verbose=False):
    
    strand = bed.loc[uaRNA].strand
    
    if gene[:5] != 'ncRNA':
    
    
        if strand == '+':
            # 3' end of uaRNA on + strand
            uaRNA_three_prime = int(bed.loc[uaRNA].end)
            # Position of uaRNA of - strand gene
            tss_loc = int(tss_bed.loc[gene].end)
            # assert that uaRNA ends after tss
            true_ua = uaRNA_three_prime > tss_loc


        else:
            # 3' end of uaRNA on - strand
            uaRNA_three_prime = int(bed.loc[uaRNA].start)
            # Position of uaRNA of + strand gene
            tss_loc = int(tss_bed.loc[gene].start)
            # assert that uaRNA ends after tss (reverse strand)
            true_ua = uaRNA_three_prime < tss_loc

        if (not true_ua) and verbose:
            print(strand, uaRNA_three_prime, tss_loc)

    else:

        
        if strand == '+':
            # 5' end of gene on + strand
            positive_five_prime = int(bed.loc[uaRNA].start)

            # 3' end of gene on - strand
            negative_three_prime = int(bed.loc[gene].start)

            # 3' end of gene on + strand
            positive_three_prime = int(bed.loc[uaRNA].end)

            # 5' end of gene on - strand
            negative_five_prime = int(bed.loc[gene].end)

        else:
            positive_five_prime = int(bed.loc[gene].start)
            negative_three_prime = int(bed.loc[uaRNA].start)

            positive_three_prime = int(bed.loc[gene].end)
            negative_five_prime = int(bed.loc[uaRNA].end)

        # 5' end of + gene must be after 3' end of - gene. Otherwise 3' gene is contained or colliding
        true_ua = int(positive_five_prime) > int(negative_three_prime)

        # 3' end of + gene must be after 5' end of - gene. Otherwise 5' gene is contained or colliding
        true_ua = true_ua and (int(positive_three_prime) > int(negative_five_prime))

        if (not true_ua) and verbose:
            print(positive_five_prime, positive_three_prime)
            print(negative_three_prime, negative_five_prime)

    return true_ua


def make_pair_dict(pair_list):
    pair_dict = {}
    for pair in pair_list:
        ncRNA, gene = pair.split(':')
        
        if ncRNA in pair_dict.keys():
            pair_dict[ncRNA].append(gene)
        else:
            pair_dict.update({ncRNA:[gene]})
    return pair_dict

def make_tss_dict(ua_pair_list):
    tss_dict = dict({})
    for nc, gene in ua_pair_list:
        if gene not in tss_dict.keys():
            tss_dict.update({gene:[nc]})
        else:
            tss_dict[gene].append(nc)
    return tss_dict
    
def get_closest(allTranscripts, gene, ncRNA_list):
    
    gene_name = gene#.split(':')[0]###########################################################
    strand = allTranscripts.loc[gene_name, 'strand']
    
    closest = ncRNA_list[0]
    closest_distance = 1e10
    
    for nc in ncRNA_list:
        if strand == '-':
            distance = np.abs(int(allTranscripts.loc[gene_name].end) - int(allTranscripts.loc[nc].start))
            if distance < closest_distance:
                closest = nc
                closest_distance = distance
        else:
            distance = np.abs(int(allTranscripts.loc[gene_name].start) - int(allTranscripts.loc[nc].end))
            
    return closest

def get_histone_annotation(histone_file, index_col, annotation):
    histone_idx = pd.read_csv(histone_file, 
                      sep='\t',  names = ['chrom', 'start', 'end', 'name', 'strand'], index_col=index_col).index
    histone_annotation = [1 if x in histone_idx else 0 for x in annotation.index]
    
    return histone_annotation


def main():
    uaRNA_annot = read_bed('NonCodingRNA/annotation/tmp/uaRNA.annotated.bed.gz', '-wo')
    chRPKM = pd.read_csv('RPKM_tables/chRNA.RPKM.bed.gz', sep='\t', index_col=0)
    rpkm_median = chRPKM.median(axis=1)
    annotation = pd.read_csv('NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz', sep='\t', index_col=0)
    
    ncRNA = read_bed('NonCodingRNA/annotation/ncRNA.bed.gz', '-wa')
    ncRNA.index = ncRNA.gene_name

    tss = read_bed('NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz', '-wa')
    tss.index = tss.id

    bed = pd.concat([tss, ncRNA], axis=0)

    print('process uaRNA pairs')
    ua_pair_list = []
    for idx, row in uaRNA_annot.iterrows():
        ua = row.id
        gene = row.id_

        if gene in uaRNA_annot.id:
            score_1 = float(rpkm_median.loc[ua.split(':')[0]])
            score_2 = float(rpkm_median.loc[gene.split(':')[0]])
            #score_2 = float(bed.loc[gene.split(':')[0]].id)###############################

            if score_1 > score_2:
                continue

        if check_true_uaRNA(bed, tss, ua, gene):
            if [ua, gene] not in ua_pair_list:
                ua_pair_list.append([ua, gene])
                
    tss_dict_annot = make_tss_dict(ua_pair_list)
    

    print('filter uaRNAs')
    ua_pair_filtered_annot = dict({})
    for k in tss_dict_annot.keys():
        pc = k#.split(':')[0]  ######################################
        if len(tss_dict_annot[k])>=2:
            ua = get_closest(bed, k, tss_dict_annot[k])
        else:
            ua = tss_dict_annot[k][0]

        if ua in ua_pair_filtered_annot.keys():
            if pc not in ua_pair_filtered_annot[ua]:
                ua_pair_filtered_annot[ua].append(pc)

        else:
            ua_pair_filtered_annot.update({ua:[pc]})
            
    uarnas_annot = pd.Index(ua_pair_filtered_annot.keys())
    
    uarna_list = []
    tss_to_map = []
    
    print('filter tss')

    for nc in uarnas_annot:
        list_of_tss = sorted(ua_pair_filtered_annot[nc])

        dir_of_genes = {}

        for t in list_of_tss:
            tg = t.split(':')[0]
            if tg in dir_of_genes.keys():
                dir_of_genes[tg].append(t)
            else:
                dir_of_genes[tg] = [t]

        list_of_closest = []
        for tg in dir_of_genes.keys():
            closest_tss = get_closest(bed, nc, dir_of_genes[tg])
            list_of_closest.append(closest_tss)

        ua = '|'.join(list_of_closest)#join(ua_pair_filtered[nc])
        uarna_list.append(ua)
        tss_to_map.extend(list_of_closest)
        
    y = [x.split('|') for x in annotation.uaRNA]
    uarnas_new = pd.Index([item for sublist in y for item in sublist])
    
    print('get annotation table')
    
    uarna_annot_df = pd.DataFrame()
    uarna_annot_df['lncRNA_tss'] = uarnas_annot
    uarna_annot_df['uaRNA'] = uarna_list
    
    uarnas_annot_idx = uarnas_annot.difference(pd.Index(uarnas_new))
    uarna_annot_df['lncRNA'] = [x.split(':')[0] for x in uarnas_annot]
    
    uarna_annot_df = uarna_annot_df.loc[uarna_annot_df.lncRNA_tss.isin(uarnas_annot_idx)]
    
    uarna_annotated_df = uarna_annot_df.groupby(['lncRNA'])['uaRNA'].apply(lambda x: '|'.join(x)).reset_index()
    uarna_annotated_df = uarna_annotated_df.set_index('lncRNA')
    
    
    qqnorm_nc = pd.read_csv('QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz', 
                            sep='\t', index_col=4)
    qqnorm_pc = pd.read_csv('QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz', 
                            sep='\t', index_col=4)
    
    uarna_annotated_df = uarna_annotated_df.loc[uarna_annotated_df.index.intersection(
        qqnorm_nc.index.union(qqnorm_pc.index)
    )]
    
    uarna_ = [x.split('|') for x in uarna_annotated_df.uaRNA]
    uarnas_tss = pd.Index([item for sublist in uarna_ for item in sublist])
    
    print('getting SAF')
    
    SAF_ncRNAs=pd.read_csv('NonCodingRNA/annotation/tmp/tss.saf', sep='\t', index_col=0)
    uaRNAs_no_anotados = uarnas_tss.difference(SAF_ncRNAs.index)
    
    tss_saf = bed.loc[uaRNAs_no_anotados, ['id', 'gene_name', 'chrom', 'start', 'end', 'strand']]
    
    GeneID = []
    for idx, row in tss_saf.iterrows():
        if 'ncRNA_' in str(row.gene_name):
            GeneID.append(str(row.gene_name))
        else:
            GeneID.append(str(row.id))
            
    tss_saf['GeneID'] = GeneID
    tss_saf = tss_saf[['GeneID', 'chrom', 'start', 'end', 'strand']]
    
    tss_saf.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
    
    tss_saf_start = [int(x)-200 for x in list((tss_saf.End + tss_saf.Start)/2)]
    tss_saf_end = [int(x)+200 for x in list((tss_saf.End + tss_saf.Start)/2)]
    tss_saf['Start'] = tss_saf_start
    tss_saf['End'] = tss_saf_end
    
    tss_saf = tss_saf.set_index('GeneID')

    tss_all_saf = pd.concat([tss_saf, SAF_ncRNAs], axis=0)
    
    tss_all_saf.to_csv('NonCodingRNA/annotation/tmp/tss_all.saf', sep='\t', index=True, header=True)
    
    uarna_annotated_df.to_csv('NonCodingRNA/annotation/Gencode.uaRNA.annotation.tab.gz', sep='\t',
                                           index=True, header=True)




if __name__ == '__main__':
    
    main()
    
#     uaRNA_annot = read_bed('NonCodingRNA/annotation/tmp/uaRNA.annotated.bed.gz', '-wo')
#     chRPKM = pd.read_csv('RPKM_tables/chRNA.RPKM.bed.gz', sep='\t', index_col=0)
#     rpkm_median = chRPKM.median(axis=1)
#     annotation = pd.read_csv('NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz', sep='\t', index_col=0)
    
#     ncRNA = read_bed('NonCodingRNA/annotation/ncRNA.bed.gz', '-wa')
#     ncRNA.index = ncRNA.gene_name

#     tss = read_bed('NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz', '-wa')
#     tss.index = tss.id

#     bed = pd.concat([tss, ncRNA], axis=0)

#     qqnorm = pd.read_csv('QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz', 
#                   sep='\t', index_col=4)
     
#     print('process uaRNA pairs')
#     ua_pair_list = []
#     for idx, row in uaRNA_annot.iterrows():
#         ua = row.id
#         gene = row.id_

#         if gene in uaRNA_annot.id:
#             score_1 = float(rpkm_median.loc[ua.split(':')[0]])
#             score_2 = float(rpkm_median.loc[gene.split(':')[0]])
#             #score_2 = float(bed.loc[gene.split(':')[0]].id)###############################

#             if score_1 > score_2:
#                 continue

#         if check_true_uaRNA(bed, tss, ua, gene):
#             if [ua, gene] not in ua_pair_list:
#                 ua_pair_list.append([ua, gene])
                
#     tss_dict_annot = make_tss_dict(ua_pair_list)
    

#     print('filter uaRNAs')
#     ua_pair_filtered_annot = dict({})
#     for k in tss_dict_annot.keys():
#         pc = k#.split(':')[0]  ######################################
#         if len(tss_dict_annot[k])>=2:
#             ua = get_closest(bed, k, tss_dict_annot[k])
#         else:
#             ua = tss_dict_annot[k][0]

#         if ua in ua_pair_filtered_annot.keys():
#             if pc not in ua_pair_filtered_annot[ua]:
#                 ua_pair_filtered_annot[ua].append(pc)

#         else:
#             ua_pair_filtered_annot.update({ua:[pc]})
            
#     uarnas_annot = pd.Index(ua_pair_filtered_annot.keys())
    
#     uarna_list = []
#     tss_to_map = []
    
#     print('filter tss')

#     for nc in uarnas_annot:
#         list_of_tss = sorted(ua_pair_filtered_annot[nc])

#         dir_of_genes = {}

#         for t in list_of_tss:
#             tg = t.split(':')[0]
#             if tg in dir_of_genes.keys():
#                 dir_of_genes[tg].append(t)
#             else:
#                 dir_of_genes[tg] = [t]

#         list_of_closest = []
#         for tg in dir_of_genes.keys():
#             closest_tss = get_closest(bed, nc, dir_of_genes[tg])
#             list_of_closest.append(closest_tss)

#         ua = '|'.join(list_of_closest)#join(ua_pair_filtered[nc])
#         uarna_list.append(ua)
#         tss_to_map.extend(list_of_closest)
        
#     y = [x.split('|') for x in annotation.uaRNA]
#     uarnas_new = pd.Index([item for sublist in y for item in sublist])
    
#     print('get annotation table')
    
#     uarna_annot_df = pd.DataFrame()
#     uarna_annot_df['lncRNA_tss'] = uarnas_annot
#     uarna_annot_df['uaRNA'] = uarna_list
    
#     uarnas_annot_idx = uarnas_annot.difference(pd.Index(uarnas_new))
#     uarna_annot_df['lncRNA'] = [x.split(':')[0] for x in uarnas_annot]
    
#     uarna_annot_df = uarna_annot_df.loc[uarna_annot_df.lncRNA_tss.isin(uarnas_annot_idx)]
    
#     uarna_annotated_df = uarna_annot_df.groupby(['lncRNA'])['uaRNA'].apply(lambda x: '|'.join(x)).reset_index()
#     uarna_annotated_df = uarna_annotated_df.set_index('lncRNA')
    
#     uarna_annotated_df = uarna_annotated_df.loc[uarna_annotated_df.index.intersection(qqnorm.index)]
    
#     uarna_ = [x.split('|') for x in uarna_annotated_df.uaRNA]
#     uarnas_tss = pd.Index([item for sublist in uarna_ for item in sublist])
    
#     print('getting SAF')
    
#     SAF_ncRNAs=pd.read_csv('NonCodingRNA/annotation/tmp/tss.saf', sep='\t', index_col=0)
#     uaRNAs_no_anotados = uarnas_tss.difference(SAF_ncRNAs.index)
    
    

#     tss_saf = bed.loc[uaRNAs_no_anotados, ['id', 'gene_name', 'chrom', 'start', 'end', 'strand']]
    
#     GeneID = []
#     for idx, row in tss_saf.iterrows():
#         if 'ncRNA_' in str(row.gene_name):
#             GeneID.append(str(row.gene_name))
#         else:
#             GeneID.append(str(row.id))
            
#     tss_saf['GeneID'] = GeneID
#     tss_saf = tss_saf[['GeneID', 'chrom', 'start', 'end', 'strand']]
    
#     tss_saf.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
    
#     tss_saf_start = [int(x)-200 for x in list((tss_saf.End + tss_saf.Start)/2)]
#     tss_saf_end = [int(x)+200 for x in list((tss_saf.End + tss_saf.Start)/2)]
#     tss_saf['Start'] = tss_saf_start
#     tss_saf['End'] = tss_saf_end
    
#     tss_saf = tss_saf.set_index('GeneID')

#     tss_all_saf = pd.concat([tss_saf, SAF_ncRNAs], axis=0)
    
#     tss_all_saf.to_csv('NonCodingRNA/annotation/tmp/tss_all.saf', sep='\t', index=True, header=True)
    
    
#     uarna_annotated_df.to_csv('NonCodingRNA/annotation/Gencode.uaRNA.annotation.tab.gz', sep='\t',
#                                            index=True, header=True)


    
