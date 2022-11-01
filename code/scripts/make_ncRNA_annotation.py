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

if __name__ == '__main__':
    
    print("loading files")
    
    uaRNA = read_bed('NonCodingRNA/annotation/tmp/uaRNA.bed.gz', '-wo')
    
    allGenes = read_bed('NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz', '-wa')
    allGenes.index = allGenes.gene_name
    ncRNA = read_bed('NonCodingRNA/annotation/ncRNA.bed.gz', '-wa')
    ncRNA.index = ncRNA.gene_name
    lncrnas = read_bed("NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz", '-wo')
    
    print(lncrnas.gene_name)
    incRNA = read_bed('NonCodingRNA/annotation/tmp/incRNA.bed.gz', '-wa')

#    bed = pd.concat([allGenes, ncRNA], axis=0)##############################
    
    tss = read_bed('NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz', '-wa')
    tss.index = tss.id
 
    bed = pd.concat([tss, ncRNA], axis=0)
   
    rtRNA = read_bed('NonCodingRNA/annotation/tmp/rtRNA.bed.gz', '-wo')
    srtRNA = read_bed('NonCodingRNA/annotation/tmp/srtRNA.bed.gz', '-wo')

    
    
    ua_pair_list = []
    for idx, row in uaRNA.iterrows():
        ua = row.id
        gene = row.id_

        if gene[:5] == 'ncRNA':
            score_1 = float(bed.loc[ua].id)
            score_2 = float(bed.loc[gene].id)
            #score_2 = float(bed.loc[gene.split(':')[0]].id)###############################

            if score_1 > score_2:
                continue

        if check_true_uaRNA(bed, tss, ua, gene):
            if [ua, gene] not in ua_pair_list:
                ua_pair_list.append([ua, gene])

    tss_dict = make_tss_dict(ua_pair_list)
    
    
    
    ua_pair_filtered = dict({})
    for k in tss_dict.keys():
        pc = k#.split(':')[0]  ######################################
        if len(tss_dict[k])>=2:
            ua = get_closest(bed, k, tss_dict[k])
        else:
            ua = tss_dict[k][0]

        if ua in ua_pair_filtered.keys():
            if pc not in ua_pair_filtered[ua]:
                ua_pair_filtered[ua].append(pc)

        else:
            ua_pair_filtered.update({ua:[pc]})
            
    ncrnas = ncRNA.index
    uarnas = pd.Index(ua_pair_filtered.keys())
    incrnas = pd.Index(incRNA.gene_name).difference(uarnas)
    srtrnas = pd.Index(srtRNA.gene_name).difference(uarnas.union(incrnas))
    rtrnas = ncrnas.difference(uarnas).difference(incrnas).difference(srtrnas)
    
    idx = list(uarnas) + list(incrnas) + list(srtrnas) + list(rtrnas)
    
    annotation = pd.DataFrame(index=idx)
    
    
    rna_type = (['uaRNA']*len(uarnas)) + (['incRNA']*len(incrnas)) + (['rtRNA']*(len(srtrnas)+len(rtrnas)))
    annotation['rna_type'] = rna_type
    
    all_reverse = pd.concat([srtRNA, rtRNA], axis=0, ignore_index=True)
    
    uarna_list = []
    rtrna_list = []
    lncrna_list = []

    tss_to_map = []
    
    for nc in annotation.index:
        if nc in ua_pair_filtered.keys():
            print('uaRNA')
            list_of_tss = sorted(ua_pair_filtered[nc])
            
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
        else:
            uarna_list.append('.')
            
        if nc in list(all_reverse.gene_name):
            print('reverse')
            rev = '|'.join(all_reverse.loc[all_reverse.gene_name == nc].gene_name_.unique())
            rtrna_list.append(rev)
        else:
            rtrna_list.append('.')
            
        if nc in list(lncrnas.gene_name):
            print('lncRNA')
            rev = '|'.join(lncrnas.loc[lncrnas.gene_name == nc].gene_name_.unique())
            lncrna_list.append(rev)
        else:
            lncrna_list.append('.')
            
    annotation['uaRNA'] = uarna_list
    annotation['rtRNA'] = rtrna_list
    annotation['lncRNA'] = lncrna_list

    tss_to_map_filtered = sorted(set(tss_to_map))
    tss_saf = bed.loc[tss_to_map, ['id', 'chrom', 'start', 'end', 'strand']]
    tss_saf.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand']
    
    tss_saf_start = [int(x)-200 for x in list((tss_saf.End + tss_saf.Start)/2)]
    tss_saf_end = [int(x)+200 for x in list((tss_saf.End + tss_saf.Start)/2)]
    tss_saf['Start'] = tss_saf_start
    tss_saf['End'] = tss_saf_end
    
    tss_saf.to_csv('NonCodingRNA/annotation/tmp/tss.saf', sep='\t', index=False, header=True)
    
    k4me1 = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K4ME1.overlaps.bed.gz',
                                   3, annotation)
    k4me3 = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K4ME3.overlaps.bed.gz', 
                                   3, annotation)
    k27ac = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K27AC.overlaps.bed.gz', 
                                   3, annotation)

    annotation['H3K4ME1'] = k4me1
    annotation['H3K4ME3'] = k4me3
    annotation['H3K27AC'] = k27ac
    
    k4me1 = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K4ME1.TSS_overlaps.bed.gz', 
                                   3, annotation)
    k4me3 = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K4ME3.TSS_overlaps.bed.gz', 
                                   3, annotation)
    k27ac = get_histone_annotation('NonCodingRNA/annotation/histone_marks/ncRNA.H3K27AC.TSS_overlaps.bed.gz',
                                   3,  annotation)

    annotation['H3K4ME1_TSS'] = k4me1
    annotation['H3K4ME3_TSS'] = k4me3
    annotation['H3K27AC_TSS'] = k27ac
    
    
    
    
    marks = []
    for idx, row in annotation.iterrows():
        if row.H3K4ME3_TSS + row.H3K27AC_TSS == 2:
            marks.append('promoter')
        elif (row.H3K4ME1_TSS == 1) and (row.H3K27AC_TSS == 1):
            marks.append('enhancer')
        elif row[['H3K4ME1_TSS', 'H3K4ME3_TSS', 'H3K27AC_TSS']].sum() > 0:
            marks.append('other')
        else:
            marks.append('no mark')

    annotation['histone_marks'] = marks
    
    annotation['RPKM'] = ncRNA.loc[annotation.index].id
    
    ncRNA_filtered = ncRNA.loc[annotation.loc[(annotation.histone_marks != 'no mark') | (annotation.RPKM >= 10)].index]
    
#     ncRNA_filtered = ncRNA.loc[annotation.loc[(annotation.histone_marks != 'no mark') | (annotation.RPKM >= 10)].index]
    
    ncRNA_filtered.to_csv('NonCodingRNA/annotation/NonCodingRNA.bed.gz', sep='\t',
                          index=False, header=False)
    
    
    
    annotation.to_csv('NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz', sep='\t', index=True, header = True)
    
    
    allgenes = pd.read_csv('NonCodingRNA/annotation/allGenes.TSS.bed.gz', sep='\t',  
                           index_col=4, names = ['chrom', 'start', 'end', 'name', 'strand'])
    
    k4me1_genes = get_histone_annotation('NonCodingRNA/annotation/histone_marks/allGenes.H3K4ME1.overlaps.bed.gz', 4,
                                   allgenes)

    k4me3_genes = get_histone_annotation('NonCodingRNA/annotation/histone_marks/allGenes.H3K4ME3.overlaps.bed.gz', 4,
                                   allgenes)

    k27ac_genes = get_histone_annotation('NonCodingRNA/annotation/histone_marks/allGenes.H3K27AC.overlaps.bed.gz', 4,
                                   allgenes)

    allgenes['H3K4ME1'] = k4me1_genes
    allgenes['H3K4ME3'] = k4me3_genes
    allgenes['H3K27AC'] = k27ac_genes
    
    k4me1_genes = get_histone_annotation(
        'NonCodingRNA/annotation/histone_marks/allGenes.H3K4ME1.TSS_overlaps.bed.gz', 
        4, allgenes)

    k4me3_genes = get_histone_annotation(
        'NonCodingRNA/annotation/histone_marks/allGenes.H3K4ME3.TSS_overlaps.bed.gz', 
        4, allgenes)

    k27ac_genes = get_histone_annotation(
        'NonCodingRNA/annotation/histone_marks/allGenes.H3K27AC.TSS_overlaps.bed.gz', 
        4, allgenes)

    allgenes['H3K4ME1_TSS'] = k4me1_genes
    allgenes['H3K4ME3_TSS'] = k4me3_genes
    allgenes['H3K27AC_TSS'] = k27ac_genes

    allgenes_out = allgenes[['name', 'H3K4ME1', 'H3K4ME3', 'H3K27AC', 'H3K4ME1_TSS', 'H3K4ME3_TSS', 'H3K27AC_TSS']]
    
    allgenes_out.to_csv('NonCodingRNA/annotation/allGenes.annotation.tab.gz', sep='\t', index=True, header = True)






    
    
    
    
    
