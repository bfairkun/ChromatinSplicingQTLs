import numpy as np
import pandas as pd


def get_distance(start, start_gene):
    return abs(int(start) - int(start_gene))



def get_closest_gene(ua, sp_gene_list, ncRNA_bed, pcQTLs_permutation, ncQTLs_permutation, 
                     ncRNA_match_bed, tss_bed):
    
    strand = ncRNA_bed.loc[ua].strand
    if strand == '+':
        uaGene_strand = '-'
    else:
        uaGene_strand = '+'
    chrom = ncRNA_bed.loc[ua].chrom
    
    closest_distance = 1e100
    
    for sp in sp_gene_list:
        
        gene = sp.split(':')[0]
        
        if strand == '+':
            ncRNA_start = int(ncRNA_bed.loc[ua].start)
#             if gene in pcQTLs_permutation.index:
#                 gene_start = int(pcQTLs_permutation.loc[gene].phe_to)
#             elif gene in ncQTLs_permutation.index:
#                 gene_start = int(ncQTLs_permutation.loc[gene].phe_to)
            if sp in tss.index:
                gene_start = int(tss.loc[sp].end)
            elif sp in ncRNA_match_bed.index:
                gene_start = int(ncRNA_match_bed.loc[sp].end)
            else:
                raise Exception(sp)

        else:
            ncRNA_start = int(ncRNA_bed.loc[ua].end)
#             if gene in pcQTLs_permutation.index:
#                 gene_start = int(pcQTLs_permutation.loc[gene].phe_from)
#             elif gene in ncQTLs_permutation.index:
#                 gene_start = int(ncQTLs_permutation.loc[gene].phe_from)
            if sp in tss.index:
                gene_start = int(tss.loc[sp].start)
            elif sp in ncRNA_match_bed.index:
                gene_start = int(ncRNA_match_bed.loc[sp].start)
            else:
                raise Exception(sp)
            
        current_distance = np.abs(ncRNA_start - gene_start)
        if current_distance < closest_distance:
            closest_distance = current_distance
            closest_gene = gene
            closest_sp = sp
            closest_gene_start = gene_start
            closest_ncRNA_start = ncRNA_start
            
        gene_tss = closest_gene_start
        ncRNA_tss = closest_ncRNA_start
        
    return gene, gene_tss, ncRNA_tss, chrom, uaGene_strand


def is_ncRNA(x):
    return x[:5]=='ncRNA'


def get_gene(df, ua, ncQTLs_permutation, pcQTLs_permutation):
    if list(df.index).count(ua) > 1:
        
        if any([is_ncRNA(x) for x in df.loc[ua].gene]):
            tss = next(x for x in df.loc[ua].gene if is_ncRNA(x))
        else:
            tss = get_closest(ncQTLs_permutation, pcQTLs_permutation, ua, 
                                      list(df.loc[ua].gene))
    else:
        tss = df.loc[ua].gene
        
    gene_to_add = tss.split(':')[0]
    
    return gene_to_add, tss


def get_uaRNA_bed(df, ua, ncQTLs_permutation, pcQTLs_permutation,
                 ncQTLs, pcQTLs):
    gene_to_add, tss = get_gene(df, ua, ncQTLs_permutation, pcQTLs_permutation)
    
    gene_list.append(gene_to_add)
    
    if gene_to_add in pcQTLs_permutation.index:
        xQTLs = pcQTLs_permutation
    else:
        xQTLs = ncQTLs_permutation
    
    strand = xQTLs.loc[gene_to_add].phe_strd 
    chrom = xQTLs.loc[gene_to_add].phe_chr
    var_id = ncQTLs_permutation.loc[ua].var_id
    
    var_id = 'chr' + var_id
    
    try:
        signed_effect = str(float(df.loc[(df.uaRNA == ua) & (df.gene == tss)].gene_beta))
    except:
        raise Exception(ua, tss)

        
    pc_start = str(xQTLs.loc[gene_to_add].phe_from)
    pc_end = str(xQTLs.loc[gene_to_add].phe_to)
    
    nc_start = str(ncQTLs_permutation.loc[ua].phe_from)
    nc_end = str(ncQTLs_permutation.loc[ua].phe_to)
    beta = float(df.loc[(df.uaRNA == ua) & (df.gene == tss)].beta)
    
    if beta > 0:
        direction = 1
    else:
        direction = -1
    
    return chrom, pc_start, pc_end, nc_start, nc_end, var_id, signed_effect, strand, direction, gene_to_add, tss


def get_gene_and_nc_lists(annotation_, ncQTLs_permutation, pcQTLs_permutation, bed, ncRNA_match_bed, tss_bed):

    genes_list = []
    nc_list = []
    gene_tss_list = []
    ncRNA_tss_list = []
    gene_tts_list = []
    ncRNA_tts_list = []
    chrom_list = []
    uaGene_strand_list = []
    for ua in annotation_.index:
        print(ua)
        splitted = annotation_.loc[ua].uaRNA.split('|')
        gene_id_list = []
        nc_id_list = []
        sp_gene_list = []
        sp_nc_list = []
        for sp in splitted:
            gene_name = sp.split(':')[0]
            if gene_name in pcQTLs_permutation.index:
                gene_id_list.append(gene_name)
                sp_gene_list.append(sp)
            elif gene_name in ncQTLs_permutation.index:
                nc_id_list.append(gene_name)
                sp_nc_list.append(sp)

        if len(gene_id_list) >= 1:
            gene_ids, gene_tss, ncRNA_tss, chrom, uaGene_strand = get_closest_gene(
                ua, sp_gene_list, bed, pcQTLs_permutation, 
                                        ncQTLs_permutation, ncRNA_match_bed, tss_bed)
            
#             gene_tss = pcQTLs_permutation.loc[gene_ids].phe_from
#             ncRNA_tss = ncQTLs_permutation.loc[ua].phe_from
#             gene_tts = pcQTLs_permutation.loc[gene_ids].phe_to
#             ncRNA_tts = ncQTLs_permutation.loc[ua].phe_to
    
            ncRNA_tss = ncQTLs_permutation.loc[ua].phe_from
            ncRNA_tts = ncQTLs_permutation.loc[ua].phe_to
            if uaGene_strand == '+':
                gene_tts = pcQTLs_permutation.loc[gene_ids].phe_to
                
            else:
                gene_tts = gene_tss
                gene_tss = pcQTLs_permutation.loc[gene_ids].phe_from
                
            if int(gene_tts) < int(gene_tss):
                print('Gene TTS < TSS. This should be rare')
                gene_tss = pcQTLs_permutation.loc[gene_ids].phe_from
                gene_tts = pcQTLs_permutation.loc[gene_ids].phe_to
                
                
            nc_ids = '.'
            
        elif len(nc_id_list) >= 1:
            nc_ids, gene_tss, ncRNA_tss, chrom, uaGene_strand = get_closest_gene(
                ua, sp_nc_list, bed, ncQTLs_permutation, 
                                      ncQTLs_permutation, ncRNA_match_bed, tss_bed)
            gene_ids = '.'  
            if uaGene_strand == '+':
                gene_tts = ncQTLs_permutation.loc[nc_ids].phe_to
                
            else:
                gene_tts = gene_tss
                gene_tss = ncQTLs_permutation.loc[nc_ids].phe_from
                
            if int(gene_tts) < int(gene_tss):
                print('ncRNA TTS < TSS. This should be rare')
                gene_tss = ncQTLs_permutation.loc[nc_ids].phe_from
                gene_tts = ncQTLs_permutation.loc[nc_ids].phe_to
#             if uaGene_strand == '+':
#                 gene_tts = ncQTLs_permutation.loc[nc_ids].phe_to
#                 ncRNA_tts = ncQTLs_permutation.loc[ua].phe_from
#             else:
#             gene_tss = ncQTLs_permutation.loc[nc_ids].phe_from
            ncRNA_tss = ncQTLs_permutation.loc[ua].phe_from
#             gene_tts = ncQTLs_permutation.loc[nc_ids].phe_to
            ncRNA_tts = ncQTLs_permutation.loc[ua].phe_to
        else:
            nc_ids = '.'
            gene_ids = '.'
            gene_tss = '.'
            ncRNA_tss = '.'
            gene_tts = '.'
            ncRNA_tts = '.'
            chrom = '.'
            uaGene_strand = '.'
            
        genes_list.append(gene_ids)
        nc_list.append(nc_ids)
#         if uaGene_strand == '+':
        gene_tss_list.append(gene_tss)
        ncRNA_tss_list.append(ncRNA_tss)
        gene_tts_list.append(gene_tts)
        ncRNA_tts_list.append(ncRNA_tts)
#         else:
#             gene_tss_list.append(gene_tts)
#             ncRNA_tss_list.append(ncRNA_tss)
#             gene_tts_list.append(gene_tss)
#             ncRNA_tts_list.append(ncRNA_tts)
        chrom_list.append(chrom)
        uaGene_strand_list.append(uaGene_strand)
        
    return genes_list, nc_list, gene_tss_list, ncRNA_tss_list, gene_tts_list, ncRNA_tts_list, chrom_list, uaGene_strand_list


def read_bed(bed_file, flag='-wa'):
    if flag == '-wa':
        bed = pd.read_csv(bed_file, sep='\t',
                    names = ['chrom', 'start', 'end', 'gene_name', 'id', 'strand'])
    elif flag == '-wo':
        bed = pd.read_csv(bed_file, sep='\t', 
                         names = ['chrom', 'start', 'end', 'gene_name', 'id', 'strand',
                                  'chrom_', 'start_', 'end_', 'gene_name_', 'id_', 'strand_','overlap'])
        
    return bed



if __name__ == '__main__':
    
    print('loading data')

    annotation = pd.read_csv('NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz', sep='\t', index_col=0)
    gencode_annotation = pd.read_csv('NonCodingRNA/annotation/Gencode.uaRNA.annotation.tab.gz', 
                                     index_col=0,  sep='\t')

    ncRNA_bed = pd.read_csv('NonCodingRNA/annotation/ncRNA.bed.gz', sep='\t', index_col=3,
                           names = ['chrom', 'start', 'end', 'score', 'strand'])

    ncQTLs = pd.read_csv('QTLs/QTLTools/chRNA.Expression_ncRNA/NominalPass.txt.gz', sep=' ')
    pcQTLs = pd.read_csv('QTLs/QTLTools/chRNA.Expression.Splicing/NominalPass.txt.gz', sep=' ')

    ncQTLs_permutation = pd.read_csv('QTLs/QTLTools/chRNA.Expression_ncRNA/PermutationPass.FDR_Added.txt.gz', sep=' ')
    ncQTLs_permutation.index = ncQTLs_permutation.phe_id
    sig_ncRNAs = ncQTLs_permutation.loc[ncQTLs_permutation.q <= 0.1]
    pcQTLs_permutation = pd.read_csv('QTLs/QTLTools/chRNA.Expression.Splicing/PermutationPass.FDR_Added.txt.gz', sep=' ')
    pcQTLs_permutation.index = pcQTLs_permutation.phe_id

    uaRNAs = annotation.loc[annotation.uaRNA != '.'].index
    incRNAs = annotation.loc[annotation.rna_type == 'incRNA'].index
    rtRNAs = annotation.index.difference(incRNAs.union(uaRNAs))


    ncRNA = read_bed('NonCodingRNA/annotation/ncRNA.bed.gz', '-wa')
    ncRNA.index = ncRNA.gene_name

    tss = read_bed('NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz', '-wa')
    tss.index = tss.id

    tss_bed = pd.concat([tss, ncRNA], axis=0)

    ncRNA_bed = pd.read_csv('NonCodingRNA/annotation/NonCodingRNA.bed.gz', sep='\t', index_col=3,
                           names = ['chrom', 'start', 'end', 'score', 'strand'])

    ncQTLs.index = ncQTLs.phe_id + '|' + ncQTLs.var_id
    pcQTLs.index = pcQTLs.phe_id + '|' + pcQTLs.var_id


    all_bed = pd.read_csv('NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz', sep='\t', index_col=3,
                           names = ['chrom', 'start', 'end', 'score', 'strand'])

    rpkm = pd.read_csv('RPKM_tables/chRNA.RPKM.bed.gz', sep='\t', index_col=0)


    print('processing hmm annotation')

    genes_list, nc_list, gene_tss_list, ncRNA_tss_list, gene_tts_list, ncRNA_tts_list, chrom_list, uaGene_strand_list = get_gene_and_nc_lists(
        annotation.loc[annotation.index.intersection(ncRNA_bed.index)], 
                          ncQTLs_permutation, pcQTLs_permutation, ncRNA_bed, ncRNA_bed, tss_bed)

    annotation_slice = annotation.loc[annotation.index.intersection(ncRNA_bed.index)].copy()
    annotation_slice['chrom'] = chrom_list
    annotation_slice['uaGene_strand'] = uaGene_strand_list
    annotation_slice['protein_coding'] = genes_list
    annotation_slice['ncRNA'] = nc_list
    
    annotation_slice['ncRNA_start'] = ncRNA_tss_list
    annotation_slice['ncRNA_end'] = ncRNA_tts_list
    annotation_slice['uaGene_start'] = gene_tss_list
    annotation_slice['uaGene_end'] = gene_tts_list
    

    annotation_slice = annotation_slice.loc[annotation_slice.index.intersection(ncQTLs_permutation.index)]

    annotation_uaRNA = annotation_slice.loc[
        (annotation_slice.protein_coding != '.') | (annotation_slice.ncRNA != '.')
    ]


    print('processing Gencode annotation')
    
    gencode_annotation = gencode_annotation.loc[
        gencode_annotation.index.intersection(ncQTLs_permutation.index)
    ]

    
    genes_list, nc_list, gene_tss_list, ncRNA_tss_list, gene_tts_list, ncRNA_tts_list, chrom_list, uaGene_strand_list = get_gene_and_nc_lists(
        gencode_annotation, ncQTLs_permutation, pcQTLs_permutation, all_bed, ncRNA_bed, tss_bed)

    gencode_annotation['chrom'] = chrom_list
    gencode_annotation['uaGene_strand'] = uaGene_strand_list
    gencode_annotation['protein_coding'] = genes_list
    gencode_annotation['ncRNA'] = nc_list
    
    gencode_annotation['ncRNA_start'] = ncRNA_tss_list
    gencode_annotation['ncRNA_end'] = ncRNA_tts_list
    gencode_annotation['uaGene_start'] = gene_tss_list
    gencode_annotation['uaGene_end'] = gene_tts_list
    
    

    gencode_uaRNA = gencode_annotation.loc[
        (gencode_annotation.protein_coding != '.') | (gencode_annotation.ncRNA != '.')
    ]

    annotation_pc_hmm = pd.DataFrame(annotation_uaRNA.loc[
        (annotation_uaRNA.protein_coding.isin(pcQTLs_permutation.index)), 
                                                 ['chrom', 'uaGene_strand', 
                                                  'protein_coding', 'ncRNA_start', 'ncRNA_end',
                                                  'uaGene_start', 'uaGene_end'
                                                 ]])

    annotation_pc_hmm['uaRNA'] = annotation_pc_hmm.index
    annotation_pc_hmm = annotation_pc_hmm.set_index('protein_coding')

    annotation_pc_lncRNA = pd.DataFrame(gencode_uaRNA.loc[
        (gencode_uaRNA.protein_coding.isin(pcQTLs_permutation.index)), 
                                        ['chrom', 'uaGene_strand', 
                                         'protein_coding', 'ncRNA_start', 'ncRNA_end',
                                                  'uaGene_start', 'uaGene_end'
                                        ]])

    annotation_pc_lncRNA['uaRNA'] = annotation_pc_lncRNA.index
    annotation_pc_lncRNA = annotation_pc_lncRNA.set_index('protein_coding')


    uaRNAs_pc = pd.concat([annotation_pc_hmm, 
               annotation_pc_lncRNA.loc[annotation_pc_lncRNA.index.difference(annotation_pc_hmm.index)]])

    uaRNAs_pc['nc_RPKM'] = list(rpkm.loc[uaRNAs_pc.uaRNA].median(axis=1))
    uaRNAs_pc['pc_RPKM'] = list(rpkm.loc[uaRNAs_pc.index].median(axis=1))



    annotation_nc_hmm = pd.DataFrame(annotation_uaRNA.loc[
        (annotation_uaRNA.ncRNA.isin(ncQTLs_permutation.index)), 
                                                 ['chrom', 'uaGene_strand', 
                                                  'ncRNA', 'ncRNA_start', 'ncRNA_end',
                                                  'uaGene_start', 'uaGene_end']])

    annotation_nc_hmm['uaRNA'] = annotation_nc_hmm.index
    annotation_nc_hmm = annotation_nc_hmm.set_index('ncRNA')

    annotation_nc_lncRNA = pd.DataFrame(gencode_uaRNA.loc[
        (gencode_uaRNA.ncRNA.isin(ncQTLs_permutation.index)), 
                                         ['chrom', 'uaGene_strand', 
                                          'ncRNA', 'ncRNA_start', 'ncRNA_end',
                                                  'uaGene_start', 'uaGene_end']])

    annotation_nc_lncRNA['uaRNA'] = annotation_nc_lncRNA.index
    annotation_nc_lncRNA = annotation_nc_lncRNA.set_index('ncRNA')


    uaRNAs_nc = pd.concat([annotation_nc_hmm, 
               annotation_nc_lncRNA.loc[annotation_nc_lncRNA.index.difference(annotation_nc_hmm.index)]])

    uaRNAs_nc['nc_RPKM'] = list(rpkm.loc[uaRNAs_nc.uaRNA].median(axis=1))
    uaRNAs_nc['pc_RPKM'] = list(rpkm.loc[uaRNAs_nc.index].median(axis=1))

    df = pd.concat([uaRNAs_pc, uaRNAs_nc], axis=0)
    df['uaGene'] = df.index
    df = df.reset_index(drop=True)




    row_list = []
    for idx, row in df.iterrows():
        uaRNA = row.uaRNA
        uaGene = row.uaGene

        if uaRNA in ncQTLs_permutation.index:
            ncQTL_var = ncQTLs_permutation.loc[uaRNA].var_id
            ncQTL_var_slope = ncQTLs_permutation.loc[uaRNA].slope
            ncQTL_nom_pval = ncQTLs_permutation.loc[uaRNA].nom_pval
            ncQTL_adj_beta_pval = ncQTLs_permutation.loc[uaRNA].adj_beta_pval
            ncQTL_q = ncQTLs_permutation.loc[uaRNA].q
        else:
            raise Exception('missing uaRNA ncQTL')

        if uaGene in pcQTLs_permutation.index:
            eQTL_var = pcQTLs_permutation.loc[uaGene].var_id
            eQTL_var_slope = pcQTLs_permutation.loc[uaGene].slope
            eQTL_nom_pval = pcQTLs_permutation.loc[uaGene].nom_pval
            eQTL_adj_beta_pval = pcQTLs_permutation.loc[uaGene].adj_beta_pval
            eQTL_q = pcQTLs_permutation.loc[uaGene].q
            use_pc = True
        elif uaGene in ncQTLs_permutation.index:
            eQTL_var= ncQTLs_permutation.loc[uaGene].var_id
            eQTL_var_slope = ncQTLs_permutation.loc[uaGene].slope
            eQTL_nom_pval = ncQTLs_permutation.loc[uaGene].nom_pval
            eQTL_adj_beta_pval = ncQTLs_permutation.loc[uaGene].adj_beta_pval
            eQTL_q = ncQTLs_permutation.loc[uaGene].q
            use_pc = False
        else:
            raise Exception('missing uaGene eQTL or ncQTL')

        ncQTL_idx = uaRNA + '|' + eQTL_var
        eQTL_idx = uaGene + '|' + ncQTL_var

        row_out = [ncQTL_var, ncQTL_var_slope, ncQTL_nom_pval, ncQTL_adj_beta_pval, ncQTL_q, 
               eQTL_var, eQTL_var_slope, eQTL_nom_pval, eQTL_adj_beta_pval, eQTL_q, 
               ncQTL_idx, eQTL_idx
              ]
        row_list.append(row_out)


    df_extension = pd.DataFrame(row_list,
                                columns = ['ncQTL_var', 'ncQTL_slope', 'ncQTL_nom_pval', 
                                         'ncQTL_adj_beta_pval', 'ncQTL_q', 
                                         'eQTL_var', 'eQTL_slope', 'eQTL_nom_pval', 'eQTL_adj_beta_pval', 
                                         'eQTL_q', 'ncQTL_idx', 'eQTL_idx'])

    df_extension['uaRNA_eQTL_nom_pval'] = list(ncQTLs.reindex(pd.Index(df_extension.ncQTL_idx)).nom_pval)
    df_extension['uaRNA_eQTL_slope'] = list(ncQTLs.reindex(pd.Index(df_extension.ncQTL_idx)).slope)

    uaGene_ncQTL_nom_pval = []
    uaGene_ncQTL_slope = []

    for idx in pd.Index(df_extension.eQTL_idx):
        if idx in pcQTLs.index:
            uaGene_ncQTL_nom_pval.append(float(pcQTLs.loc[idx].nom_pval))
            uaGene_ncQTL_slope.append(float(pcQTLs.loc[idx].slope))
        elif idx in ncQTLs.index:
            uaGene_ncQTL_nom_pval.append(float(ncQTLs.loc[idx].nom_pval))
            uaGene_ncQTL_slope.append(float(ncQTLs.loc[idx].slope))
        else:

            uaGene_ncQTL_nom_pval.append(np.nan)
            uaGene_ncQTL_slope.append(np.nan)


    df_extension['uaGene_ncQTL_nom_pval'] = uaGene_ncQTL_nom_pval
    df_extension['uaGene_ncQTL_slope'] = uaGene_ncQTL_slope


    uaRNA_df = pd.merge(df, df_extension, left_index=True, right_index=True)

    uaRNA_df['ncQTL_beta'] = uaRNA_df.ncQTL_slope * uaRNA_df.uaGene_ncQTL_slope
    uaRNA_df['eQTL_beta'] = uaRNA_df.eQTL_slope * uaRNA_df.uaRNA_eQTL_slope

    uaRNA_df['ncQTL_sign'] = [1 if (x > 0) else -1 for x in uaRNA_df.ncQTL_beta]
    uaRNA_df['eQTL_sign'] = [1 if (x > 0) else -1 for x in uaRNA_df.eQTL_beta]

    uaRNA_df['ncQTL_significance'] = [np.mean([-np.log10(uaRNA_df.loc[i].ncQTL_nom_pval),
                                              -np.log10(uaRNA_df.loc[i].uaGene_ncQTL_nom_pval)]
                                             ) for i in uaRNA_df.index]

    uaRNA_df['eQTL_significance'] = [np.mean([-np.log10(uaRNA_df.loc[i].eQTL_nom_pval),
                                              -np.log10(uaRNA_df.loc[i].uaRNA_eQTL_nom_pval)]
                                             ) for i in uaRNA_df.index]

    uaRNA_df['ncQTL_color'] = (uaRNA_df.ncQTL_sign * uaRNA_df.ncQTL_significance)/np.max(np.abs((uaRNA_df.ncQTL_sign * uaRNA_df.ncQTL_significance)))
    uaRNA_df['eQTL_color'] = (uaRNA_df.eQTL_sign * uaRNA_df.eQTL_significance)/np.max(np.abs((uaRNA_df.eQTL_sign * uaRNA_df.eQTL_significance)))

    uaRNA_df = uaRNA_df.set_index('uaRNA')
    
    
    uaRNA_df.to_csv('NonCodingRNA/QTLs/summary.uaRNA.tab.gz', sep='\t', index=True, header=True)






