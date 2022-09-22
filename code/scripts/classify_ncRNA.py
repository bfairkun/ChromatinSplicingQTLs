import numpy as np
import pandas as pd
import argparse
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

def check_true_uaRNA(bed, ncRNA1, ncRNA2):

    strand1 = bed.loc[ncRNA1].strand
    if strand1 == '+':
        try:
            pos_five = int(bed.loc[ncRNA1].start)
        except:
            pos_five = int(bed.loc[ncRNA1].start.iloc[0])
            print(pos_five)
        try:
            neg_three = int(bed.loc[ncRNA2].start)
        except:
            neg_three = int(bed.loc[ncRNA2].start.iloc[0])
            print(neg_three)
    else:
        try:
            pos_five = int(bed.loc[ncRNA2].start)
        except:
            pos_five = int(bed.loc[ncRNA2].start.iloc[0])
            print(pos_five)
        try:
            neg_three = int(bed.loc[ncRNA1].start)
            
        except:
            neg_three = int(bed.loc[ncRNA1].start.iloc[0])
            print(neg_three)
     
    try:
        true_ua = int(pos_five) > int(neg_three)
    except:
        print(pos_five)
        print(neg_three)
    
    return true_ua

def is_uaRNA(allTranscripts, nc_pair, srtRNA, rtRNA_rev, coRNA):
    
    token = (nc_pair in srtRNA)
    token = token or (nc_pair in rtRNA_rev)
    token = token or (nc_pair in coRNA)
    
    if token:
        return False
    
    n1, n2 = nc_pair.split(':')
    is_ua = check_true_uaRNA(allTranscripts, n1, n2)
    
    return is_ua

def filter_uaRNA(allTranscripts, uaRNA_idx, srtRNA_idx, rtRNA_rev_idx, coRNA_idx):
    true_uaRNA = {}
    for ua in uaRNA_idx:
        is_ua = is_uaRNA(allTranscripts, ua, srtRNA_idx, rtRNA_rev_idx, coRNA_idx)
        if is_ua:
            
            ncRNA_id, gene_id = ua.split(':')
            if gene_id in true_uaRNA.keys():
                true_uaRNA[gene_id].append(ncRNA_id)
                
            else:
                true_uaRNA.update({gene_id:[ncRNA_id]})
            
#             true_uaRNA.append(ua)
            
    
    true_uaRNA_list = []
    for g in true_uaRNA.keys():
        nc_list = true_uaRNA[g]
        if len(nc_list) > 1:
            uaRNA_ = get_closest(allTranscripts, g, nc_list)
        else:
            uaRNA_ = true_uaRNA[g][0]
            
        true_uaRNA_list.append(uaRNA_ + ':' + g)
            
    return pd.Index(true_uaRNA_list)

def get_closest(allTranscripts, gene, ncRNA_list):
    strand = allTranscripts.loc[gene, 'strand']
    
    closest = ncRNA_list[0]
    closest_distance = 1e10
    
    for nc in ncRNA_list:
        if strand == '-':
            distance = np.abs(int(allTranscripts.loc[gene].end) - int(allTranscripts.loc[nc].start))
            if distance < closest_distance:
                closest = nc
                closest_distance = distance
        else:
            distance = np.abs(int(allTranscripts.loc[gene].start) - int(allTranscripts.loc[nc].end))
            
    return closest

def make_pair_dict(pair_list):
    pair_dict = {}
    for pair in pair_list:
        ncRNA, gene = pair.split(':')
        
        if ncRNA in pair_dict.keys():
            pair_dict[ncRNA].append(gene)
        else:
            pair_dict.update({ncRNA:[gene]})
    return pair_dict


parser = argparse.ArgumentParser()
parser.add_argument('--merge', action='store_true', required=False)


if __name__ == '__main__':
    
    args = parser.parse_args()
    merge = args.merge
    
    if merge:
        dir_name = "NonCodingRNA_merged"
    else:
        dir_name = "NonCodingRNA_annotation"
    
    print("loading files")
    
    allGenes = read_bed(dir_name + '/annotation/tmp/allGenes.Gencode.bed.gz', '-wa')
    ncRNA = read_bed(dir_name + '/annotation/ncRNA.bed.gz', '-wa')
    incRNA = read_bed(dir_name + '/annotation/tmp/incRNA.bed.gz', '-wa')
    uaRNA = read_bed(dir_name + '/annotation/tmp/uaRNA.bed.gz', '-wo')
    coRNA = read_bed(dir_name + '/annotation/tmp/coRNA.bed.gz', '-wo')
    ctRNA = read_bed(dir_name + '/annotation/tmp/ctRNA.bed.gz', '-wo')
    rtRNA = read_bed(dir_name + '/annotation/tmp/rtRNA.bed.gz', '-wo')
    srtRNA = read_bed(dir_name + '/annotation/tmp/srtRNA.bed.gz', '-wo')
    lncRNA = read_bed(dir_name + '/annotation/tmp/lncRNA.ncRNA.bed.gz', '-wo')
    snoRNA = read_bed(dir_name + '/annotation/tmp/snoRNA.ncRNA.bed.gz', '-wo')
    pseudogenes = read_bed(dir_name + '/annotation/tmp/pseudogenes.ncRNA.bed.gz', '-wo')
    
    print("ncRNA:gene")
    
    uaRNA_gene = uaRNA.loc[[x[:5] != 'ncRNA' for x in uaRNA.gene_name_]]
    uaRNA_ncRNA = uaRNA.loc[[x[:5] == 'ncRNA' for x in uaRNA.gene_name_]]

    uaRNA_idx_gene = pd.Index((uaRNA_gene['gene_name'] + ':' + uaRNA_gene['id_']).unique())
    uaRNA_idx_ncRNA = pd.Index((uaRNA_ncRNA['gene_name'] + ':' + uaRNA_ncRNA['gene_name_']).unique())

    uaRNA_idx = uaRNA_idx_gene.union(uaRNA_idx_ncRNA)
    
    #uaRNA_idx = pd.Index((uaRNA['gene_name'] + ':' + uaRNA['id_']).unique())
    rtRNA_idx = pd.Index((rtRNA['gene_name'] + ':' + rtRNA['gene_name_']).unique())
    srtRNA_idx = pd.Index((srtRNA['gene_name'] + ':' + srtRNA['gene_name_']).unique())
    rtRNA_rev_idx = pd.Index([x.split(':')[1] +':' + x.split(':')[0] for x in srtRNA_idx])
    coRNA_idx = pd.Index((coRNA['gene_name'] + ':' + coRNA['gene_name_']).unique())
    ctRNA_idx = pd.Index((ctRNA['gene_name'] + ':' + ctRNA['gene_name_']).unique())
    
    lncRNA_idx = pd.Index((lncRNA['gene_name'] + ':' + lncRNA['gene_name_']).unique())
    snoRNA_idx = pd.Index((snoRNA['gene_name'] + ':' + snoRNA['gene_name_']).unique())
    pseudo_idx = pd.Index((pseudogenes['gene_name'] + ':' + pseudogenes['gene_name_']).unique())
    
    allGenes.index = allGenes.gene_name
    ncRNA.index = ncRNA.gene_name
    allTranscripts = pd.concat([allGenes, ncRNA], axis=0)
    
    print("Filter classifications")
    
    true_uaRNA = filter_uaRNA(allTranscripts, uaRNA_idx, srtRNA_idx, rtRNA_rev_idx, coRNA_idx)
    true_coRNA = coRNA_idx.difference(srtRNA_idx).difference(rtRNA_rev_idx)
    uaRNAs = pd.Index([x.split(':')[0] for x in true_uaRNA])
    true_incRNA = pd.Index(incRNA.gene_name).difference(uaRNAs)
    # true_coRNA = coRNA_idx.difference(srtRNA_idx).difference(rtRNA_rev_idx)
    
    
    uaRNA_dict = make_pair_dict(true_uaRNA)
    coRNA_dict = make_pair_dict(true_coRNA)
    rtRNA_dict = make_pair_dict(rtRNA_idx.difference(srtRNA_idx))
    srtRNA_dict = make_pair_dict(srtRNA_idx)
    lncRNA_dict = make_pair_dict(lncRNA_idx)
    snoRNA_dict = make_pair_dict(snoRNA_idx)
    pseudo_dict = make_pair_dict(pseudo_idx)
    ctRNA_dict = make_pair_dict(ctRNA_idx)
    
    print("make annotation table")
    
    annotation_df = pd.DataFrame(index = ncRNA.gene_name)
    
    uaRNA_annot = []
    coRNA_annot = []
    rtRNA_annot = []
    srtRNA_annot = []
    lncRNA_annot = []
    snoRNA_annot = []
    pseudo_annot = []
    ctRNA_annot = []

    annot = []
    for idx in annotation_df.index:
        rna_type = []
        if ((idx in uaRNA_dict.keys()) and not (idx in srtRNA_dict.keys())):
            uaRNA_annot.append('|'.join(uaRNA_dict[idx]))
            rna_type.append('uaRNA')
        else:
            uaRNA_annot.append('.')
        if ((idx in coRNA_dict.keys()) and not (idx in srtRNA_dict.keys())):
            coRNA_annot.append('|'.join(coRNA_dict[idx]))
            rna_type.append('coRNA')
        else:
            coRNA_annot.append('.')

        if idx in srtRNA_dict.keys():
            rtRNA_annot.append('|'.join(srtRNA_dict[idx]))
            rna_type.append('srtRNA')
        else:
            if (idx in rtRNA_dict.keys()) and (idx not in uaRNA_dict.keys()) and (idx not in coRNA_dict.keys()):
                rtRNA_annot.append('|'.join(rtRNA_dict[idx]))
                rna_type.append('rtRNA')
            else:
                rtRNA_annot.append('.')
            #rtRNA_annot.append('.')
        if idx in lncRNA_dict.keys():
            lncRNA_annot.append('|'.join(lncRNA_dict[idx]))
            rna_type.append('lncRNA')
        else:
            lncRNA_annot.append('.')
            
        if idx in snoRNA_dict.keys():
            snoRNA_annot.append('|'.join(snoRNA_dict[idx]))
            rna_type.append('snoRNA')
        else:
            snoRNA_annot.append('.')
            
        if idx in pseudo_dict.keys():
            pseudo_annot.append('|'.join(pseudo_dict[idx]))
            rna_type.append('pseudo')
        else:
            pseudo_annot.append('.')
        if (idx in  pd.Index(incRNA.gene_name)) and ((idx not in coRNA_dict.keys()) and (idx not in ctRNA_dict.keys()) and (idx not in uaRNA_dict.keys()) and ((idx not in rtRNA_dict.keys())) and (idx not in srtRNA_dict.keys()) ):
            rna_type.append('incRNA')
        if (idx not in srtRNA_dict.keys()) and (idx not in rtRNA_dict.keys()) and (idx not in lncRNA_dict.keys()) and (idx not in pseudo_dict.keys()) and (idx in ctRNA_dict.keys()): 
            ctRNA_annot.append('|'.join(ctRNA_dict[idx]))
            rna_type.append('ctRNA')
        else:
            ctRNA_annot.append('.')
        annot.append(','.join(rna_type))
        
    annotation_df['rna_type'] = annot
    annotation_df['uaRNA'] = uaRNA_annot
    annotation_df['coRNA'] = coRNA_annot
    annotation_df['rtRNA'] = rtRNA_annot
    annotation_df['ctRNA'] = ctRNA_annot
    annotation_df['lncRNA'] = lncRNA_annot
    annotation_df['snoRNA'] = snoRNA_annot
    annotation_df['pseudogene'] =pseudo_annot
    
    annotation_df.to_csv(dir_name + '/annotation/ncRNA.annotation.tab.gz', sep='\t', index=True, header=True)






    