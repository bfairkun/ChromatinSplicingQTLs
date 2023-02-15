import pandas as pd
import sys


if __name__ == '__main__':
    
    _, annot_junc_file, clusters_file = sys.argv
    
    ### Add annotations to junction annotation data
    
    annot = pd.read_csv(annot_junc_file, 
                        sep='\t', names = ['chrom', 'start', 'end', 'junction_id', 'ensembl_id', 'strand', 
                                           'ben_annot', 'yang_annot', 'gene_id', 'seq_5p', 'pwm_5p',
                                           'seq_3p', 'pwm_3p', 'seq_5p_x', 'maxentscan_5p', 'seq_3p_ext',
                                           'maxentscan_3p'
                                          ])

    loc_5p = []
    loc_3p = []
    for idx, row in annot.iterrows():
        if row.strand == '+':
            five = row.chrom + ':' + str(row.start) + ':' + row.strand
            three = row.chrom + ':' + str(row.end) + ':' + row.strand
        else:
            five = row.chrom + ':' + str(row.end) + ':' + row.strand
            three = row.chrom + ':' + str(row.start) + ':' + row.strand
        loc_5p.append(five)
        loc_3p.append(three)

    annot ['id_5p'] = loc_5p
    annot ['id_3p'] = loc_3p

    annot['intron_length'] = annot.end - annot.start
    annot['count_exons'] = [1]*len(annot)
    
    annot.to_csv("SplicingAnalysis/NMDJunctions/Annotation/annotation_junctions.bed.gz", sep='\t',
                 index=False, header=True)
    
    # add junction id to clusters
    
    clusters = pd.read_csv(clusters_file, sep='\t',
                      names = ['chrom', 'start', 'end', 'strand', 'ben_annot', 'yang_annot',
                              'ensembl_id', 'gene_id', 'cluster'])

    clusters['junction_id'] = clusters[['chrom', 'start', 'end', 'strand']].astype(str).agg(':'.join, axis=1)
    
    clusters.to_csv("SplicingAnalysis/NMDJunctions/Annotation/annotation_clusters.bed.gz", sep='\t',
                 index=False, header=True)
    
        
        
        
        