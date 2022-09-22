import pandas as pd
import argparse 

def get_histone_annotation(histone_file, index_col, annotation):
    histone_idx = pd.read_csv(histone_file, 
                      sep='\t',  names = ['chrom', 'start', 'end', 'name', 'strand'], index_col=index_col).index
    histone_annotation = [1 if x in histone_idx else 0 for x in annotation.index]
    
    return histone_annotation

parser = argparse.ArgumentParser()
parser.add_argument('--hmm', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    hmm = args.hmm
    
    annotation = pd.read_csv('NonCodingRNA{hmm}/annotation/ncRNA.annotation.tab.gz'.format(hmm=hmm), sep='\t',  index_col=0)
    k4me1 = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K4ME1.overlaps.bed.gz'.format(hmm=hmm), 
                                   3, annotation)
    k4me3 = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K4ME3.overlaps.bed.gz'.format(hmm=hmm), 
                                   3, annotation)
    k27ac = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K27AC.overlaps.bed.gz'.format(hmm=hmm), 
                                   3, 
                                           annotation)

    annotation['H3K4ME1'] = k4me1
    annotation['H3K4ME3'] = k4me3
    annotation['H3K27AC'] = k27ac
    
    k4me1 = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K4ME1.TSS_overlaps.bed.gz'.format(hmm=hmm), 
                                   3, annotation)
    k4me3 = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K4ME3.TSS_overlaps.bed.gz'.format(hmm=hmm), 
                                   3, annotation)
    k27ac = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/ncRNA.H3K27AC.TSS_overlaps.bed.gz'.format(hmm=hmm), 
                                   3,  annotation)

    annotation['H3K4ME1_TSS'] = k4me1
    annotation['H3K4ME3_TSS'] = k4me3
    annotation['H3K27AC_TSS'] = k27ac

    annotation.to_csv('NonCodingRNA{hmm}/annotation/ncRNA.histone.tab.gz'.format(hmm=hmm), sep='\t', index=True, header = True)
    
    allgenes = pd.read_csv('NonCodingRNA{hmm}/annotation/allGenes.TSS.bed.gz'.format(hmm=hmm), sep='\t',  index_col=4,
                      names = ['chrom', 'start', 'end', 'name', 'strand'])

    k4me1_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K4ME1.overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    k4me3_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K4ME3.overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    k27ac_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K27AC.overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    allgenes['H3K4ME1'] = k4me1_genes
    allgenes['H3K4ME3'] = k4me3_genes
    allgenes['H3K27AC'] = k27ac_genes
    
    k4me1_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K4ME1.TSS_overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    k4me3_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K4ME3.TSS_overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    k27ac_genes = get_histone_annotation('NonCodingRNA{hmm}/annotation/histone_marks/allGenes.H3K27AC.TSS_overlaps.bed.gz'.format(hmm=hmm), 4,
                                   allgenes)

    allgenes['H3K4ME1_TSS'] = k4me1_genes
    allgenes['H3K4ME3_TSS'] = k4me3_genes
    allgenes['H3K27AC_TSS'] = k27ac_genes

    allgenes_out = allgenes[['name', 'H3K4ME1', 'H3K4ME3', 'H3K27AC', 'H3K4ME1_TSS', 'H3K4ME3_TSS', 'H3K27AC_TSS']]
    
    allgenes_out.to_csv('NonCodingRNA{hmm}/annotation/allGenes.histone.tab.gz'.format(hmm=hmm), sep='\t', index=True, header = True)
    
    
    
    