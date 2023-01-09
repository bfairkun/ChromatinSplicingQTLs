import numpy as np
import pandas as pd


def get_bed12_slice(bed12, gene_dir, gene_name):
    bed12_slice = bed12.loc[bed12.gene.isin(
        gene_dir.loc[gene_dir.ensembl == gene_name].gene_symbol)
             ]
    bed12_slice = bed12_slice[bed12_slice.columns[:12]]
    return bed12_slice


def make_QTL_bed(ua_list, QTL_name, bed12, gene_dir):
    gene_list = []

    chrom_list = []
    start_list = []
    end_list = []
    var_id_list = []
    signed_effect_list = []
    strand_list = []
    direction_list = []

    fh_apQTLs_df = pd.DataFrame()
    fh_diQTLs_df = pd.DataFrame()

    fh_apQTLs_gene = open('NonCodingRNA/QTLs/apQTLs.{QTL_name}.uaGene.bed'.format(QTL_name=QTL_name), 'w')
    fh_diQTLs_gene = open('NonCodingRNA/QTLs/diQTLs.{QTL_name}.uaGene.bed'.format(QTL_name=QTL_name), 'w')

    fh_apQTLs_ncRNA = open('NonCodingRNA/QTLs/apQTLs.{QTL_name}.uaRNA.bed'.format(QTL_name=QTL_name), 'w')
    fh_diQTLs_ncRNA = open('NonCodingRNA/QTLs/diQTLs.{QTL_name}.uaRNA.bed'.format(QTL_name=QTL_name), 'w')

    for ua in ua_list:

        uaGene = uaRNA_df.loc[ua, 'uaGene']

        chrom = uaRNA_df.loc[ua, 'chrom']
        strand = uaRNA_df.loc[ua, 'uaGene_strand']

        nc_end = uaRNA_df.loc[ua, 'ncRNA_end']
        pc_end = uaRNA_df.loc[ua, 'uaGene_end']
        
        nc_start = uaRNA_df.loc[ua, 'ncRNA_start']
        pc_start = uaRNA_df.loc[ua, 'uaGene_start']

        if QTL_name == 'ncQTL':
            signed_effect = uaRNA_df.loc[ua, 'uaGene_ncQTL_slope']
        else:
            signed_effect = uaRNA_df.loc[ua, 'eQTL_slope']

        QTL_var = 'chr' + uaRNA_df.loc[ua, QTL_name + '_var']
        QTL_sign = uaRNA_df.loc[ua, QTL_name + '_sign']
        pc_slice = get_bed12_slice(bed12, gene_dir, uaGene)

        uaGene_row_list = [chrom, pc_start, pc_end, QTL_var, signed_effect, strand]
        uaRNA_row_list = [chrom, nc_start, nc_end, QTL_var, signed_effect, strand]
        uaGene_row = '\t'.join([str(x) for x in uaGene_row_list]) + '\n'
        uaRNA_row = '\t'.join([str(x) for x in uaRNA_row_list]) + '\n'

        if QTL_sign > 0:
            fh_apQTLs_df = pd.concat([fh_apQTLs_df, pc_slice], axis=0)
            fh_apQTLs_gene.write(uaGene_row)
            fh_apQTLs_ncRNA.write(uaRNA_row)
        else:
            fh_diQTLs_df = pd.concat([fh_diQTLs_df, pc_slice], axis=0)
            fh_diQTLs_gene.write(uaGene_row)
            fh_diQTLs_ncRNA.write(uaRNA_row)

    fh_apQTLs_df.to_csv('NonCodingRNA/QTLs/apQTLs.{QTL_name}.uaGene.bed12'.format(QTL_name=QTL_name), 
                        sep='\t', index=False, header=False)
    fh_diQTLs_df.to_csv('NonCodingRNA/QTLs/diQTLs.{QTL_name}.uaGene.bed12'.format(QTL_name=QTL_name), 
                        sep='\t', index=False, header=False)

    fh_apQTLs_gene.close()
    fh_diQTLs_gene.close()
    fh_apQTLs_ncRNA.close()
    fh_diQTLs_ncRNA.close()

    
if __name__ == '__main__':

    uaRNA_df = pd.read_csv('NonCodingRNA/QTLs/summary.uaRNA.tab.gz', sep='\t', index_col=0)

    bed12 = pd.read_csv('scripts/GenometracksByGenotype/PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz', 
                        sep='\t',
                        names = ['chrom', 'start', 'end', 'gene', 'score', 'strand', 'start_thick', 
                                'end_thick', 'rgb', 'chunks', 'chunk_len', 'chunk_start', 'gene_name']
                        )

    gene_dir = pd.read_csv('NonCodingRNA/annotation/genes.txt.gz', sep='\t',
                           names = ['ensembl', 'gene_symbol'])




    ua_eQTLs = uaRNA_df.index[(uaRNA_df.eQTL_q <= 0.1) & (uaRNA_df.uaRNA_eQTL_nom_pval <= 0.05)]
    ua_ncQTLs = uaRNA_df.index[(uaRNA_df.ncQTL_q <= 0.1) & (uaRNA_df.uaGene_ncQTL_nom_pval <= 0.05)]
    make_QTL_bed(ua_eQTLs, 'eQTL', bed12, gene_dir)
    make_QTL_bed(ua_ncQTLs, 'ncQTL', bed12, gene_dir)

