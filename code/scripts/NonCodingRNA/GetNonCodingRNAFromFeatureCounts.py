import pandas as pd
import argparse 
from gtfparse import read_gtf

parser = argparse.ArgumentParser()
parser.add_argument('--phenotype', type=str, required=True)
# parser.add_argument('--ncRNA', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    phenotype = args.phenotype
    
    gtf = read_gtf('ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf')
    
    ncRNA_genes = gtf.loc[(gtf.gene_type.isin(['snoRNA', 'snRNA',
                                               'lncRNA', 'unprocessed_pseudogene', 
                                               'transcribed_unprocessed_pseudogene',
                                               'pseudogene', 'rRNA_pseudogene',
                                               'transcribed_processed_pseudogene',
                                               'transcribed_unitary_pseudogene',
                                               'transcribed_unprocessed_pseudogene',
                                               'translated_processed_pseudogene',
                                               'translated_unprocessed_pseudogene', 
                                               'unprocessed_pseudogene'
                                              ])) & (gtf.feature == 'gene')].gene_id
    
    counts = pd.read_csv('featureCounts/{phenotype}/Counts.txt'.format(phenotype=phenotype), 
                                 sep='\t', skiprows=1, index_col=0)
    
    ncRNA_counts = counts.loc[ncRNA_genes]
    
    ncRNA_counts.to_csv('featureCounts/{phenotype}_annotated_ncRNA/Counts.txt'.format(phenotype=phenotype),
                        sep='\t', index=True, header=True)
