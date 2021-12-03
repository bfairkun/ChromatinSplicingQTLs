import pandas as pd
import argparse 


parser = argparse.ArgumentParser()

parser.add_argument('--phenotype', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    phenotype = args.phenotype
        
    protein_coding = pd.read_csv('featureCounts/{phenotype}/Counts.txt'.format(phenotype=phenotype), 
                                 sep='\t', skiprows=1)

    eRNA = pd.read_csv('featureCounts/{phenotype}_eRNA/Counts.txt'.format(phenotype=phenotype),
                       sep='\t', skiprows=1)

    cheRNA = pd.read_csv('featureCounts/{phenotype}_cheRNA/Counts.txt'.format(phenotype=phenotype), 
                         sep='\t', skiprows=1)
    
    RNA_counts = pd.concat([protein_coding, eRNA, cheRNA])
    RNA_counts = RNA_counts[[x for x in RNA_counts.columns if x[-5:-1] != 'bam.']]
    RNA_cols = [x for x in RNA_counts.columns if (('/2/' not in x) and ('/3/' not in x))]
    RNA_counts = RNA_counts[RNA_cols]
    
    RNA_counts.columns = [x if x[-3:] != 'bam' else x.split('/')[-3] for x in RNA_counts.columns]
    
    RNA_cols = list(RNA_counts.columns[:6])
    samples = RNA_counts.columns[6:]
    
    
    Fastq_samples = pd.read_csv('config/samples.tsv', sep='\t', comment='#')
    igsr_samples = pd.read_csv('../data/igsr_samples.tsv.gz', sep='\t', index_col=0)
    
    
    YRI_samples = list(igsr_samples.loc[samples].loc[igsr_samples.loc[samples]['Population code'] == 'YRI'].index)
    
    RNA_columns = RNA_cols + YRI_samples
    
    RNA_counts = RNA_counts[RNA_columns]
    
    
    RNA_counts.to_csv('QTLs/QTLTools/{phenotype}.AllRNA/Counts.txt'.format(phenotype=phenotype), sep='\t',
                      index=False, header=True)