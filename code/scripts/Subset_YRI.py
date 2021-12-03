import pandas as pd
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--phenotype', type=str, required=True)
parser.add_argument('--ncRNA', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    phenotype = args.phenotype
    ncRNA = args.ncRNA
        
    ncRNA_counts = pd.read_csv('featureCounts/{phenotype}_{ncRNA}/Counts.txt'.format(phenotype=phenotype, ncRNA=ncRNA), 
                                 sep='\t', skiprows=1, index_col=0)
    
    ncRNA_counts = ncRNA_counts[[x for x in ncRNA_counts.columns if x[-5:-1] != 'bam.']]
    RNA_cols = [x for x in ncRNA_counts.columns if (('/2/' not in x) and ('/3/' not in x))]
    ncRNA_counts = ncRNA_counts[RNA_cols]
    
    ncRNA_counts.columns = [x if x[-3:] != 'bam' else x.split('/')[-3] for x in ncRNA_counts.columns]
    
    RNA_cols = list(ncRNA_counts.columns[:6])
    samples = ncRNA_counts.columns[6:]
    
    
    Fastq_samples = pd.read_csv('config/samples.tsv', sep='\t', comment='#')
    igsr_samples = pd.read_csv('../data/igsr_samples.tsv.gz', sep='\t', index_col=0)
    
    
    samples = pd.Index(samples).intersection(igsr_samples.index)
    
    print(len(samples))
    
    igsr_samples = igsr_samples.loc[samples]
    
    YRI_samples = list(igsr_samples.loc[igsr_samples['Population code'] == 'YRI'].index)
    
    RNA_columns = RNA_cols + YRI_samples
    
    ncRNA_counts = ncRNA_counts[RNA_columns]
    
    
    ncRNA_counts.to_csv('featureCounts/{phenotype}_{ncRNA}.Subset_YRI/Counts.txt'.format(phenotype=phenotype, ncRNA=ncRNA), sep='\t',
                      index=True, header=True)