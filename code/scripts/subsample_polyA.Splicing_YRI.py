import pandas as pd

qqnorm_splicing = pd.read_csv('QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.sorted.qqnorm.bed.gz', 
                             sep='\t')

samples = pd.read_csv('../data/igsr_samples.tsv.gz', sep='\t', index_col=0)

YRI_samples = qqnorm_splicing.columns.intersection(samples.loc[samples['Population code'] == 'YRI'].index)

YRI_columns = list(qqnorm_splicing.columns[:6]) + list(YRI_samples)
qqnorm_splicing[YRI_columns].to_csv('QTLs/QTLTools/polyA.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz', sep='\t',
                                index=False, header=True)
