import pandas as pd

ProCap = pd.read_csv('QTLs/QTLTools/ProCap/OnlyFirstReps.qqnorm.bed.gz', sep='\t')
annotation = [x.split('_chr')[0] for x in ProCap.pid]
ProCap['annotation'] = annotation
ProCap_enhancer = ProCap.loc[ProCap.annotation != 'promoter']

eRNA_SAF = pd.DataFrame()
eRNA_SAF['GeneID'] = list(ProCap_enhancer.loc[ProCap_enhancer.annotation != 'promoter'].pid)
eRNA_SAF['Chr'] = list(ProCap_enhancer['#Chr'])
eRNA_SAF['Start'] = list(ProCap_enhancer.start - 200)
eRNA_SAF['End'] = list(ProCap_enhancer.end + 200)
eRNA_SAF['Strand'] = list(ProCap_enhancer.strand)

eRNA_SAF.to_csv('../data/eRNA.saf', index=False, header=True, sep='\t')
