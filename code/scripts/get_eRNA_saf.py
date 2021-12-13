import pandas as pd

ProCap = pd.read_csv('ProCapAnalysis/CountTable.hg38.features.bed', sep='\t',
                   names = ['chrom', 'start', 'end', 'names', 'annotation', 'strand'])

ProCap_enhancer = ProCap.loc[ProCap.annotation != 'promoter']

eRNA_SAF = pd.DataFrame()
eRNA_SAF['GeneID'] = list([str(list(x)[1].annotation) + '_' + str(list(x)[1].names) for x in ProCap_enhancer.iterrows()])
eRNA_SAF['Chr'] = list(ProCap_enhancer.chrom)
eRNA_SAF['Start'] = list(ProCap_enhancer.start - 200)
eRNA_SAF['End'] = list(ProCap_enhancer.end + 200)
eRNA_SAF['Strand'] = list(ProCap_enhancer.strand)

eRNA_SAF.to_csv('../data/eRNA.saf', index=False, header=True, sep='\t')
