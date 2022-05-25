import pandas as pd

ProCap = pd.read_csv('ProCapAnalysis/ProCap.CountTable.hg38.bed', sep='\t')

ProCap_enhancer = ProCap.loc[[x != 'promoter' for x in ProCap.pid]]

eRNA_SAF = pd.DataFrame()
eRNA_SAF['GeneID'] = list(ProCap_enhancer.pid)
eRNA_SAF['Chr'] = list(ProCap_enhancer['#Chr'])
eRNA_SAF['Start'] = list(ProCap_enhancer.start - 200)
eRNA_SAF['End'] = list(ProCap_enhancer.end + 200)
eRNA_SAF['Strand'] = list(ProCap_enhancer.strand)

eRNA_SAF.to_csv('../data/eRNA.saf', index=False, header=True, sep='\t')

eRNA_SAF = pd.DataFrame()

promoter_name = [x + '_+' for x in ProCap_enhancer.pid]
promoter_name += [x + '_-' for x in ProCap_enhancer.pid]

eRNA_SAF['GeneID'] = promoter_name
eRNA_SAF['Chr'] = list(ProCap_enhancer['#Chr']) + list(ProCap_enhancer['#Chr'])
eRNA_SAF['Start'] = list(ProCap_enhancer.start - 200) + list(ProCap_enhancer.start - 200)
eRNA_SAF['End'] = list(ProCap_enhancer.end + 200) + list(ProCap_enhancer.end + 200)
eRNA_SAF['Strand'] = (['+']*len(ProCap_enhancer.pid)) + (['-']*len(ProCap_enhancer.pid))

eRNA_SAF.to_csv('../data/eRNA_both_strands.saf', index=False, header=True, sep='\t')
