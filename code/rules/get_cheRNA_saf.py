import numpy as np
import pandas as pd
import subprocess as sp

sp.run('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83531/suppl/GSE83531_SupplementaryTable_1.xlsx',
      shell=True)

xls = pd.ExcelFile("GSE83531_SupplementaryTable_1.xlsx") 
sp.run('rm GSE83531_SupplementaryTable_1.xlsx', shell=True)

# first table contains identifier cheRNA annotations for K562 cells; hg38
cheRNA = xls.parse(0) 
cheRNA['geneID'] = [x.strip() for x in cheRNA.geneID]
cheRNA = cheRNA.loc[cheRNA['RNA category'] == 'cheRNA']

GeneID = []
Chr = []
Start = []
End = []

cheSAF = pd.DataFrame()

for idx, row in cheRNA.iterrows():
    
    GeneID.append(row.geneID.strip())
    chrom, loc = row.location.strip().split(':')
    start, end = loc.split('-')
    Chr.append(chrom)
    Start.append(start)
    End.append(end)
    
cheSAF['GeneID'] = GeneID
cheSAF['Chr'] = Chr
cheSAF['Start'] = Start
cheSAF['End'] = End
cheSAF['Strand'] = ['.']*len(GeneID)


cheSAF= cheSAF.groupby(['Chr', 'Start', 'End']).agg({'GeneID':'first', 
                             'Strand':'first' }).reset_index().sort_values(['Chr', 'Start', 'End'])

cheSAF = cheSAF[['GeneID', 'Chr', 'Start', 'End', 'Strand']]

cheSAF.to_csv('../data/cheRNA_K562_GSE83531.saf', sep='\t', index=False, header=True)
