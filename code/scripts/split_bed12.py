import numpy as np
import pandas as pd

bed12_names = ['chrom', 'start', 'end', 'gene', 'score', 'strand',
               'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
               'blockSizes', 'blockStarts'
              ]
X = pd.read_csv('Metaplots/AssayProfiles/References/ExpressedGeneList.bed12', sep='\t', names = bed12_names)

RPKM = pd.read_csv('QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstRepsUnstandardized.qqnorm.bed.gz', 
                   sep='\t', index_col=3)
samples = [x for x in RPKM.columns[5:] if (('.bam' not in x) and (x != 'NA18855'))]

n = 3500

CDS_length = []
for idx, row in X.iterrows():
    blockSizes = np.sum([int(x) for x in row.blockSizes.rstrip(',').split(',')])
    CDS_length.append(blockSizes)
    
X['CDS_length'] = CDS_length
X['gene_length'] = X.end - X.start

for i in range(4):
    
    Qname = 'Q' + str(i+1)
    
    eQi = RPKM[samples].median(axis=1).sort_values()[i*n:((i+1)*n)].index
    X.loc[X.gene.isin(eQi), bed12_names].to_csv(
        'Metaplots/AssayProfiles/References/ExpressedGeneList.expression.{Qname}.bed12'.format(Qname=Qname),
        index = False, header=False, sep='\t'
    )
    
    X.sort_values('CDS_length')[i*n:((i+1)*n)][bed12_names].to_csv(
        'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.{Qname}.bed12'.format(Qname=Qname),
        index = False, header=False, sep='\t'
    )
    
    X.sort_values('gene_length')[i*n:((i+1)*n)][bed12_names].to_csv(
        'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.{Qname}.bed12'.format(Qname=Qname),
        index = False, header=False, sep='\t'
    )

# Q1 = RPKM[samples].median(axis=1).sort_values()[:3500].index
# Q2 = RPKM[samples].median(axis=1).sort_values()[3500:7000].index
# Q3 = RPKM[samples].median(axis=1).sort_values()[7000:10500].index
# Q4 = RPKM[samples].median(axis=1).sort_values()[10500:]

# X.loc[X.gene.isin(Q1)].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.expression.Q1.bed12',
#     index = False, header=False, sep='\t'
# )




# X.sort_values('CDS_length')[:n][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q1.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('CDS_length')[n:n*2][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q2.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('CDS_length')[n*2:n*3][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q3.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('CDS_length')[n*3:][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q4.bed12',
#     index = False, header=False, sep='\t'
# )

# X.sort_values('gene_length')[:n][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q1.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('gene_length')[n:n*2][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q2.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('gene_length')[n*2:n*3][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q3.bed12',
#     index = False, header=False, sep='\t'
# )
# X.sort_values('gene_length')[n*3:][bed12_names].to_csv(
#     'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q4.bed12',
#     index = False, header=False, sep='\t'
# )