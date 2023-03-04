import numpy as np
import pandas as pd

bed12_names = ['chrom', 'start', 'end', 'gene', 'score', 'strand',
               'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
               'blockSizes', 'blockStarts'
              ]
X = pd.read_csv('Metaplots/AssayProfiles/References/ExpressedGeneList.bed12', sep='\t', names = bed12_names)

CDS_length = []
for idx, row in X.iterrows():
    blockSizes = np.sum([int(x) for x in row.blockSizes.rstrip(',').split(',')])
    CDS_length.append(blockSizes)
    
X['CDS_length'] = CDS_length
X['gene_length'] = X.end - X.start

n = 3500
X.sort_values('CDS_length')[:n][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q1.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('CDS_length')[n:n*2][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q2.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('CDS_length')[n*2:n*3][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q3.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('CDS_length')[n*3:][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q4.bed12',
    index = False, header=False, sep='\t'
)

X.sort_values('gene_length')[:n][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q1.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('gene_length')[n:n*2][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q2.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('gene_length')[n*2:n*3][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q3.bed12',
    index = False, header=False, sep='\t'
)
X.sort_values('gene_length')[n*3:][bed12_names].to_csv(
    'Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q4.bed12',
    index = False, header=False, sep='\t'
)