import numpy as np
import pandas as pd
# from tqdm import tqdm
import sys

def process_long_reads(longreads_file, df, split=False):
    gm = pd.read_csv(longreads_file, sep='\t', 
                  names = ['chrom', 'start', 'end', 'junction_id', 'read', 'strand'])
    
    gm.start = gm.start.astype(int)
    gm.end = gm.end.astype(int)
    
    gm.end -= 1
    
    exon_psi_concat = []
    
    chroms = ['chr' + str(i) for i in range(1, 23)]
    
    for chrom in chroms:
        
        print(chrom)
        
        exon_psi = merge_annot_longreads(df, gm, chrom)
        exon_psi_concat.append(exon_psi)
        
            
    exon_psi = pd.concat(exon_psi_concat)

    return exon_psi


def merge_annot_longreads(df, gm, chrom):
    
    df_slice = df.loc[df.chrom == chrom]
        
    gm_slice = gm.loc[gm.chrom == chrom]

    included = df_slice.merge(gm_slice, left_on=['chrom', 'donor_c1', 'receptor_c1'], right_on=['chrom', 'start', 'end'])
    
#     included = included.merge(gm_slice, left_on=['chrom', 'donor_c2', 'receptor_c2'], right_on=['chrom', 'start', 'end'])
    included = included.merge(gm_slice, left_on=['chrom', 'donor_c2', 'receptor_c2', 'read'], right_on=['chrom', 'start', 'end', 'read'])
    
    included = included.groupby('event').size().reset_index()
    included.columns = ['exon', 'included_juncs']

    skipped = df_slice.merge(gm_slice, left_on=['chrom', 'donor_c1', 'receptor_c2'], right_on=['chrom', 'start', 'end'])
    skipped = skipped.groupby('event').size().reset_index()
    skipped.columns = ['exon', 'skipped_juncs']

    

    exon_psi = pd.merge(skipped, included, how='outer', left_on='exon', right_on='exon').fillna(0)
    exon_psi['PSI'] = exon_psi.included_juncs/(exon_psi.included_juncs + exon_psi.skipped_juncs)
    exon_psi = exon_psi.set_index('exon')

    df__ = df_slice.set_index('event')
    exon_psi['LE_n'] = df__.loc[exon_psi.index].LE_n
    
    return exon_psi

if __name__ == "__main__":
    args = sys.argv[1:3]
    input_file = args[0]
    output_file = args[1]
    
    exons = pd.read_csv('../../EVENT_INFO-hg38.tab.gz', sep='\t')[['EVENT', 'COMPLEX', 'CO_C1', 'CO_A', 'CO_C2', 'LE_n']].dropna()
    exons = exons.loc[~exons.COMPLEX.isin(['MIC_S', 'IR', 'Alt3', 'Alt5', 'MIC-M'])]
    exons = exons.set_index('EVENT')

    impact = pd.read_csv('../../PROT_IMPACT-hg38-v3.tab.gz', sep='\t')
    impact = impact.set_index('EventID')
    impact = impact.loc[~impact.ONTO.isin(['NonCoding', "5' UTR", "3' UTR"])]


    exon = exons.index
    LE = exons.LE_n
    c1 = exons.CO_C1.str.split('-', expand=True)
    chrom = c1[0].str.split(':', expand=True)[0]
    donor_c1 = c1[1].astype(int)

    ca = exons.CO_A.str.split('-', expand=True)
    receptor_c1 = ca[0].str.split(':', expand=True)[1].astype(int) - 1
    donor_c2 = ca[1].astype(int)

    c2 = exons.CO_C2.str.split('-', expand=True)
    receptor_c2 = c2[0].str.split(':', expand=True)[1].astype(int) - 1

    exon = exons.index
    LE = exons.LE_n
    c1 = exons.CO_C1.str.split('-', expand=True)
    chrom = c1[0].str.split(':', expand=True)[0]
    donor_c1 = c1[1].astype(int)

    ca = exons.CO_A.str.split('-', expand=True)
    receptor_c1 = ca[0].str.split(':', expand=True)[1].astype(int) - 1
    donor_c2 = ca[1].astype(int)

    c2 = exons.CO_C2.str.split('-', expand=True)
    receptor_c2 = c2[0].str.split(':', expand=True)[1].astype(int) - 1


    df = pd.DataFrame()
    df['event'] = list(exon)
    df['chrom'] = list(chrom)
    df['donor_c1'] = list(donor_c1)
    df['receptor_c1'] = list(receptor_c1)
    df['donor_c2'] = list(donor_c2)
    df['receptor_c2'] = list(receptor_c2)
    df['LE_n'] = list(LE)
    df = df.loc[[x[:5]=='HsaEX' for x in df.event]]
    df = df.loc[df.event.isin(impact.loc[~impact.ONTO.isin(['NonCoding', "5' UTR", "3' UTR"])].index)]

    exon_psi = process_long_reads(input_file, df, split=True)
    exon_psi.to_csv(output_file, sep='\t', index=True, header=True)
