import numpy as np
import pandas as pd
import os

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--distance', default = '-1', type=str, required=False)
parser.add_argument('--nstates', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    distance = args.distance
    nstates = args.nstates
    
#     annotation_list = os.listdir("NonCodingRNA/tables/")
#     annotation_list = [x for x in annotation_list if ("." + nstates + "states.ncRNA.bed.gz") in x]
    
#     df = pd.DataFrame([])
    
#     for annot in annotation_list:
#         df_annot = pd.read_csv("NonCodingRNA/tables/" + annot, 
#                                sep='\t', names = ['chrom', 'start', 'end', 'name', 'score', 'strand'])
        
#         df = pd.concat([df, df_annot], axis=0)
        
#     df.to_csv('NonCodingRNA/annotation/ncRNA.' + nstates + 'states.bed.gz', sep='\t', header=False, index=False)
    
#     df.loc[(df.end - df.start) >= distance].to_csv('NonCodingRNA/annotation/ncRNA_filtered.' + nstates + 'states.bed.gz', 
#                                                 sep='\t', header=False, index=False)
    
    
    
    
    annotation_list = os.listdir("NonCodingRNA/tables/")
    if distance == '-1':
        annotation_list = [x for x in annotation_list if (".merged_" + nstates + "states.bed.gz") in x]
    else:
        annotation_list = [x for x in annotation_list if (".merged_" + nstates + "states." + distance + "bp_merge.bed.gz") in x]
    
    df = pd.DataFrame([])
    
    for annot in annotation_list:
        df_annot = pd.read_csv("NonCodingRNA/tables/" + annot, 
                               sep='\t', names = ['chrom', 'start', 'end', 'name', 'score', 'strand'])
        
        df = pd.concat([df, df_annot], axis=0)
      
    if distance == '-1':
        df.to_csv('NonCodingRNA/annotation/allHMM.' + nstates + 'states.bed.gz', sep='\t', header=False, index=False)
    else:
        df.to_csv('NonCodingRNA/test/allHMM.' + nstates + 'states.' + distance + 'bp_merge.bed.gz', sep='\t', header=False, index=False)
    
    