import numpy as np
import pandas as pd
import os

import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--length', type=int, required=True)
parser.add_argument('--nstates', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    length = args.length
    nstates = args.nstates
    
    annotation_list = os.listdir("NonCodingRNA/tables/")
    annotation_list = [x for x in annotation_list if ("." + nstates + "states.ncRNA.bed.gz") in x]
    
    df = pd.DataFrame([])
    
    for annot in annotation_list:
        df_annot = pd.read_csv("NonCodingRNA/tables/" + annot, 
                               sep='\t', names = ['chrom', 'start', 'end', 'name', 'score', 'strand'])
        
        df = pd.concat([df, df_annot], axis=0)
        
    df.to_csv('NonCodingRNA/annotation/ncRNA.' + nstates + 'states.bed.gz', sep='\t', header=False, index=False)
    
    df.loc[(df.end - df.start) >= length].to_csv('NonCodingRNA/annotation/ncRNA_filtered.' + nstates + 'states.bed.gz', 
                                                sep='\t', header=False, index=False)
    
    
    
    
    annotation_list = os.listdir("NonCodingRNA/tables/")
    annotation_list = [x for x in annotation_list if (".merged_" + nstates + "states.tab.gz") in x]
    
    df = pd.DataFrame([])
    
    for annot in annotation_list:
        df_annot = pd.read_csv("NonCodingRNA/tables/" + annot, 
                               sep='\t', names = ['chrom', 'start', 'end', 'name', 'score', 'strand'])
        
        df = pd.concat([df, df_annot], axis=0)
        
    df.to_csv('NonCodingRNA/annotation/allHMM.' + nstates + 'states.bed.gz', sep='\t', header=False, index=False)
    
    