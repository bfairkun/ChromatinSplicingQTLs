import pandas as pd
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--saf',  type=str, required=True)
parser.add_argument('--output',  type=str, required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    saf_file = args.saf
    output = args.output
    
    saf = pd.read_csv(saf_file, sep='\t')
    
    saf_neg = saf.copy()
    saf_neg.Strand = ['-' if x == '+' else '+' for x in saf.Strand]
    
    saf_neg.GeneID = [x + '_reverse' for x in saf.GeneID]
    
    saf_out = pd.concat([saf, saf_neg], axis=0, ignore_index=True)

    saf_out.to_csv(output, sep='\t', index=False, header=True)