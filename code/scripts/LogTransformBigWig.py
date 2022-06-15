import numpy as np
import pandas as pd
import argparse

def process_bedgraph(bedgraph, bedgraph_out):
    with open(bedgraph_out, 'w') as out:
        with open(bedgraph, 'r') as fh:
            for line in fh:
                if line[0] == '#':
                    print(line)
                    out.write(line)
                else:
                    bed_line = line.rstrip().split('\t')
                    bed_line[3] = str(round(np.log1p(float(bed_line[3])), 2))
                    bed_out = '\t'.join(bed_line) + '\n'
                    out.write(bed_out)
                    
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    
    process_bedgraph(input_file, output_file)
