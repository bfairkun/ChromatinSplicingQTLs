import sys
import gzip

def process_vcf_line(line):
    if line.startswith('##'):
        return None
    
    elif line.startswith('#CHROM'):
        line_split = line.split('\t')
        line_header = '\t'.join(line_split[:4] + line_split[8:])
        return line_header
    else:
        line = line[3:].replace('0|0', '0').replace('1|1', '2').replace('0|1', '1').replace('1|0', '0')
        line_split = line.split('\t')
        line_out = '\t'.join(line_split[:4] + line_split[8:])
        return line_out
        
    
if __name__ == '__main__':
    _, f_in, f_out = sys.argv
    
    with gzip.open(f_out, 'wt') as fh:
        with gzip.open(f_in, 'rt') as vcf:
            for line in vcf:
                line_ = process_vcf_line(line)
                if line_:
                    fh.write(line_)
