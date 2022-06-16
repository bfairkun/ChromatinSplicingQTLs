import pysam
import argparse
import os

def read_has_long_junc(read, max_junc_length):
    cigar = read.cigartuples
    N = [x[1] for x in cigar if (x[0] == 3)]
    long_junc = any([x > max_junc_length for x in N])
    return long_junc


def process_bam(bam_input, bam_output, max_junc_length):
    print('processing bam file in:')
    print(bam_input)
    samfile = pysam.AlignmentFile(bam_input, "rb")
    sam_filtered = pysam.AlignmentFile(bam_output, "wb", template=samfile)
    counter = 0
    counter_long = 0
    for read in samfile.fetch():

        counter += 1
        if (counter % 1000000) == 0:
            print('processed ' + str(counter) + ' reads')
        long_junc = read_has_long_junc(read, max_junc_length)

        if long_junc:
            counter_long += 1
            continue

        else:
            sam_filtered.write(read)

    sam_filtered.close()
    samfile.close()
    
    print('Removed long junctions successfully')
    print('removed ' + str(counter_long) + ' reads')
    print('Stored results in:')
    print(bam_output)
    print('')
    

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--max_junc_length', type=int, required=True)

if __name__ == '__main__':
    
    if not os.path.isdir('NonCodingRNA'):
        os.mkdir('NonCodingRNA/')
        
    if not os.path.isdir('NonCodingRNA/bam/'):
        os.mkdir('NonCodingRNA/bam/')
    
    args = parser.parse_args()
    bam_input = args.input
    bam_output = args.output
    max_junc_length = args.max_junc_length
    
    process_bam(bam_input, bam_output, max_junc_length)
    