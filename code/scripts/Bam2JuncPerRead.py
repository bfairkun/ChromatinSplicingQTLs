######################################################################
# @author      : cnajar (cnajar@midway2-login2.rcc.local)
# @file        : Bam2JuncPerRead
# @created     : Monday Jan 23, 2023 16:50:00 CST
#
# @description : Takes a BAM file  as input. Outputs a bed file with all 
# the junctions found for each read. Designed for long read sequencing, 
# but it should also work for short reads
######################################################################

import pandas as pd
import pysam
import argparse 

def getSpliceSites(entry, min_len = 70, max_len = 100000):
    '''
    Takes a read and looks for coordinates of splice junctions using
    the CIGAR tuple. 
    '''
    cigarList = entry.cigartuples
    X=entry.pos
    cummulative=0
    coordinates=[]
    for (code,seqlen) in cigarList:
        if (code == 0) or (code == 2):
            # match or deletions
            cummulative+=seqlen ### add lengh of match or deletion
        elif code == 3:
            # skip (junction)
            five_prime_ss = str(X+cummulative)
            cummulative+=seqlen  ### add the intron length
            three_prime_ss = str(X+cummulative+1) ### 3' exon start (prior exon splice-site + intron length)
            if seqlen > max_len: ### ignore very long introns
                return None
            elif seqlen > min_len: ### only report introns longer than min length
                coordinates.append([five_prime_ss,three_prime_ss])
    return coordinates

def getChromAndStrand(bamf, entry):
    '''
    Get chromosome and strand information for read
    '''
    chrom = bamf.getrname( entry.rname )
    if entry.is_reverse:
        strand = '-'
    else:
        strand = '+'
    return chrom, strand

def getEntrySpliceSitesBed(bamf, entry, read_name):
    '''
    Create BED file rows for a read. Each entry corresponds to 
    a splice junction.
    '''
    chrom, strand = getChromAndStrand(bamf, entry)
    coord = getSpliceSites(entry)
    if not coord:
        return None
    bed_rows = [[chrom, coord[i][0], coord[i][1], 
                 read_name + ':' + str(i+1), read_name, strand] for i in range(len(coord))]
    return bed_rows
    
def bed_rows_to_str(bed_rows):
    '''
    BED rows list to str
    '''
    bed_rows_str = '\n'.join(['\t'.join(x) for x in bed_rows]) + '\n'
    return bed_rows_str

def BamFileToJuncBed(bam_path, output):
    '''
    Open BAM file and write BED file with all junctions
    per read. Each read can have multiple junctions. Each
    junction can appear repeated in multiple reads.
    '''
    bamf = pysam.Samfile(bam_path, "rb" )
    autosomes = ['chr' + str(i) for i in range(1, 23)]
    with open(output, 'w') as fh:
        read_counts = 1
        for entry in bamf.fetch():
            chrom = bamf.getrname( entry.rname )
            if chrom not in autosomes:
                continue
            if 'N' in entry.cigarstring:
                read_name = 'LongRead.' + str(read_counts)
                bed_rows = getEntrySpliceSitesBed(bamf, entry, read_name)
                if not bed_rows:
                    continue
                bed_str = bed_rows_to_str(bed_rows)
                fh.write(bed_str)
                read_counts += 1
    
def main(args):
    bam_path = args.input
    output = args.output
    BamFileToJuncBed(bam_path, output)
    
    
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)

if __name__ == '__main__':
    
    args = parser.parse_args()
    main(args)
    