#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TranslateLongReads
# @created     : Friday Mar 01, 2024 11:18:03 CST
#
# @description : 
######################################################################

import sys
import gzip
# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Arg1" ,"Arg2"]

_, MyArg1, MyArg2 = sys.argv

from Bio import SeqIO
from collections import deque

def flatten_ranges(range_list):
    flattened_list = []
    for r in range_list:
        flattened_list.extend(range(r.start, r.stop))
    return flattened_list

def read_gtf(gtf_file):
    start_codons_dict = {}
    opener = gzip.open if gtf_file.endswith('.gz') else open
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Skip invalid lines
            chrom = parts[0]
            feature_type = parts[2]
            strand = parts[6]
            if feature_type != 'start_codon':
                continue  # Skip if not a start codon feature
            if strand == '+':
                start_pos = int(parts[3]) - 1
            elif strand == '-':
                start_pos = int(parts[3]) + 2
            if chrom not in start_codons_dict:
                start_codons_dict[chrom] = set()
            start_codons_dict[chrom].add(start_pos)
    return start_codons_dict


def read_bed12(bed_file):
    transcripts = []
    opener = gzip.open if bed_file.endswith('.gz') else open
    with opener(bed_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            score = int(parts[4])
            strand = parts[5]
            block_sizes = tuple(list(map(int, parts[10].split(','))))
            block_starts = tuple(list(map(int, parts[11].split(','))))
            exons = tuple([(start + block_start, start + block_start + block_sizes[i]) for i, block_start in enumerate(block_starts)])
            transcripts.append({'chrom': chrom, 'start': start, 'end': end, 'name': name, 'score': score, 'strand': strand, 'exons': exons, 'block_sizes' : block_sizes, 'block_starts' : block_starts})
    return transcripts
transcripts = read_bed12("LongReads/FullLengthReadsBed12/GM1.bed.gz")
transcripts[0]

# import pdb; pdb.set_trace()

def get_relative_pos(transcript, absolute_pos):
    relative_pos = None
    exons = transcript['exons']
    if transcript['start'] <= absolute_pos <= transcript['end']:
        relative_pos = 0
        if transcript['strand'] == '+':
            for exon_start, exon_end in exons:
                if exon_start <= absolute_pos <= exon_end:
                    relative_pos += absolute_pos - exon_start
                    break
                elif absolute_pos > exon_end:
                    relative_pos += exon_end - exon_start
                elif absolute_pos < exon_start:
                    relative_pos = None
        elif transcript['strand'] == '-':
            for exon_start, exon_end in list(reversed(exons)):
                if exon_start <= absolute_pos <= exon_end:
                    relative_pos += exon_end - absolute_pos
                    break
                elif absolute_pos < exon_start:
                    relative_pos += exon_end - exon_start
                elif absolute_pos > exon_end:
                    relative_pos = None
    return relative_pos
print(transcripts[0])
print(get_relative_pos(transcripts[0], 959260))


def extract_exonic_sequence(transcript, seq_dict):
    chrom = transcript['chrom']
    strand = transcript['strand']
    exons = transcript['exons']
    exonic_sequence = ''
    intron_junctions = []
    if strand == '-':
        exons = list(reversed(exons))
    for i, (exon_start, exon_end) in enumerate(exons):
        exon_seq = seq_dict[chrom].seq[exon_start:exon_end]
        if strand == '-':
            exon_seq = exon_seq.reverse_complement()
        exonic_sequence += str(exon_seq)
        if i < len(transcript['exons']) - 1:
            intron_junctions.append((exon_end, transcript['exons'][i+1][0]))
            exonic_sequence += '|'
    return exonic_sequence, intron_junctions
print(transcripts[0])
extract_exonic_sequence(transcripts[0], seq_dict)

        # import pdb; pdb.set_trace()


def write_bed12_with_orf(transcripts, output_bed_file, start_codons_dict):
    with open(output_bed_file, 'w') as f:
        for transcript in transcripts:
            exonic_sequence, intron_junctions = extract_exonic_sequence(transcript, fasta_file)
            start_codon_pos = identify_start_codon(transcript, start_codons_dict)
            orf_sequence = 'M' if start_codon_pos is not None else ''
            for i in range(len(exonic_sequence)):
                if exonic_sequence[i] == '|':
                    orf_sequence += '|'
                elif i == start_codon_pos:
                    orf_sequence += 'M'
                elif i == len(exonic_sequence) - 1:
                    orf_sequence += '*'
            thick_start = transcript['start'] + start_codon_pos if start_codon_pos is not None else transcript['start']
            thick_end = transcript['end']
            f.write(f"{transcript['chrom']}\t{transcript['start']}\t{transcript['end']}\t{transcript['name']}\t{transcript['score']}\t{transcript['strand']}\t{thick_start}\t{thick_end}\t0\t{len(transcript['exons'])}\t{','.join(str(e[1] - e[0]) for e in transcript['exons'])}\t{','.join(str(e[0] - transcript['start']) for e in transcript['exons'])}\n")

def read_start_codons(start_codons_bed):
    start_codons_dict = {}
    with open(start_codons_bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            if chrom not in start_codons_dict:
                start_codons_dict[chrom] = set()
            start_codons_dict[chrom].add(pos)
    return start_codons_dict

# if __name__ == '__main__':
# bed_file = 'your_bed_file.bed'
# fasta_file = 'your_fasta_file.fasta'
# output_bed_file = 'output_bed_file.bed'
# start_codons_bed = 'start_codons.bed'
# transcripts = read_bed12(bed_file)
# start_codons_dict = read_start_codons(start_codons_bed)
seq_dict = SeqIO.to_dict(SeqIO.parse("ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa", 'fasta'))
transcripts = read_bed12("LongReads/FullLengthReadsBed12/GM1.bed.gz")
start_codons_dict = read_gtf("ReferenceGenome/Annotations/gencode.v37.chromasomal.annotation.gtf.gz")

def count_bars_until_n_position(sequence, n):
    bar_count = 0
    non_bar_count = 0
    for i,char in enumerate(sequence):
        if char == "|":
            bar_count += 1
        else:
            non_bar_count += 1
        if non_bar_count == n:
            break
    return bar_count

def identify_start_codon(transcript, start_codons_dict, seq_dict=None):
    start_codon_pos = None
    relative_start_codon_pos = None
    relative_start_codon_pos_with_offset = None
    start_codon_seq = None
    chrom = transcript['chrom']
    strand = transcript['strand']
    exons = transcript['exons']
    exonic_region = set(flatten_ranges([range(start, stop) for (start, stop) in exons]))
    if chrom in start_codons_dict: 
        valid_start_codons = sorted(list(exonic_region.intersection(start_codons_dict[chrom])))
        if valid_start_codons:
            if strand == '+':
                start_codon_pos = valid_start_codons[0]
            elif strand == '-':
                start_codon_pos = valid_start_codons[-1]
            relative_start_codon_pos = get_relative_pos(transcript, start_codon_pos)
            if seq_dict:
                sequence, _ = extract_exonic_sequence(transcript, seq_dict)
                relative_start_codon_pos_with_offset = count_bars_until_n_position(sequence, relative_start_codon_pos) + relative_start_codon_pos
                start_codon_seq = sequence.replace("|", "")[relative_start_codon_pos:relative_start_codon_pos+3]
    return start_codon_pos, relative_start_codon_pos, relative_start_codon_pos_with_offset, start_codon_seq

transcripts[0]['exons']
extract_exonic_sequence(transcripts[0], seq_dict)
identify_start_codon(transcripts[0], start_codons_dict, seq_dict)
identify_start_codon(transcripts[0], start_codons_dict, seq_dict)

### scratch code
for i, transcript in enumerate(transcripts[0:100]):
    start_codon_pos, relative_start_codon_pos, relative_start_codon_pos_with_offset, start_codon_seq = identify_start_codon(transcript, start_codons_dict, seq_dict)
    print(i, transcript['chrom'],transcript['start'], transcript['end'], transcript['strand'], start_codon_pos, relative_start_codon_pos, relative_start_codon_pos_with_offset, start_codon_seq)
    # print(extract_exonic_sequence(transcript, seq_dict))
    # if relative_start_codon_pos != relative_start_codon_pos_with_offset:
    #     print(i, transcript['strand'], start_codon_pos, relative_start_codon_pos, relative_start_codon_pos_with_offset, start_codon_seq) 
    

import re
def find_orf(dna_sequence):
    # Find ORF using regex
    orf_match = re.search(r"(?:[^|]*)(A|?T|?G(?:[^|]{3})*?)(?:TAA|TAG|TGA)", dna_sequence)
    if orf_match:
        return orf_match.group(1)
    else:
        return None

# Test the function
dna_sequence = "ATC|GATCGAT|GAT|GTTAAAT|GCTAGCTAG"
orf = find_orf(dna_sequence)
print("ORF:", orf)

