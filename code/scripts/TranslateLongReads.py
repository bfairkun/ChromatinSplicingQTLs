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
from Bio import SeqIO
import re

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa","LongReads/FullLengthReadsBed12/GM1.bed.gz","ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz","SplicingAnalysis/regtools_annotate_combined/comprehensive.bed.gz","scratch/MarkedORFs_first.bed" ,"scratch/MarkedORFs_longest.bed"]

_, f_fasta, f_bed12, f_gtf, f_intronsbed, file1name, file2name = sys.argv
print(sys.argv)

# file1name = "scratch/MarkedORFs_first.bed"
# file2name = "scratch/MarkedORFs_longest.bed"
# f_intronsbed = "SplicingAnalysis/regtools_annotate_combined/comprehensive.bed.gz"
# f_fasta = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
# f_bed12 = "LongReads/FullLengthReadsBed12/GM1.bed.gz"
# f_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz"


def flatten_ranges(range_list):
    flattened_list = []
    for r in range_list:
        flattened_list.extend(range(r.start, r.stop))
    return flattened_list

def read_bed_file(filename):
    chromosome_dict = {}
    if filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(filename, 'rt') as bed_file:
        for line in bed_file:
            if line.startswith('#'):  # Skip comment lines
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:  # Skip malformed lines
                continue
            chromosome = fields[0]
            start = fields[1]
            end = fields[2]
            if not start.isdigit() or not end.isdigit():  # Check if start and end are integers
                continue
            start = int(start)
            end = int(end)
            if chromosome not in chromosome_dict:
                chromosome_dict[chromosome] = set()
            chromosome_dict[chromosome].add((start, end))
    return chromosome_dict

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
# transcripts = read_bed12("LongReads/FullLengthReadsBed12/GM1.bed.gz")
# transcripts[0]

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
# print(transcripts[0])
# print(get_relative_pos(transcripts[0], 959260))

def get_absolute_pos(transcript, relative_pos, replace_None=True):
    if replace_None:
        absolute_pos=transcript['start']
    else:
        absolute_pos = None
    exons = transcript['exons']
    if 0 <= relative_pos < transcript['end'] - transcript['start']:
        cumulative_block_size = 0
        if transcript['strand'] == '+':
            for i, block_size in enumerate(transcript['block_sizes']):
                remainder = relative_pos - cumulative_block_size
                cumulative_block_size += block_size
                if cumulative_block_size >= relative_pos:
                    absolute_pos = transcript['start'] + remainder + transcript['block_starts'][i]
                    break
        elif transcript['strand'] == '-':
            for i, block_size in enumerate(list(reversed(transcript['block_sizes']))):
                remainder = relative_pos - cumulative_block_size
                cumulative_block_size += block_size
                if cumulative_block_size >= relative_pos:
                    absolute_pos = transcript['start'] + list(reversed(transcript['block_starts']))[i] + block_size - remainder
                    break
    return absolute_pos
# print(transcripts[0])
# print(get_absolute_pos(transcripts[0], 100))
# print(get_absolute_pos(transcripts[0], 0))
# print(get_absolute_pos(transcripts[10], 100))
# print(get_absolute_pos(transcripts[10], 1000))


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
            if strand == '+':
                intron_junctions.append((exon_end, exons[i+1][0]+1))
            else:
                intron_junctions.append((exons[i+1][1], exon_start+1))
            exonic_sequence += '|'
    return exonic_sequence, intron_junctions
# print(transcripts[0])
# extract_exonic_sequence(transcripts[0], seq_dict)

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

def identify_start_codon(transcript, start_codons_dict):
    start_codon_pos = None
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
    return start_codon_pos

def insert_marks_for_longset_ORF(sequence):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for longest ORF
    """
    try:
        longest_orf_match = max(re.findall(r"(?=(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)))", sequence, flags=re.IGNORECASE), key = len)
    except ValueError:
        longest_orf_match = None
    if longest_orf_match:
        start_codon_pos = sequence.find(longest_orf_match)
        stop_codon_pos = start_codon_pos + len(longest_orf_match)
        return sequence[0:start_codon_pos] + "^" + longest_orf_match + "*" + sequence[stop_codon_pos:]
    else:
        return sequence
# insert_marks_for_longset_ORF("GGGATGAAAGGGAAA|GGG|T|AA|GGGAAA")

def insert_marks_for_defined_ORF(sequence, start_codon_pos=None):
    """
    return sequence with "^" to mark start codon at provided string position, "*" to mark stop for longest ORF. If no stop codon is found, then no * will be inserted.
    """
    if start_codon_pos and start_codon_pos < len(sequence):
        # orf_match = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        # orf_match_nostop = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match_nostop = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        if orf_match:
            # orf match with start and stop codon
            return sequence[:start_codon_pos] + "^" + orf_match.group(0) + "*" + sequence[orf_match.span(0)[1] + start_codon_pos:]
        elif orf_match_nostop:
            # orf match with no stop
            return sequence[:start_codon_pos] + "^" + sequence[start_codon_pos:]
        else:
            return sequence
    else:
        return sequence
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|T|AA|GGGAAA", 3)
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|TT|AA|GGGAAA", 3)

def get_NMD_detective_B_classification(sequence):
    """
    sequence should be marked with '^' for start, '*' for stop, and '|' for splice juncs
    """
    CDS = re.search(r"\^(\w+)\*", sequence.replace("|", ""))
    InternalStopExon = re.search(r"(^|\|)([\^ACGTNacgtn]*\*[ACGTNacgtn]*)\|", sequence)
    if "^" not in sequence:
        return "No CDS"
    elif re.search(r"\^[\w|]+$", sequence):
        return "No stop"
    elif re.search(r"(^|\|)[\^ACGTNacgtn]*\*[ACGTNacgtn]*$", sequence):
        return "Last exon"
    elif len(CDS.group(1)) <= 125:
        return "Start proximal"
    elif len(InternalStopExon.group(2)) >= 407:
        return "Long exon"
    elif re.search(r"\*[ACGTNacgtn]{0,50}\|[ACGTNacgtn]+$", sequence):
        return "50 nt rule"
    else:
        return "Trigger NMD"

def get_thickStart_thickStop(transcript, sequence):
    """
    given transcript and sequence with ORF marked if present, output thickStart and thickEnd
    """
    CDS = re.search(r"\^([ATCGN|]+)\*?", sequence.replace("|", ""), flags=re.IGNORECASE)
    if CDS:
        CDS_relative_start, CDS_relative_stop = CDS.span(0)
    else:
        CDS_relative_start, CDS_relative_stop = 0, 0
    start_pos = get_absolute_pos(transcript, CDS_relative_start)
    stop_pos = get_absolute_pos(transcript, CDS_relative_stop)
    return tuple(sorted([start_pos, stop_pos]))
# get_thickStart_thickStop(transcripts[0], "GGG^ATGAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAA|GGG|T|AA*|GGGAA")
# get_thickStart_thickStop(transcripts[0], "GGG^ATGAAAGGGAAA|GGG|T|AA|GGGAAAA")
# get_thickStart_thickStop(transcripts[0], "GGGATGAAAGGGAAA|GGG|T|AA|GGGAA")

def get_CDS_marked_transcript_seqs(transcript, seq_dict, start_codons_dict):
    sequence, introns = extract_exonic_sequence(transcript, seq_dict)
    start_codon_pos = identify_start_codon(transcript, start_codons_dict)
    if start_codon_pos:
        relative_start_codon_pos = get_relative_pos(transcript, start_codon_pos)
        relative_start_codon_pos_with_offset = count_bars_until_n_position(sequence, relative_start_codon_pos) + relative_start_codon_pos
        start_codon_seq = sequence.replace("|", "")[relative_start_codon_pos:relative_start_codon_pos+3]
        sequence_for_first_ORF = insert_marks_for_defined_ORF(sequence, relative_start_codon_pos_with_offset)
    else:
        sequence_for_first_ORF = sequence
    sequence_for_longest_ORF = insert_marks_for_longset_ORF(sequence)
    return sequence_for_first_ORF, sequence_for_longest_ORF


seq_dict = SeqIO.to_dict(SeqIO.parse(f_fasta, 'fasta'))
transcripts = read_bed12(f_bed12)
introns_dict = read_bed_file(f_intronsbed)
start_codons_dict = read_gtf(f_gtf)

NMD_classification_to_rgb = {'No CDS':'150,150,150', 'No stop':'152,78,163', 'Last exon':'55,126,184', 'Start proximal':'255,127,0', 'Long exon':'255,255,51', '50 nt rule':'77,175,74', 'Trigger NMD':'228,26,28'}

with open(file1name, 'w') as f_annotatedStart, open(file2name, 'w') as f_longestORF:
    for i, transcript in enumerate(transcripts):
        if i >= 0:
            # if transcript['name'] == 'm54304_190613_163530/37946044/ccs':
            if True:
            # if i < 10000:
                # print(i)
                sequence, introns = extract_exonic_sequence(transcript, seq_dict)
                try:
                    identifiable_introns = [intron in introns_dict[transcript['chrom']] for intron in introns]
                except KeyError:
                    identifiable_introns = [None]
                sequence_for_first_ORF, sequence_for_longest_ORF = get_CDS_marked_transcript_seqs(transcript, seq_dict, start_codons_dict)
                firstStart_NMD = get_NMD_detective_B_classification(sequence_for_first_ORF)
                longestORF_NMD = get_NMD_detective_B_classification(sequence_for_longest_ORF)
                thickStart_firstStart, thickStop_firstStart = get_thickStart_thickStop(transcript, sequence_for_first_ORF)
                thickStart_longestORF, thickStop_longestORF = get_thickStart_thickStop(transcript, sequence_for_longest_ORF)
                # if transcript['name']=='m54304_190611_190806/61145831/ccs':
                #     print(transcript)
                #     break
                introns_str = ','.join([f"{transcript['chrom']}_{start}_{stop}_{transcript['strand']}" for start, stop in introns])
                if introns_str == '':
                    introns_str = None
                _ = f_annotatedStart.write(f"{transcript['chrom']}\t{transcript['start']}\t{transcript['end']}\t{transcript['name']}\t{transcript['score']}\t{transcript['strand']}\t{thickStart_firstStart}\t{thickStop_firstStart}\t{NMD_classification_to_rgb[firstStart_NMD]}\t{len(transcript['exons'])}\t{','.join(str(e) for e in transcript['block_sizes'])}\t{','.join(str(e) for e in transcript['block_starts'])}\t{sequence_for_first_ORF}\t{firstStart_NMD}\t{all(identifiable_introns)}\t{introns_str}\n")
                _ = f_longestORF.write(f"{transcript['chrom']}\t{transcript['start']}\t{transcript['end']}\t{transcript['name']}\t{transcript['score']}\t{transcript['strand']}\t{thickStart_longestORF}\t{thickStop_longestORF}\t{NMD_classification_to_rgb[longestORF_NMD]}\t{len(transcript['exons'])}\t{','.join(str(e) for e in transcript['block_sizes'])}\t{','.join(str(e) for e in transcript['block_starts'])}\t{sequence_for_longest_ORF}\t{longestORF_NMD}\t{all(identifiable_introns)}\t{introns_str}\n")



