#!/usr/bin/env sh

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : BamToBigwig
# @created     : Monday Jul 13, 2020 14:06:48 CDT
#
# @description : Convert BAM to bigwig raw coverage
######################################################################

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0`  <in.fa.fai> <in.bam> <out.bw> <Other bedtools genomcov parameters>"
  echo "Parameters:"
  echo "<in.fa.fai> Required. Fasta index file containing chrome sizes. See samtools faidx"
  echo "<in.bam> Required. position sorted bam file"
  echo "<out.bw> Required. output bigwig file"
  echo "<Other bedtools genomecov paramaters> Optional. Other paramaters to pass to bedtools genomecov"
  echo ""
  echo "This program will create bigwig file from bam alignments. Makes use of bedtools genomecov (https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) to create an temporary intermediate bedgraph file. This bedtools genomecov step has some useful parameters, which can be specified in this script with the  <Other bedtools genomecov parameters>argument. For example, use -split for RNA-seq to ignore intron alignment gaps, and use -pc to fill in paired-end read gaps from ChIP-seq data"
  exit 0
fi

set -xe

fai=$1
bam=$2
bw_plus=$3
bw_minus=$4
shift 4
ExtraFlags=$*

temp_bamfile_plus=$(mktemp)
temp_bamfile_minus=$(mktemp)
temp_bgfile_plus=$(mktemp)
temp_bgfile_minus=$(mktemp)

# Separate reads to two bams based on strand of R1
samtools merge -f $temp_bamfile_plus <(samtools view -bh -F144 $bam) <(samtools view -bh -f146 $bam)
samtools merge -f $temp_bamfile_minus <(samtools view -bh -F128 -f16 $bam) <(samtools view -bh -f162  $bam)

# Create bedgraphs
samtools view -bh -F256 $temp_bamfile_plus | bedtools genomecov -ibam - -bga $ExtraFlags | bedtools sort -i - > $temp_bgfile_plus
samtools view -bh -F256 $temp_bamfile_minus | bedtools genomecov -ibam - -bga $ExtraFlags | bedtools sort -i - > $temp_bgfile_minus

#Create bigwigs
bedGraphToBigWig $temp_bgfile_plus $fai $bw_plus
bedGraphToBigWig $temp_bgfile_minus $fai $bw_minus

rm ${temp_bamfile_plus}
rm ${temp_bamfile_minus}
rm ${temp_bgfile_plus}
rm ${temp_bgfile_minus}
