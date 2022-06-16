#!/usr/bin/env bash

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : BamToBigwig
# @created     : Monday Jul 13, 2020 14:06:48 CDT
#
# @description : Convert BAM to bigwig raw coverage
######################################################################

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0`  <in.fa.fai> <in.bam> <out.bw> [KEY=VALUE]"
  echo "Parameters:"
  echo "<in.fa.fai> Required. Fasta index file containing chrome sizes. See samtools faidx"
  echo "<in.bam> Required. position sorted bam file"
  echo "<out.bw> Required. output bigwig file"
  echo "[KEY=VALUE] Optional. Other named args. Possible KEYs include SORT_ARGS, SAMTOOLSVIEW_ARGS, GENOMECOV_ARGS, MKTEMP_ARGS, REGION, and bw_minus. VALUES that contain spaces or other weird symbols should be quoted to avoid bash expanding things"
  echo "If bw_minus is supplied, then the <out.bw> will be for + strand coverage and bw_minus VALUE output file will be for - strand coverage"
  echo "REGION can be particularly useful for outputing just region; example: REGION=chr11:65,495,143-65,509,111"
  echo "GENOMECOV_ARGS can be particularly useful for adding parameters like -split to the bedtools genomecov command for RNA-seq spliced read coverage (GENOMECOV_ARGS=-split)"
  echo ""
  echo "This program will create bigwig file from sorted indexed bam alignments. Makes use of bedtools genomecov (https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) to create an temporary intermediate bedgraph file. This bedtools genomecov step has some useful parameters, which can be specified in this script with the optional [GENOMECOV_ARGS=VALUE] argument. For example, set VALUE to '-split' for RNA-seq to ignore intron alignment gaps, and use '-pc' to fill in paired-end read gaps from ChIP-seq data"
  exit 0
fi

set -xe

fai=$1
bam=$2
bw=$3
shift 3

# Parse optional named args
## Define defaults for optional named args
# Usefule for filtering bam with some command
export EXTRA_BAM_FILTER_PIPE=""
# Useful for adding sort options: SORT_ARGS="-T /path/to/different/tmp/dir""
export SORT_ARGS=""
# Useful for setting samtools view flags to filter for only R1 or something: SAMTOOLSVIEW_ARGS=-f64. For stranded mode these will be applied after merge stranded bams
export SAMTOOLSVIEW_ARGS="-F256"
# Useful for setting bedtools genomecov flag to pileup split alignments or something: GENOMECOV_ARGS=-split
export GENOMECOV_ARGS=""
# Useful for setting mktemp parameters: MKTEMP_ARGS="-p /scratch/midway2/bjf79"
export MKTEMP_ARGS=""
# Used in samtools command to only output reads in a region
export REGION=""
# If supplied, then $bw output will be plus strand and $bw_minus output will be minus strand
export bw_minus=""
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)
   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"
   export "$KEY"="$VALUE"
done

if [[ $bw_minus == "" ]]
then
    temp_bgfile=$(mktemp $MKTEMP_ARGS)
    samtools view -h $SAMTOOLSVIEW_ARGS $bam $REGION | perl -lane 'print if (( $F[0] =~ /^@/ ) || ( abs($F[8]) <= 10000  ))' | samtools view -bh | bedtools genomecov -ibam - -bga $GENOMECOV_ARGS | LC_COLLATE=C sort $SORT_ARGS -k1,1 -k2,2n > $temp_bgfile
    bedGraphToBigWig $temp_bgfile $fai $bw
    rm ${temp_bgfile}
else
    temp_bamfile_plus=$(mktemp $MKTEMP_ARGS)
    temp_bamfile_minus=$(mktemp $MKTEMP_ARGS)
    temp_bgfile_plus=$(mktemp $MKTEMP_ARGS)
    temp_bgfile_minus=$(mktemp $MKTEMP_ARGS)

    # Separate reads to two bams based on strand of R1
    samtools merge -f $temp_bamfile_plus <(samtools view -bh -F144 $bam $REGION) <(samtools view -bh -f146 $bam $REGION)
    samtools merge -f $temp_bamfile_minus <(samtools view -bh -F128 -f16 $bam $REGION) <(samtools view -bh -f162  $bam $REGION)

    # Create bedgraphs
    samtools view -h $SAMTOOLSVIEW_ARGS $temp_bamfile_plus | perl -lane 'print if (( $F[0] =~ /^@/ ) || ( abs($F[8]) <= 10000  ))' | samtools view -bh | bedtools genomecov -ibam - -bga $GENOMECOV_ARGS |  LC_COLLATE=C sort $SORT_ARGS -k1,1 -k2,2n > $temp_bgfile_plus
    samtools view -h $SAMTOOLSVIEW_ARGS $temp_bamfile_minus | perl -lane 'print if (( $F[0] =~ /^@/ ) || ( abs($F[8]) <= 10000  ))' | samtools view -bh | bedtools genomecov -ibam - -bga $GENOMECOV_ARGS | LC_COLLATE=C sort $SORT_ARGS -k1,1 -k2,2n > $temp_bgfile_minus

    #Create bigwigs
    bedGraphToBigWig $temp_bgfile_plus $fai $bw
    bedGraphToBigWig $temp_bgfile_minus $fai $bw_minus

    rm ${temp_bamfile_plus}
    rm ${temp_bamfile_minus}
    rm ${temp_bgfile_plus}
    rm ${temp_bgfile_minus}
fi
