#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : AnnotateIntrons
# @created     : Tuesday Dec 06, 2022 16:41:39 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz SplicingAnalysis/regtools_annotate_combined/basic.bed.gz ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf_all_introns.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


numers_f_in <- args[1]
all_ints_f_in <- args[2]
leafviz_f_in <- args[3]

library(tidyverse)
library(data.table)

numers <- fread(numers_f_in, select=1, sep=' ')
introns <- fread(all_ints_f_in, sep='\t')
introns.leafvizscript <- fread(leafviz_f_in, sep='\t', col.names=c("chrom", "start", "stop", "gene", "gene_id","strand","transcript", "intron_num", "transcript_tag", "tag"))

IntronAnnotations <- introns.leafvizscript %>%
  group_by(chrom, start, stop, strand, gene_id) %>%
  summarise(Annotation = case_when(
    all(transcript_tag=="nonsense_mediated_decay") ~ "Unique to nonsense_mediated_decay",
    all(transcript_tag=="non_stop_decay") ~ "Unique to non_stop_decay",
    all(transcript_tag=="processed_transcript") ~ "Unique to processed_transcript",
    all(transcript_tag=="retained_intron") ~ "Unique to retained_intron",
    any(transcript_tag=="protein_coding") ~ "In protein_coding",
    TRUE ~ "Other"
  )) %>%
  ungroup()


