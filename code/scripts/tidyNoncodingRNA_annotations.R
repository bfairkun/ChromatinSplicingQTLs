#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : tidyNoncodingRNA_annotations
# @created     : Friday Oct 07, 2022 11:54:49 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "NonCodingRNA_annotation/annotation/ncRNA.bed.gz NonCodingRNA_annotation/annotation/ncRNA.annotation.tab.gz NonCodingRNA_annotation/annotation/ncRNA.categorized.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

bedin_fn = args[1]
annotations_fn = args[2]
bedout_fn = args[3]

library(tidyverse)

bed <- read_tsv(bedin_fn, col_names=c("chrom", "start", "stop", "gene_name", "score", "strand"))
annotations <- read_tsv(annotations_fn) %>%
    dplyr::select(gene_name, rna_type) %>%
    separate_rows(rna_type, sep=",")

bed %>%
    inner_join(annotations) %>%
    mutate(rna_type=paste0("ncRNA_", rna_type)) %>%
    dplyr::select(1:3, rna_type, gene_name, strand) %>%
    arrange(chrom, start, stop) %>%
    write_tsv(bedout_fn, col_names=F)
