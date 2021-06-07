#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PreparePhenotypeTable_ProCap
# @created     : Tuesday Jun 01, 2021 13:03:37 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "ProCapAnalysis/SupplementDat2.csv scratch/test.procap.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


FileIn <- args[1]
FileOut <- args[2]

library(tidyverse)

dat <- read_csv(FileIn)

dat %>%
    filter(allele_mappable & variable_expression) %>%
    mutate(pid=paste(chromosome, position, dominant_strand, refseq_classification,  sep="."),
           start = floor(position-nTSS_distance/2),
           end = ceiling(position+nTSS_distance/2)) %>%
    mutate(gid=pid) %>%
    select(`#Chr`=chromosome, start, end, pid, gid, strand=dominant_strand, matches("^GM\\d+$")) %>%
    rename_with(~str_replace(., "GM", "NA"), starts_with("GM")) %>%
    write_tsv(FileOut)


