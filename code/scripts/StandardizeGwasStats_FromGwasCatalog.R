#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : StandardizeGwasStats_MS
# @created     : Monday Apr 18, 2022 12:05:37 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "gwas_summary_stats/StatsForColoc/GCST004599.unstandardized.txt.gz scratch/standardized2.tst.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

f_in <- args[1]
f_out <- args[2]

dat <- read_tsv(f_in)

dat %>%
    select(loci, chrom=hm_chrom, start=hm_pos, beta=hm_beta, SE=standard_error, A1=hm_effect_allele, A2=hm_other_allele) %>%
    mutate(chrom=paste0("chr", chrom)) %>%
    write_tsv(f_out)

