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
                 "gwas_summary_stats/StatsForColoc/IMSGC2019.unstandardized.txt.gz scratch/standardized.tst.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

f_in <- args[1]
f_out <- args[2]

dat <- read_tsv(f_in)

dat %>%
    select(loci, chrom, start, P, OR, A1, A2) %>%
    # filter(P<0.000001) %>%
    # head(10000) %>%
    mutate(beta = log(OR)) %>%
    mutate(SE=abs(beta/qnorm(P/2))) %>%
    select(loci, chrom, start, beta, SE, A1, A2) %>%
    write_tsv()

