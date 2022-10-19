#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Subset_YRI
# @created     : Wednesday Oct 19, 2022 09:35:22 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTLs/QTLTools/polyA.Splicing/OnlyFirstRepsUnstandardized.qqnorm.bed.gz scratch/YRI.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f.in <- args[1]
f.out <- args[2]

library(tidyverse)

igsr <- read_tsv("../data/igsr_samples.tsv.gz")

YRI <- igsr %>% filter(`Population code` == "YRI") %>% pull(`Sample name`)

dat <- read_tsv(f.in)

dat %>% dplyr::select(1:6, intersect(YRI, colnames(dat))) %>%
    write_tsv(f.out, na="NA")
