#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : SplitLeafcutterPerindCounts
# @created     : Tuesday May 18, 2021 17:05:57 CDT
#
# @description : Split leafcutter perind.counts.gz file based on sample names
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

file_in <- args[1] #FileiIn, also prefix for output file(s)

library(tidyverse)

dat.in <- read_delim(file_in, delim=' ')

Phenotypes <- colnames(dat.in[,-1]) %>% str_replace("(.+?)_.+", "\\1") %>% unique()

for (i in seq_along(Phenotypes)){
    print(paste("Processing", i))
    dat.in %>%
        select(chrom, starts_with(paste0(Phenotypes[i], "_"))) %>%
        # rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
        select(chrom, matches("\\_1$")) %>%
        rename_with(~str_replace(., '.+?_(.+?)_1', '\\1')) %>%
        write_delim(delim=' ', paste0(file_in, ".", Phenotypes[i], ".gz"))
}
