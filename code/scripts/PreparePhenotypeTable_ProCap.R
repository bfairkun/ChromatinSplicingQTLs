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
                 "ProCapAnalysis/ProCap.CountTable.hg38.bed scratch/test.procap.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


FileIn <- args[1]
FileOut <- args[2]

library(tidyverse)
library(RNOmni)

dat <- read_tsv(FileIn) %>%
    distinct(pid, .keep_all=T) %>%
    select(1:6, matches("^GM\\d+$")) %>%
    rename_with(~str_replace(., "GM", "NA"), starts_with("GM"))


dat.cpm <- dat %>%
    select(-c(1:3,5,6)) %>%
    column_to_rownames("pid") %>%
    as.matrix()

dat.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)


Out <- dat %>%
    select(1:6) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("pid")),
               by = "pid") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    filter(`#Chr` %in% paste0("chr", 1:22)) %>%
    arrange(`#Chr`, start)

# Write all samples out
write_tsv(Out, FileOut)

