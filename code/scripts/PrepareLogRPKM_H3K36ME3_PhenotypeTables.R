#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : PrepareLogRPKM_PhenotypeTables
# @created     : Saturday Oct 15, 2022 15:19:49 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "MiscCountTables/H3K36ME3.bed QTLs/QTLTools/H3K36ME3/OnlyFirstReps.sorted.qqnorm.bed.gz scratch/TestRPKM_2.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

counts.fn <- args[1]
standardized.fn <- args[2]
out.fn <- args[3]

library(tidyverse)
library(edgeR)



standardized <- read_tsv(standardized.fn)

counts <- read_tsv(counts.fn) %>%
  mutate(Length = end-start)

rpkm <- counts %>%
  column_to_rownames("pid") %>%
  dplyr::select(-c(1:5), -Length) %>%
  DGEList() %>%
  calcNormFactors() %>%
  rpkm(prior.count=0.1, log=T, gene.length=counts$Length) %>%
  as.data.frame() %>%
  rownames_to_column("pid")

standardized <- standardized[standardized$pid %in% rpkm$pid,]

standardized %>%
  dplyr::select(1:6) %>%
  inner_join(
    rpkm
  ) %>%
  write_tsv(out.fn)

