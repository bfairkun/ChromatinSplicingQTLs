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
                 "featureCounts/MetabolicLabelled.30min/Counts.txt QTLs/QTLTools/MetabolicLabelled.30min/OnlyFirstReps.sorted.qqnorm.bed.gz scratch/TestRPKM.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

counts.fn <- args[1]
standardized.fn <- args[2]
out.fn <- args[3]

library(tidyverse)
library(edgeR)

read_and_rename <- function(fn){
  read_tsv(fn, comment="#") %>%
    dplyr::select(Geneid, Length, contains("Alignments")) %>%
    rename_at(vars(contains("Alignments")), ~str_replace(., "Alignments/STAR_Align/.+?/(.+?)/1/Filtered.bam", "\\1")) %>%
    return()
}

standardized <- read_tsv(standardized.fn)

counts <- read_and_rename(counts.fn) %>%
  dplyr::select(Geneid, Length, !matches(".+?_\\d+")) %>%
  filter(Geneid %in% standardized$pid)

rpkm <- counts %>%
  dplyr::select(-Length) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  rpkm(prior.count=0.1, log=T, gene.length=counts$Length) %>%
  as.data.frame() %>%
  rownames_to_column("pid")

standardized %>%
  dplyr::select(1:6) %>%
  inner_join(
    rpkm
  ) %>%
  write_tsv(out.fn)

