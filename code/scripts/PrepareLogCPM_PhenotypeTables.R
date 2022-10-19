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
                 "featureCounts/H3K27AC/Counts.txt QTLs/QTLTools/H3K27AC/OnlyFirstReps.sorted.qqnorm.bed.gz scratch/TestCPM.bed.gz", what='character')
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
    dplyr::select(Geneid, contains("Alignments")) %>%
    rename_at(vars(contains("Alignments")), ~str_replace(., "Alignments/Hisat2_Align/.+?/(.+?).1.wasp_filterd.markdup.sorted.bam", "\\1")) %>%
    return()
}

standardized <- read_tsv(standardized.fn)

counts <- read_and_rename(counts.fn) %>%
  filter(Geneid %in% standardized$pid)

cpm <- counts %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  cpm(prior.count=0.1, log=T) %>%
  as.data.frame() %>%
  rownames_to_column("pid")

standardized %>%
  dplyr::select(1:6) %>%
  inner_join(
    cpm
  ) %>%
  write_tsv(out.fn)

