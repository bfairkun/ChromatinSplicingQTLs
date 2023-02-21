#!/usr/bin/env Rscript

######################################################################
# @author      : cnajar (cnajar@midway2-login2.rcc.local)
# @file        : Normalize_ReadLengthMapExperimentResults_Counts
# @created     : Saturday Feb 20, 2023 12:18:49 CDT
#
# @description :
######################################################################

library(tidyverse)
library(edgeR)

counts.fn <- 'ReadLengthMapExperimentResults/featureCounts/Counts.txt'
standardized.fn <- 'QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.sorted.qqnorm.bed.gz'
out.fn <- 'ReadLengthMapExperimentResults/tables/AllRNASeq.Normalized.RPKM.bed.gz'

read_and_rename <- function(fn){
  read_tsv(fn, comment="#") %>%
    dplyr::select(Geneid, Length, contains("ReadLengthMapExperiment")) %>%
    rename_at(vars(contains("chRNA.Expression.Splicing")), ~str_replace(., 
                "ReadLengthMapExperiment/chRNA.Expression.Splicing/(.+?)/1/Filtered.bam", "chRNA.\\1")) %>%
    rename_at(vars(contains("Expression.Splicing")), ~str_replace(., 
                "ReadLengthMapExperiment/Expression.Splicing/(.+?)/1/Filtered.bam", "polyA.\\1")) %>%
    rename_at(vars(contains("MetabolicLabelled.30min")), ~str_replace(., 
                "ReadLengthMapExperiment/MetabolicLabelled.30min/(.+?)/1/Filtered.bam", "RNA30min.\\1")) %>%
    rename_at(vars(contains("MetabolicLabelled.60min")), ~str_replace(., 
                "ReadLengthMapExperiment/MetabolicLabelled.60min/(.+?)/1/Filtered.bam", "RNA60min.\\1")) %>%
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