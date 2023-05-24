#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : SumJuncCountsAcrossYRISamples
# @created     : Wednesday Jan 18, 2023 12:13:57 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/regtools_annotate_combined/Comprehensive.ALL SplicingAnalysis/regtools_annotate_combined/Comprehensive.ALL SplicingAnalysis/regtools_annotate_combined/comprehensive.bed.gz SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.Expression.Splicing.bed.gz SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.MetabolicLabelled.30min.bed.gz SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.MetabolicLabelled.60min.bed.gz SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.chRNA.Expression.Splicing.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(data.table)

output_fn <- args[1]
output_longtable_fn <- args[2]
juncslist_fn_in <- args[3]
juncstables_fns <- args[-c(1:3)]

SpliceJunctionCountTables <- setNames(juncstables_fns, str_replace(juncstables_fns, "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.(.+?).bed.gz", "\\1")) %>%
  lapply(fread)

GetSamples <- function(df){
  df %>%
    dplyr::select(-c(1:6)) %>%
    colnames() %>%
    as.data.frame()
}

SamplesInAllDatasets <- lapply(SpliceJunctionCountTables, GetSamples) %>%
  bind_rows(.id = "Dataset") %>%
  dplyr::select(Dataset, Sample=".") %>%
  add_count(Sample) %>%
  filter(n==4)

ReformatJuncTable <- function(df){
  df %>%
    dplyr::select(1:6, SamplesInAllDatasets$Sample)
}

SpliceJunctionCountTables$chRNA.Expression.Splicing %>% ReformatJuncTable() %>% head()
SpliceJunctionCountTables$chRNA.Expression.Splicing  %>% dim()

# Long.table <- lapply(SpliceJunctionCountTables, ReformatJuncTable) %>%
#   lapply(pivot_longer,names_to="Sample", values_to="Count", -c(1:6)) %>%
#   bind_rows(.id="Dataset")

Long.table <- lapply(SpliceJunctionCountTables, pivot_longer,names_to="Sample", values_to="Count", -c(1:6)) %>%
  bind_rows(.id="Dataset")

Long.table %>%
    write_tsv(output_longtable_fn)

JuncSums <- Long.table %>%
    filter(Sample %in% SamplesInAllDatasets) %>%
    group_by(Dataset, junc, .drop=F) %>%
    summarise(JuncSum = sum(Count)) %>%
    pivot_wider(names_from = "Dataset", values_from="JuncSum") %>%
    write_tsv(output_fn)
