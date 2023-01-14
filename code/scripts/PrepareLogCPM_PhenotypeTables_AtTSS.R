#!/usr/bin/env Rscript

######################################################################
# @author      : cnajar (cnajar@midway2-login1.rcc.local)
# @file        : PrepareCPMTablesFromFeatureCountsForChromatinAtTSS
# @created     : Friday January 1, 2023 15:39:08 CDT
#
# @description : from featureCountsAtRegion (AtTSS) output, prepare CPM 
# phenotype tables for chromatin marks covering the TSS of expressed genes.
# The indended function for these tables is for making plots.
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/AtTSS/H3K27AC/Counts.txt ExpressionAnalysis/polyA/ExpressedGeneList.txt QTLs/QTLTools/H3K27AC/OnlyFirstRepsUnstandardized_AtTSS.qqnorm.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

counts.fn <- args[1]
Genes_bed_f_in <- args[2]
out.fn <- args[3]

library(tidyverse)
library(edgeR)


read_and_rename <- function(fn){
  read_tsv(fn, comment="#") %>%
    dplyr::select(Geneid, contains("Alignments")) %>%
    rename_at(vars(contains("Alignments")), ~str_replace(., "Alignments/Hisat2_Align/.+?/(.+?).1.wasp_filterd.markdup.sorted.bam", "\\1")) %>%
    return()
}

gene.list <- read_tsv(Genes_bed_f_in, col_names=c("Chr", "Start", "End", "Geneid", "score", "Strand"))

counts <- read_and_rename(counts.fn) %>%
  filter(Geneid %in% gene.list$Geneid)

cpm <- counts %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  cpm(prior.count=0.1, log=T) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid")

#standardized %>%
#  dplyr::select(1:6) %>%
#  inner_join(
#    cpm
#  ) %>%
#  write_tsv(out.fn)


gene.list %>%
    dplyr::select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               cpm,
               by = "Geneid"
    ) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start) %>%
    write_tsv(out.fn)

