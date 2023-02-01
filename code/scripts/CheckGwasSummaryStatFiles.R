#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CheckGwasSummaryStatFiles
# @created     : Tuesday Jan 31, 2023 11:22:48 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "config/gwas_table.tsv", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

gwas.table.fn <- args[1]

library(tidyverse)
library(data.table)

GetNamedFileVector <- function(gwas_accessions.vec, wildcard_path){
    wildcard_path.splt <- strsplit(wildcard_path, "[{}]")
    paste0(wildcard_path.splt[[1]][1], gwas_accessions.vec, wildcard_path.splt[[1]][3]) %>%
        setNames(gwas_accessions.vec) %>%
        return()
}

gwas.table <- read_tsv(gwas.table.fn, comment="#")

gwas.table$gwas

bed_full_wildcard <- "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz"
leadsnps_wildcard <- "gwas_summary_stats/leadSnps/{accession}.bed"
unstandardized_wildcard <- "gwas_summary_stats/StatsForColoc/{accession}.unstandardized.txt.gz"
standardized_wildcard <- "gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz"

bed_full.sampledat <- GetNamedFileVector(gwas.table$gwas, unstandardized_wildcard) %>%
    lapply(read_tsv, n_max=1)
