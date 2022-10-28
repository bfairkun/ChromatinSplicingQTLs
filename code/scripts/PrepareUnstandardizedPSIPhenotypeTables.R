#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : PrepareUnstandardizedPSIPhenotypeTables
# @created     : Tuesday Oct 18, 2022 16:02:36 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.Expression.Splicing.bed.gz QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.sorted.qqnorm.bed.gz scratch/PSI.table.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

PSI_f.in <- args[1]
standardized_f.in <- args[2]
f.out <- args[3]

standardized <- read_tsv(standardized_f.in)

psi <- read_tsv(PSI_f.in)


standardized %>%
    dplyr::select(1:6) %>%
    inner_join(
               psi %>% dplyr::select(-c(4:5)) %>% dplyr::rename("#Chr"="#Chrom")
    ) %>%
    dplyr::select(1:3, pid, gid, strand, everything()) %>%
    arrange(`#Chr`, start, end) %>%
    write_tsv(f.out,na="NA" )


