#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : CombinePermutationPassSNPsIntoBed
# @created     : Wednesday Apr 05, 2023 14:25:25 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/test.aggregate.bed.gz QTLs/QTLTools/APA_Nuclear/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/APA_Total/PermutationPass.FDR_Added.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(data.table)


dat <- args[-1] %>%
    setNames(str_replace(., "QTLs/QTLTools/(.+?)/PermutationPass.FDR_Added.txt.gz", "\\1")) %>%
    lapply(fread) %>%
    bind_rows(.id="PC") %>%
    dplyr::select(phe_id, PC, var_chr, var_from,var_id, nom_pval, q, slope, slope_se) %>%
    mutate(var_to = var_from + 1,
           strand = ".",
           name = paste(PC, phe_id, sep=";")) %>%
    filter(!is.na(var_chr)) %>%
    dplyr::select(var_chr, var_from, var_to, name, nom_pval, strand, slope, slope_se, q, var_id) %>%
    arrange(var_chr, var_from) %>%
    write_tsv(args[1], col_names=F)
