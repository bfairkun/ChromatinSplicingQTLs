#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : GetSitesFromPermutationPass
# @created     : Wednesday Oct 19, 2022 21:33:53 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTLs/QTLTools/Expression.Splicing.Subset_YRI/PermutationPass.txt.gz scratch/Test.sites", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

read_delim(args[1], delim=' ') %>%
    dplyr::select(var_id) %>%
    write_tsv(args[2], col_names=F)
