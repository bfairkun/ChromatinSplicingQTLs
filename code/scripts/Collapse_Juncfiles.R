#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : Collapse_Juncfiles
# @created     : Tuesday Sep 06, 2022 15:01:58 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "my extended list", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

fout <- args[1]
f_in_list <- args[-1]

library(tidyverse)
library(data.table)

dat <- lapply(f_in_list, fread) %>%
    bind_rows()

dat %>%
    separate(V11, into=c("bSize1", "bSize2"), sep=",", convert=T, remove=F) %>%
    mutate(JuncStart = V2 + bSize1) %>%
    mutate(JuncEnd = V3 - bSize2) %>%
    distinct(JuncStart, JuncEnd, .keep_all=T) %>%
    select(-c("bSize1", "bSize2", "JuncStart", "JuncEnd")) %>%
    arrange(V1, V2, V3) %>%
    write_tsv(fout, col_names=F)

