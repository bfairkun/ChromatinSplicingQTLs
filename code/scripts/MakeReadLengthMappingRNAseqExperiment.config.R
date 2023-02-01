#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : MakeReadLengthMappingRNAseqExperiment.config
# @created     : Monday Jan 30, 2023 12:12:18 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

dat <- read_tsv("config/samples.tsv", comment="#")

dat %>%
    filter(Assay == "RNA-seq") %>%
    filter(!Phenotype == "ProCap") %>%
    filter(RepNumber == 1) %>%
    distinct(Phenotype, RepNumber, IndID) %>%
    add_count(IndID) %>%
    filter(n==4) %>%
    write_tsv("config/RNASeqReadLengthMapExperiment.samples.tsv")
