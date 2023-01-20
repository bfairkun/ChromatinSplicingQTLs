#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : AddQvalueToQtlToolsOutput
# @created     : Monday May 17, 2021 19:36:28 CDT
#
# @description : Add qvalue to QTLtools output
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 " QTLs/QTLTools/chRNA.IR/PermutationPass.txt.gz scratch/Qvals.txt.gz   ", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
FileOut <- args[2]
Pass <- args[3]

library(tidyverse)
library(qvalue)

dat.in <- read_delim(FileIn, delim=' ') %>% na.omit()

if (Pass %in% c('PermutationPass', 'GroupedPermutationPass', 'PermutationPass500kb', 'PermutationPass250kb')){
dat.in$q <- signif(qvalue(dat.in$adj_beta_pval)$qvalues, 5)
} else {
dat.in$q <- signif(qvalue(dat.in$nom_pval)$qvalues, 5)
}


write_delim(dat.in, FileOut, delim=' ')
