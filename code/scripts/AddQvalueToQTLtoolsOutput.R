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
                 "QTLs/QTLTools/Expression.Splicing/PermutationPass.txt.gz QTLs/QTLTools/Expression.Splicing/PermutationPass.QvalsAdded.txt.gz QTLs/QTLTools/Expression.Splicing/PermutationPass.PermutationTestPvalsHist.pdf", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
FileOut <- args[2]
PlotOut <- args[3]

library(tidyverse)
library(qvalue)

dat.in <- read_delim(FileIn, delim=' ')

#QQ Plot
# dat.in %>%
#     mutate(expected.p = -log10(percent_rank(adj_beta_pval))) %>%
#     mutate(observed.p = case_when(
#                                   -log10(adj_beta_pval) > 20 ~ 20,
#                                   TRUE ~ -log10(adj_beta_pval)
#                  )) %>%
#     ggplot(aes(x=expected.p, y=observed.p)) +
#     geom_point() +
#     geom_abline() +
#     xlab("-log10(Expected P)") +
#     ylab("-log10(Observed P)") +
#     theme_classic()

PlotOut.P <- dat.in %>%
    ggplot(aes(x=adj_beta_pval)) +
    geom_histogram() +
    ylab("Frequency") +
    xlab("P-value\nfeature-level permutation test") +
    theme_classic()
ggsave(PlotOut, PlotOut.P, height=3, width=3)

dat.in$q <- signif(qvalue(dat.in$adj_beta_pval)$qvalues, 5)

write_delim(dat.in, FileOut, delim=' ')
