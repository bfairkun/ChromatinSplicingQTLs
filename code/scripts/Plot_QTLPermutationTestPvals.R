#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Plot_QTLPermutationTestPvals
# @created     : Tuesday Jul 06, 2021 10:14:26 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args

if(interactive()){
    args <- scan(text=
                 " QTLs/QTLTools/chRNA.IR/PermutationPass.FDR_Added.txt.gz scratch/Pvals.QQ.pdf  scratch/Pvals.hist.pdf ", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
QQPlotOut <- args[2]
HistPlotOut <- args[3]

library(tidyverse)

dat.in <- read_delim(FileIn, delim=' ')


# QQ Plot
PlotOut.QQ <-
 dat.in %>%
     mutate(expected.p = -log10(percent_rank(adj_beta_pval))) %>%
     mutate(observed.p = case_when(
                                   -log10(adj_beta_pval) > 20 ~ 20,
                                   TRUE ~ -log10(adj_beta_pval)
                  )) %>%
     ggplot(aes(x=expected.p, y=observed.p)) +
     geom_point() +
     geom_abline() +
     xlab("-log10(Expected P)") +
     ylab("-log10(Observed P)") +
     theme_classic()
ggsave(QQPlotOut, PlotOut.QQ, height=3, width=3)

PlotOut.P <- dat.in %>%
    ggplot(aes(x=adj_beta_pval)) +
    geom_histogram(breaks=seq(0,1,0.02)) +
    ylab("Frequency") +
    xlab("P-value\nfeature-level permutation test") +
    theme_classic()
ggsave(HistPlotOut, PlotOut.P, height=3, width=3)


