#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Plot_QTLPermutationTestNumSig
# @created     : Tuesday Jul 06, 2021 10:37:29 CDT
#
# @description : Plot number gene level hits at varying FDR thresholds
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/Test.Out.pdf QTLs/QTLTools/CTCF/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/H3K4ME1/PermutationPass.FDR_Added.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

PermutationPassFDR.Added.Files <- Sys.glob("QTLs/QTLTools/*/PermutationPass.FDR_Added.txt.gz")
PlotOut_f <- args[1]
# PermutationPassFDR.Added.Files <- args[-1]

library(tidyverse)

dat <- sapply(PermutationPassFDR.Added.Files, read_delim, simplify=FALSE, delim=' ') %>%
    bind_rows(.id = "id") %>%
    mutate(id = str_replace(id, "QTLs/QTLTools/(.+?)/PermutationPass.FDR_Added.txt.gz", "\\1")) %>%
    select(id, q)

feats.tested <- dat %>%
    filter(!is.na(q)) %>%
    count(id) %>%
    mutate(NewLabel = paste0(id, "\n(", n, " feats tested)"))

FDR_HitsTable <-
    dat %>%
    mutate(
        `0.1`=q<0.1,
        `0.05`=q<0.05,
        `0.01`=q<0.01) %>%
    select(-q) %>%
    gather("Threshold", "IsSig", -id) %>%
    group_by(Threshold, id) %>%
    summarize(NumSignificant = sum(IsSig, na.rm=T)) %>%
    left_join(feats.tested, by="id")

PlotOut <- ggplot(FDR_HitsTable, aes(x=Threshold, y=NumSignificant, fill=NewLabel)) +
    geom_col() +
    # scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    facet_wrap(~NewLabel, scales="free_y") +
    xlab("FDR cutoff") +
    theme_classic() +
    theme(legend.position="none")

ggsave(PlotOut_f, PlotOut, height=6, width=6)
