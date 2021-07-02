#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : PlotFeatureSizes_FromPhenotypeTable
# @created     : Tuesday Jun 29, 2021 16:31:32 CDT
#
# @description : Plot feature sizes
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.qqnorm.bed.gz scratch/test.featureSizes.pdf", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

f_in <- args[1]
plot_out <- args[2]

dat <- read_tsv(f_in)

pid_lengths <- dat %>%
    mutate(length = end - start) %>%
    mutate(FeatureType = "pid") %>%
    select(length, FeatureType)

gid_lengths <- dat %>%
    group_by(gid) %>%
    summarize(
              start = min(start),
              end = max(end)
    ) %>%
    mutate(length = end - start) %>%
    mutate(FeatureType = "gid") %>%
    select(length, FeatureType)


p <- bind_rows(pid_lengths, gid_lengths) %>%
    ggplot(aes(x=length)) +
    geom_histogram() +
    scale_x_continuous(trans="log10") +
    facet_wrap(~FeatureType, scales="free_y") +
    theme_classic()

ggsave(plot_out, p, height=3, width=3)
