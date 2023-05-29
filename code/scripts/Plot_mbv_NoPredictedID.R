#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : Plot_mbv
# @created     : Thursday Jun 03, 2021 13:59:02 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QC/mbvLongReads/data/GM10.txt scratch/test.mbv.pdf", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)

f_in <- args[1]
plot_out <- args[2]
expectedSampleID <- args[3]

dat <- read_delim(f_in, ' ')

# best_het <- dat %>%
#     filter(perc_het_consistent == max(perc_het_consistent)) %>% pull(SampleID)
# best_hom <- dat %>%
#     filter(perc_het_consistent == max(perc_het_consistent)) %>% pull(SampleID)

# p<- dat %>%
#     mutate(IsSample=SampleID==expectedSampleID) %>%
#     arrange(IsSample) %>%
#     ggplot(aes(x=perc_het_consistent, y=perc_hom_consistent, color=IsSample, label=SampleID)) +
#     geom_text(size=1.5) +
#     theme_minimal() +
#     theme(legend.position="none")

p<- dat %>%
    mutate(IsSample=SampleID==expectedSampleID) %>%
    arrange(IsSample) %>%
    ggplot(aes(x=perc_het_consistent, y=perc_hom_consistent, label=SampleID)) +
    geom_text(size=1.5) +
    theme_minimal() +
    theme(legend.position="none")

ggsave(plot_out, p, height=3, width=3)

