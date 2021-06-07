#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PlotPCA_FromPhenotypeTable
# @created     : Thursday Jun 03, 2021 11:46:58 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz scratch/test.pca.pdf", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)

f_in <- args[1]
plot_out <- args[2]

dat <- read_tsv(f_in)

pca <- dat %>%
    select(-c(`#Chr`, 'start', 'end', 'gid', 'strand')) %>%
    column_to_rownames('pid') %>%
    as.matrix() %>%
    t() %>%
    prcomp() %>%
    extract2("x")

pca.to.plot <- pca %>%
    as.data.frame() %>%
    select(PC1:PC4) %>%
    rownames_to_column("sample")

p <- pca.to.plot %>%
    ggplot(aes(x=PC1, y=PC2, label=sample)) +
    geom_text(size=5)
ggsave(plot_out, p, height=3, width=3)
