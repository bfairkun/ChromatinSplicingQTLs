#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : TidyDeeptoolsComputeMatrix
# @created     : Monday Oct 10, 2022 10:41:51 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/compputeMatrix.tidy.txt.gz scratch/Out.tab scratch/out.mat.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(data.table)

f_out <- args[1]
f_in <- args[2]
f_mat_in <- args[3]
f_in_rev <- args[4]


header <- read_tsv(f_in, comment="#", n_max=0)

dat <- read_tsv(f_in, comment="#", skip=3, col_names=colnames(header)[-1])
mat.gz <- fread(f_mat_in, sep='\t', skip=1, select=1:6)

colnames(dat)

dat %>%
  mutate(gene_strand = mat.gz$V6) %>%
  mutate(rn = row_number()) %>%
  gather("sample", "NormalizedCoverage", -rn, -gene_strand) %>%
  separate(sample, into=c("sample", "bin"), convert=T, sep='_') %>%
  mutate(sample=str_replace(sample, "^(.+?)-(.+?)-(.+?)-(.+?)$", "\\1;\\2;\\3;\\4")) %>%
  separate(sample, into=c("dummy", "type", "strand", "genotype"), sep=';') %>%
  replace_na(list(bin=0)) %>%
  group_by(genotype, strand, bin, gene_strand) %>%
  summarise(mean=mean(NormalizedCoverage, na.rm=T)) %>%
  ggplot(aes(x=bin, y=mean, color=genotype)) +
  geom_line() +
  # scale_y_continuous(trans='log10') +
  facet_grid(rows = vars(strand), cols=vars(gene_strand), scales = "free_y") +
  theme_bw()
