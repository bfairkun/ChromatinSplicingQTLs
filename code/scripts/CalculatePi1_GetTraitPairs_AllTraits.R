#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CalculatePi1_GetTraitPairs
# @created     : Thursday Jun 23, 2022 10:50:15 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "10 ../code/scratch/PairwisePi1Traits.", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

Num_f_out_chunks <- args[1]
f_out_prefix <- args[2]
permutation_f_in <- args[-c(1:2)]

library(tidyverse)
# library(data.table)
library(qvalue)


PermutationPass.dat <- permutation_f_in %>%
  setNames(str_replace(., "QTLs/QTLTools/(.+?)/PermutationPassForColoc.txt.gz", "\\1")) %>%
  lapply(read_delim, delim=' ') %>%
  bind_rows(.id="Phenotype") %>%
  select(PC=Phenotype, phe_id, p_permutation=adj_beta_pval, singletrait_topvar=var_id, singletrait_topvar_chr = var_chr, singletrait_topvar_pos=var_from) %>%
  group_by(PC) %>%
  mutate(FDR = qvalue(p_permutation)$qvalues) %>%
  mutate(phe_id = str_replace(phe_id, "(.+):(.+)$", "\\1;\\2")) %>%
  separate(phe_id, into=c("P", "GeneLocus"), sep=";") %>%
  unite(Trait, PC, P, sep=";")

dat.pairs <- PermutationPass.dat %>%
  filter(FDR<0.1) %>%
  left_join(., PermutationPass.dat, by="GeneLocus") %>%
  filter(!Trait.x==Trait.y) %>%
  separate(Trait.x, into=c("PC1", "P1"), sep=";") %>%
  separate(Trait.y, into=c("PC2", "P2"), sep=";") %>%
  group_by(GeneLocus) %>%
  mutate(Random = sample(10,1))

dat.pairs %>%
  group_by(Random) %>%
  group_walk(~ write_tsv(.x, paste0(f_out_prefix, .y$Random, ".txt.gz")))



