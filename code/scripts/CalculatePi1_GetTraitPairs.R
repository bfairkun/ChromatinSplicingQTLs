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
                 "hyprcoloc/Results/ForColoc/MolColocStandard/results.txt.gz hyprcoloc/Results/ForColoc/MolColocStandard/tidy_results_OnlyColocalized.txt.gz hyprcoloc/Results/ForColoc/MolColocStandard/PairwiseColocStats.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_in <- args[1]
f_in_betas <- args[2]
f_out <- args[3]

library(tidyverse)
library(data.table)
library(qvalue)

dat <- read_tsv(f_in)
dat.betas <- read_tsv(f_in_betas)

PermutationPass.dat <- Sys.glob("../code/QTLs/QTLTools/*/PermutationPassForColoc.txt.gz") %>%
  setNames(str_replace(., "../code/QTLs/QTLTools/(.+?)/PermutationPassForColoc.txt.gz", "\\1")) %>%
  lapply(fread, sep=' ') %>%
  bind_rows(.id="Phenotype") %>%
  select(PC=Phenotype, phe_id, p_permutation=adj_beta_pval, singletrait_topvar=var_id, singletrait_topvar_chr = var_chr, singletrait_topvar_pos=var_from) %>%
  group_by(PC) %>%
  mutate(FDR = qvalue(p_permutation)$qvalues) %>%
  mutate(phe_id = str_replace(phe_id, "(.+):(.+)$", "\\1;\\2")) %>%
  separate(phe_id, into=c("P", "GeneLocus"), sep=";") %>%
  unite(Trait, PC, P, sep=";")

dat.tidy <- dat %>%
  left_join(
    dat.betas %>%
      select(GeneLocus=Locus, TopCandidateSNP=snp, Trait=phenotype_full, beta, beta_se, p)
    , by=c("GeneLocus", "TopCandidateSNP", "Trait")) %>%
  left_join(PermutationPass.dat, by=c("GeneLocus", "Trait")) %>%
  left_join(., ., by = "GeneLocus") %>%
    filter(Trait.x != Trait.y) %>%
    separate(Trait.x, into = c("PC1","P1"), sep = "[, ;]", remove=F) %>%
    separate(Trait.y, into = c("PC2","P2"), sep = "[, ;]", remove=F) %>%
    rowwise() %>%
    ungroup() %>%
    mutate(IsColocalizedPair = HyprcolocIteration.x == HyprcolocIteration.y) %>%
    replace_na(list(IsColocalizedPair = FALSE))

dat.tidy %>%
  write_tsv(f_out)



