#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : GetGeneSetsForQTLSNP_Enrichment
# @created     : Monday Jul 10, 2023 13:10:24 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                   "../code/hyprcoloc/Results/ForColoc/MolColocNonRedundantFullSplicing/tidy_results_OnlyColocalized.txt.gz scratch/ContrastingGeneSets.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}



library(tidyverse)

f_in <- args[1]
f_out <- args[2]

coloc.dat <- read_tsv(f_in) %>%
  separate(phenotype_full, into=c("PC", "P"), sep=";")

coloc.dat$PC %>% unique()

eGenesWith_molQTL_Coloc <-  coloc.dat %>%
  group_by(snp, Locus) %>%
  filter(any(str_detect(PC, "^Expression.Splicing"))) %>%
  ungroup() %>%
  filter(str_detect(PC, "^Expression.Splicing") & P==Locus) %>%
  pull(P) %>% unique()

eGenesWith_hQTL_Coloc <- coloc.dat %>%
  group_by(snp, Locus) %>%
  filter(any(str_detect(PC, "^Expression.Splicing"))) %>%
  filter(any(PC %in% c("H3K27AC", "H3K4ME3", "H3K4ME1", "H3K36ME3", "ProCap"))) %>%
  ungroup() %>%
  filter(str_detect(PC, "^Expression.Splicing") & P==Locus) %>%
  pull(P) %>% unique()


eGenesWith_Non_hQTL_Coloc <- setdiff(eGenesWith_molQTL_Coloc, eGenesWith_hQTL_Coloc)

bind_rows(
  hQTL_eGenes = data.frame(genes = eGenesWith_hQTL_Coloc),
  Other_eGenes = data.frame(genes = eGenesWith_Non_hQTL_Coloc),
  .id="Class"
) %>%
  write_tsv(f_out)

