#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PreSnakemake_MakeConfigs_v2
# @created     : Monday Jun 07, 2021 10:29:46 CDT
#
# @description : Version 2 of making samples.tsv config file for snakemake.
# Read in the old samples.tsv file, and edit it to account for known sample
# swaps. Also, add columns for GEO accession download and add 4sU samples from
# Li et al
######################################################################

library(tidyverse)

samples <- read_tsv("config/OldSamplesConfig/20210607_samples_Warning_KnownSampleSwampsExist.tsv")

chRNA.SampleSwaps <- read_tsv("../data/20210604_chRNA_SampleIDs_FromBamToFix.txt")

SampleSwapKey <- chRNA.SampleSwaps %>%
    mutate(BestMatchNoNa = case_when(
                                 is.na(BestMatch) ~ ExpectedSampleID,
                                 TRUE ~ BestMatch
                                 )) %>%
    select(ExpectedSampleID, BestMatchNoNa) %>%
    deframe()

chRNA.SampleSwaps %>%
    dplyr::count(BestMatch) %>%
    print(n=Inf)


samples %>%
    filter(Phenotype == "chRNA.Expression.Splicing" & IndID=="NA19130") %>%
    select(R1_local, IndID) %>%
    pull(R1_local)

test <- samples %>%
    mutate(IndID2 = case_when(
                              Phenotype == "chRNA.Expression.Splicing" ~ recode(IndID, !!!SampleSwapKey),
                              TRUE ~ IndID
                              )) %>%
    filter(Phenotype == "chRNA.Expression.Splicing" & IndID2=="NA19130") %>%
    select(R1_local, IndID2)
test$R1_local
test$IndID2
