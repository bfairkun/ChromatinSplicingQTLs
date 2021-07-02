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
                                 is.na(BestMatch) ~ paste0(ExpectedSampleID,"Unchecked"),
                                 TRUE ~ BestMatch
                                 )) %>%
    select(ExpectedSampleID, BestMatchNoNa) %>%
    # select(ExpectedSampleID, BestMatch) %>%
    deframe()


# LiEtAlTable <- read_csv("config/ExternalFastqDataAccessions/LiSampleList_SraRunTable_4sU_PRJNA302818.txt") %>%
#     select(incubated_with_4su_for, IndID = cell_line_id, SRA_Run=Run) %>%
#     mutate(
#            Phenotype = paste0("MetabolicLabelled.", incubated_with_4su_for),
#            Assay = "RNA-seq",
#            PairedEnd = FALSE,
#            SRA_Project = "PRJNA302818",
#            RepNumber = 1
#             ) %>%
#     select(Phenotype, Assay, PairedEnd, SRA_Run, SRA_Project, IndID, RepNumber)


LiEtAlTable <- read_tsv("config/ExternalFastqDataAccessions/LiSampleList_ENARunTable_4sU_PRJNA302818.txt") %>%
    select(study_accession, R1_ftp = fastq_ftp, R1_aspera = fastq_aspera, sample_title) %>%
    mutate(IndID.Time = str_replace(sample_title, ".+?\\[(.+?)-(\\d+)-.+", "NA\\1.\\2")) %>%
    separate(IndID.Time, into=c("IndID", "incubated_with_4su_for"), sep="\\.") %>%
    mutate(
           Phenotype = paste0("MetabolicLabelled.", incubated_with_4su_for, "min"),
           Assay = "RNA-seq",
           PairedEnd = FALSE,
           RepNumber = 1
            ) %>%
    select(Phenotype, study_accession,  Assay, PairedEnd, IndID, RepNumber, R1_ftp, R1_aspera)


samples %>%
    mutate(IndID = case_when(
                              Phenotype == "chRNA.Expression.Splicing" ~ recode(IndID, !!!SampleSwapKey),
                              TRUE ~ IndID
                              )) %>%
    mutate(PairedEnd = TRUE) %>%
    bind_rows(LiEtAlTable) %>%
    mutate(Include = T) %>%
    write_tsv("config/samples.tsv")

