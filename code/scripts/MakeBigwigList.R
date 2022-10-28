#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : MakeBigwigList
# @created     : Tuesday Oct 11, 2022 15:14:06 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)

samples <- read_tsv("config/samples.tsv", comment='#')

samples %>%
    count(Phenotype)

RNASeqPhenotypes <- samples %>%
    filter(Assay=="RNA-seq") %>%
    filter(!Phenotype=="ProCap") %>%
    distinct(Phenotype) %>% pull(Phenotype)

igsr_samples <- read_tsv("../data/igsr_samples.tsv.gz")

YRI.samples <- igsr_samples %>%
    filter(`Population code`=="YRI") %>%
    pull("Sample name")

RecodePhenotypeNames <- c("chRNA.Expression.Splicing"="chRNA",
                        "Expression.Splicing"="polyA.RNA",
                        "MetabolicLabelled.30min"="4sU.30min.RNA",
                        "MetabolicLabelled.60min"="4sU.60min.RNA")

samples.to.include <- samples %>%
    filter(!Phenotype %in% c("ProCap")) %>%
    filter(RepNumber==1) %>%
    filter( (IndID %in% YRI.samples) | (Phenotype == "CTCF") ) %>%
    distinct(IndID, RepNumber, Phenotype)

samples.to.include %>%
    mutate(BigwigFilepath = str_glue("bigwigs/{Phenotype}/{IndID}.{RepNumber}.bw")) %>%
    mutate(Group_label = recode(Phenotype, !!!RecodePhenotypeNames)) %>%
    mutate(Strand = ".") %>%
    dplyr::select(SampleID=IndID, BigwigFilepath, Group_label, Strand) %>%
    write_tsv("Metaplots/bwList.allunstranded.tsv")

bind_rows(
        samples.to.include %>%
            filter(!Phenotype == "chRNA.Expression.Splicing") %>%
            mutate(BigwigFilepath = str_glue("bigwigs/{Phenotype}/{IndID}.{RepNumber}.bw")) %>%
            mutate(Strand = "."),
        samples.to.include %>%
            filter(Phenotype == "chRNA.Expression.Splicing") %>%
            mutate(Strand = "-") %>%
            mutate(BigwigFilepath = str_glue("bigwigs/chRNA.Expression.Splicing_stranded/{IndID}.{RepNumber}.plus.bw")),
        samples.to.include %>%
            filter(Phenotype == "chRNA.Expression.Splicing") %>%
            mutate(Strand = "+") %>%
            mutate(BigwigFilepath = str_glue("bigwigs/chRNA.Expression.Splicing_stranded/{IndID}.{RepNumber}.minus.bw"))) %>%
  mutate(Group_label = recode(Phenotype, !!!RecodePhenotypeNames)) %>%
    dplyr::select(SampleID=IndID, BigwigFilepath, Group_label, Strand) %>%
    write_tsv("Metaplots/bwList.tsv")

colors.xl <- readxl::read_excel("../data/ColorsForPhenotypes.xlsx", sheet=2)

colors.xl %>%
    mutate(BedgzFile = case_when(
                                 ParentAssay %in% RNASeqPhenotypes ~ str_glue("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{ParentAssay}.bed.gz"),
                                 TRUE ~ ""
                                 )) %>%
    mutate(Group_label = recode(ParentAssay, !!!RecodePhenotypeNames)) %>%
    mutate(Supergroup = Group_label) %>%
    filter(!Group_label=="ProCap") %>%
    dplyr::select(Group_label, Group_color=Hex, BedgzFile, Supergroup) %>%
    arrange(factor(Group_label, levels=c(
                    "CTCF",
                    "H3K4ME1",
                    "H3K27AC",
                    "H3K4ME3",
                    "ProCap",
                    "H3K36ME3",
                    "chRNA",
                    "4sU.30min.RNA",
                    "4sU.60min.RNA",
                    "polyA.RNA"
                    ))) %>%
    write_tsv("Metaplots/bwGroups.tsv")
