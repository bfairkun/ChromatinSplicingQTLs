#!/usr/bin/env Rscript

if(interactive()){
  args <- scan(text=
                 paste("featureCounts/AtInternalExons/chRNA.Expression/Counts.txt", "featureCounts/chRNA.Expression/Counts.txt"), what='character')
} else{
  args <- commandArgs(trailingOnly=TRUE)
}

ExonExpressionFileIn <- args[1]
GeneExpressionFileIn <- args[2]
Pass <- args[3]


library(tidyverse)
library(edgeR)

exons <- read_tsv(ExonExpressionFileIn, comment = "#") %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "Alignments/STAR_Align/chRNA.Expression.Splicing/(.+?)/(.+?)/Filtered.bam", "\\1.\\2")) %>%
  select(1:6, ends_with('.1')) %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "(.+?)\\.1$", "\\1"))

exons.rpkm <- exons %>%
  select(-(2:6)) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  rpkm(prior.count = 0.1, gene.length=exons$Length)

ExpressedGenes <- read_tsv("ExpressionAnalysis/polyA/ExpressedGeneList.txt", col_names = c("chr", "start", "stop", "gene", "score", "strand")) %>% pull(gene)

MedExonCov <- data.frame(MedExonCoverage =
             (exons.rpkm %>% apply(1, median))) %>%
  rownames_to_column("ExonID") %>%
  separate(ExonID, into=c("Geneid", "ExonNum"), sep="_") %>%
  filter(Geneid %in% ExpressedGenes) %>%
  mutate(ExonNum = as.numeric(ExonNum))

dim(MedExonCov %>% head())




