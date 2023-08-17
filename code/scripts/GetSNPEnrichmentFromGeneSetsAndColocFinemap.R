#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : GetSNPEnrichmentFromGeneSetsAndColocFinemap
# @created     : Monday Jul 10, 2023 13:42:17 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTL_SNP_Enrichment/ContrastingGeneSets/MolColocNonRedundantFullSplicing.tsv.gz QTL_SNP_Enrichment/FinemapIntersections/MolColocTestOnlyExpressionFullGeauvadis.bed.gz hyprcoloc/Results/ForColoc/MolColocTestOnlyExpressionFullGeauvadis/tidy_results_OnlyColocalized.txt.gz scratch/test.Enrichment", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

# "hyprcoloc/Results/ForColoc/MolColocTestOnlyExpression//"
# "hyprcoloc/Results/ForColoc/MolColocTestOnlyExpressionFullGeauvadis/"

GeneSetsIn <- args[1]
FinemapIn <- args[2]
FinemapIn.coloc <- args[3]
PrefixOut <- args[4]

library(tidyverse)
library(data.table)



GeneSets <- read_tsv(GeneSetsIn)

print("Number genes in sets based on gene sets file")
count(GeneSets, Class)

dat.coloc <- read_tsv(FinemapIn.coloc) %>%
  separate(phenotype_full, into=c("PC", "P"), sep=";") %>%
  filter(str_detect(PC, "^Expression.Splicing")) %>%
  inner_join(GeneSets, by=c("P"="genes"))



dat.bed <- fread(FinemapIn, col.names=c("SNP_chrom", "SNP_start", "SNP_end", "SNP", "PosteriorPr", "Annotation_chrom", "Annotation_start", "Annotation_stop", "Annotation", "Overlap")) %>%
  separate(SNP, into=c("SNP", "HyprcolocCluster", "GeneLocus"), sep="_", convert=T)

dat.bed <- bind_rows(
  dat.bed,
  dat.bed %>%
    filter(
      str_detect(Annotation, "^Splice.+_[01]$")) %>%
    mutate(Annotation = case_when(
      str_detect(Annotation, "^Splice.+_1$") ~ "Annotated splice site",
      str_detect(Annotation, "^Splice.+_0$") ~ "Unannotated splice site",
      TRUE ~ Annotation
    )))

print("Number genes in sets after filtering for those colocalized in coloc run for finemapping")
dat.bed %>%
  inner_join(
    dat.coloc %>%
      sample_frac(replace=F),
    by=c("GeneLocus"="Locus", "HyprcolocCluster"="iteration")) %>%
  distinct(P, Class) %>%
  count(Class)


set.seed(0)
N <- 100
dat = list()
for (i in 1:N){
  print(paste("Running", i))
  dat[[i]] <- dat.bed %>%
    inner_join(
      dat.coloc %>%
        sample_frac(replace=T),
      by=c("GeneLocus"="Locus", "HyprcolocCluster"="iteration")) %>%
    group_by(Class, Annotation) %>%
    summarise(PipInAnnotation = sum(PosteriorPr)) %>%
    ungroup() %>%
    pivot_wider(names_from="Class", values_from="PipInAnnotation") %>%
    mutate(Enrichment = log2(Other_eGenes/hQTL_eGenes)) %>%
    mutate(BootstrapRep = i)
}

dat.to.plot <- bind_rows(dat)

dat.to.plot %>%
  write_tsv(paste0(PrefixOut, "dat.tsv.gz"))

dat.to.plot %>%
  group_by(Annotation) %>%
  summarise(
    mean = mean(Enrichment, na.rm=T),
    Lower = quantile(Enrichment, probs = 0.05, na.rm=T),
    Upper = quantile(Enrichment, probs = 0.95, na.rm=T)) %>%
ggplot(aes(x=Annotation, y=mean)) +
  geom_col() +
  geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
  coord_flip() +
  theme_bw()

ggsave(paste0(PrefixOut, "P.pdf"), height=8, width=7)





