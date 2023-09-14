#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : HandpickedColocExampleScatterPlots
# @created     : Thursday Sep 07, 2023 20:38:51 CDT
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
library(data.table)

phenotypes <- c("2:241743815:241750151:clu_8613_+", "ENSG00000180902.18")
gwas.accession <- "GCST007799"
gwas_locus.ofinterest <- "chr2_241759225_N_N_GCST007799"

molQTL.dat <- fread(str_glue("hyprcoloc/LociWiseSummaryStatsInput/ForGWASColoc/{gwas.accession}.txt.gz")) %>%
  mutate()
gwas.dat <- fread(str_glue("gwas_summary_stats/StatsForColoc/{gwas.accession}.standardized.txt.gz"))

molQTL.dat.atLocus <- molQTL.dat %>%
  filter(phenotype %in%  phenotypes) %>%
  # pull(phenotype) %>% unique()
  mutate(PhenotypeClass = str_replace(source_file, "^QTLs/QTLTools/(.+?)/NominalPassForGWASColocChunks/.+$", "\\1"))

molQTL.dat %>%
  filter(gwas_locus==gwas_locus.ofinterest)  %>%
  pull(phenotype) %>% unique()

molQTL.dat.atLocus %>%
  filter(gwas_locus==gwas_locus.ofinterest) %>%
  pull(phenotype) %>% unique()

coloc.plots.example.dat.Gene1 <- bind_rows(
  molQTL.dat.atLocus %>%
    filter(gwas_locus==gwas_locus.ofinterest) %>%
    mutate(Signal = -log10(p)) %>%
    dplyr::select(Signal, PhenotypeClass, snp ) %>%
    separate(snp, into=c("chrom", "pos", "ref", "alt"), sep=":", remove=F) %>%
    mutate(pos = as.numeric(pos)),
  gwas.dat %>%
    filter(loci==gwas_locus.ofinterest) %>%
    mutate(Signal = -log10(2*pnorm( abs(beta/SE), lower=F ))) %>%
    mutate(PhenotypeClass = "gwas") %>%
    dplyr::select(Signal, PhenotypeClass, chrom, pos=start, ref=A1, alt=A2) %>%
    mutate(chrom = str_replace(chrom, "^chr(.+$)", "\\1")) %>%
    mutate(pos= as.numeric(pos))
) %>%
  dplyr::select(-ref, -alt) %>%
  pivot_wider(names_from = "PhenotypeClass", values_from = c("Signal","snp"), id_cols=c("pos", "chrom")) %>%
  dplyr::rename(snp=snp_Expression.Splicing) %>%
  dplyr::select(-snp_gwas, -snp_polyA.Splicing) %>%
  mutate(IsTopSNP = chrom=="2" & pos==241759225) %>%
  mutate(gene = "D2HGDH") %>%
  filter(!is.na(snp))


### NUDT14 example

phenotypes <- c("14:105175987:105176534:clu_34934_-", "ENSG00000183828.15")
gwas.accession <- "GCST004611"
gwas_locus.ofinterest <- "chr14_105178084_N_N_GCST004611"

molQTL.dat <- fread(str_glue("hyprcoloc/LociWiseSummaryStatsInput/ForGWASColoc/{gwas.accession}.txt.gz")) %>%
  mutate()
gwas.dat <- fread(str_glue("gwas_summary_stats/StatsForColoc/{gwas.accession}.standardized.txt.gz"))

molQTL.dat.atLocus <- molQTL.dat %>%
  filter(phenotype %in%  phenotypes) %>%
  # pull(phenotype) %>% unique()
  mutate(PhenotypeClass = str_replace(source_file, "^QTLs/QTLTools/(.+?)/NominalPassForGWASColocChunks/.+$", "\\1"))

molQTL.dat %>%
  filter(gwas_locus==gwas_locus.ofinterest)  %>%
  pull(phenotype) %>% unique()

molQTL.dat.atLocus %>%
  filter(gwas_locus==gwas_locus.ofinterest) %>%
  pull(phenotype) %>% unique()

coloc.plots.example.dat.Gene2 <- bind_rows(
  molQTL.dat.atLocus %>%
    filter(gwas_locus==gwas_locus.ofinterest) %>%
    mutate(Signal = -log10(p)) %>%
    dplyr::select(Signal, PhenotypeClass, snp ) %>%
    separate(snp, into=c("chrom", "pos", "ref", "alt"), sep=":", remove=F) %>%
    mutate(pos = as.numeric(pos)),
  gwas.dat %>%
    filter(loci==gwas_locus.ofinterest) %>%
    mutate(Signal = -log10(2*pnorm( abs(beta/SE), lower=F ))) %>%
    mutate(PhenotypeClass = "gwas") %>%
    dplyr::select(Signal, PhenotypeClass, chrom, pos=start, ref=A1, alt=A2) %>%
    mutate(chrom = str_replace(chrom, "^chr(.+$)", "\\1")) %>%
    mutate(pos= as.numeric(pos))
) %>%
  dplyr::select(-ref, -alt) %>%
  distinct(chrom, pos,  PhenotypeClass, .keep_all=T) %>%
  pivot_wider(names_from = "PhenotypeClass", values_from = c("Signal","snp"), id_cols=c("pos", "chrom") ) %>%
  dplyr::rename(snp=snp_Expression.Splicing) %>%
  dplyr::select(-snp_gwas, -snp_polyA.Splicing) %>%
  mutate(IsTopSNP = chrom=="14" & pos==105172594) %>%
  mutate(gene = "NUDT14") %>%
  filter(!is.na(snp))

coloc.plots.example.dat.Gene1 %>%
    filter(Signal_Expression.Splicing>30 & Signal_gwas>15) %>%
    head(1) %>%
    print(width=Inf)

bind_rows(coloc.plots.example.dat.Gene2 %>%
              mutate(TopSNP="14:105172594:T:C"),
          coloc.plots.example.dat.Gene1 %>%
              mutate(TopSNP="2:241759875:C:T"),
          coloc.plots.example.dat.Gene1 %>%
              mutate(TopSNP="2:241759225:G:A")) %>%
    filter(!is.na(Signal_Expression.Splicing)) %>%
    write_tsv("Misc/HandpickedColocPlotsToHighlight.dat.NoLD.tsv.gz")
