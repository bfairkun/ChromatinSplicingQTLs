#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Clasify_sQTLs
# @created     : Thursday Aug 10, 2023 10:39:26 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "pi1/PairwisePi1Traits.P.all.txt.gz scratch/test.sQTLs.tsv.gz scratch/tset.sQTLs.full.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

FileIn <- args[1]
FileOut <- args[2]
FileOut_Full <- args[3]

library(data.table)
library(tidyverse)

PC1.PossibleValues <- c("polyA.Splicing")
PC2.PossibleValues <-c("Expression.Splicing",  "H3K4ME3", "H3K36ME3", "H3K27AC", "H3K4ME1")


dat <- fread(FileIn) %>%
  filter((PC1 %in% PC1.PossibleValues) & (PC2 %in% PC2.PossibleValues))

Intron.Annotations <- read_tsv("../data/IntronAnnotationsFromYang.Updated.tsv.gz") %>%
  mutate(IntronName = paste(chrom, start, end, strand, sep=":"))

dat.sQTLs.eQTLs.ForScatter <- dat %>%
  mutate(IntronName = str_replace(P1, "^(.+?:)clu_.+?([+-])$", "chr\\1\\2")) %>%
  mutate(ClusterName =  str_replace(P1, "^(.+?:).+?(clu_.+?[+-])$", "chr\\1\\2")) %>%
  left_join(Intron.Annotations)

PC1.filter = c("polyA.Splicing")
PC2.filter = c( "Expression.Splicing")
PC2.SignificanceFilter <- c("H3K4ME3", "H3K27AC", "H3K36ME3")

sQTLs.ByType.Full <- dat.sQTLs.eQTLs.ForScatter %>%
  filter(PC1 %in% PC1.filter) %>%
  group_by(P1) %>%
  filter(!any((PC2 %in% PC2.SignificanceFilter) & (trait.x.p.in.y < 0.01))) %>%
  ungroup() %>%
  filter(PC2 %in% PC2.filter) %>%
  mutate(SuperAnnotation.simplified = recode(SuperAnnotation, "UnannotatedJunc_UnproductiveCodingGene"="Unproductive", "AnnotatedJunc_ProductiveCodingGene"="Productive", "AnnotatedJunc_UnproductiveCodingGene"="Unproductive", "UnannotatedJunc_ProductiveCodingGene"="Productive")) %>%
  group_by(ClusterName) %>%
  mutate(sQTL.type = case_when(
    any(SuperAnnotation.simplified == "Unproductive")  ~ "u-sQTL",
    all(SuperAnnotation.simplified == "Productive") ~ "p-sQTL",
    TRUE ~ "Other"
    )) %>%
  ungroup() %>%
  filter(!sQTL.type=="Other")

sQTLs.ByType <- sQTLs.ByType.Full %>%
    group_by(ClusterName) %>%
    slice(which.min(p_permutation.x)) %>%
    ungroup()

count(sQTLs.ByType, sQTL.type)

sQTLs.ByType %>%
    write_tsv(FileOut)

sQTLs.ByType.Full %>%
    write_tsv(FileOut_Full)
