#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

args <- commandArgs(trailingOnly=TRUE)

featureCounts_FileIn <- args[1]
PhenotypesBedOut_OnlyFirstReps <- args[2]

print('hola')
table_genes <- read_tsv(featureCounts_FileIn) #%>% select(Chr, Start, End, Strand, Geneid)

print('como')

genes_bed <- table_genes %>% select(Chr, Start, End, Strand, Geneid)

print('estan')

dat <- table_genes %>% select(-c("Strand", "Chr", "Start", "End")) %>% inner_join(genes_bed, ., by="Geneid")

print('todos')

#genes_bed <- table_genes %>% select(Chr, Start, End, Strand, Geneid)

#print('como')

#dat <- read_tsv(featureCounts_FileIn, comment='#') %>%
#    select(-c("Strand", "Chr", "Start", "End")) %>%
#    inner_join(genes_bed, ., by="Geneid") 

print('aqui')
dat.cpm <- dat %>%
    filter(Chr %in% paste0("chr", 1:22)) %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)

print('en')
#Filter for top N autosomal genes based on median expression
MedCpm <- sort(apply(dat.cpm, 1, median), decreasing=T)
NGenesToInclude <- min(14000, dim(dat.cpm)[1])
print(NGenesToInclude)
GenesToInclude <- MedCpm[1:NGenesToInclude] %>% names()
dat.cpm.filtered <- dat.cpm[GenesToInclude,] %>% as.matrix()
print(paste("The ", NGenesToInclude, " included genes are a cutoff of", 2**MedCpm[NGenesToInclude], "cpm"))

print('este')
#Standardize across individuals (rows),
dat.standardized <- dat.cpm.filtered %>% t() %>% scale() %>% t() %>% na.omit()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)

print('hermoso')

Out <- dat %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>% 
    select(-c('Start')) %>%
    arrange(`#Chr`, start)

print('dia')

write_tsv(Out, PhenotypesBedOut_OnlyFirstReps)

