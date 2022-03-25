#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

args <- commandArgs(trailingOnly=TRUE)

featureCounts_FileIn <- args[1]
PhenotypesBedOut_OnlyFirstReps <- args[2]

rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- 'rename_STAR_alignment_samples'

table_genes <- read_tsv(featureCounts_FileIn, comment = "#") #%>% select(Chr, Start, End, Strand, Geneid)


genes_bed <- table_genes %>% select(Chr, Start, End, Strand, Geneid)


# dat <- table_genes %>% select(-c("Strand", "Chr", "Start", "End")) %>% inner_join(genes_bed, ., by="Geneid")

dat <- read_tsv(featureCounts_FileIn, comment = "#") %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(-c("Strand", "Chr", "Start", "End")) %>%
    inner_join(genes_bed, ., by="Geneid") %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))
#     mutate(Chr=paste0("chr", Chr))



#genes_bed <- table_genes %>% select(Chr, Start, End, Strand, Geneid)

#print('como')

#dat <- read_tsv(featureCounts_FileIn, comment='#') %>%
#    select(-c("Strand", "Chr", "Start", "End")) %>%
#    inner_join(genes_bed, ., by="Geneid") 

dat.cpm <- dat %>%
    filter(Chr %in% paste0("chr", 1:22)) %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)

#Filter for top N autosomal genes based on median expression
MedCpm <- sort(apply(dat.cpm, 1, median), decreasing=T)
NGenesToInclude <- min(14000, dim(dat.cpm)[1])
print(NGenesToInclude)
GenesToInclude <- MedCpm[1:NGenesToInclude] %>% names()
dat.cpm.filtered <- dat.cpm[GenesToInclude,] %>% as.matrix()
print(paste("The ", NGenesToInclude, " included genes are a cutoff of", 2**MedCpm[NGenesToInclude], "cpm"))

#Standardize across individuals (rows),
dat.standardized <- dat.cpm.filtered %>% t() %>% scale() %>% t() %>% na.omit()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)


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


write_tsv(Out, PhenotypesBedOut_OnlyFirstReps)

