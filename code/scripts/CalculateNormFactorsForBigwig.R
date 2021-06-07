#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : CalculateNormFactorsForBigwig
# @created     : Friday May 21, 2021 09:52:57 CDT
#
# @description : Calculate normalization factors (TMM method) for libraries in
# featureCounts output, and summarize in text file.
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/Counts.txt scratch/NormFactors.tsv QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz 4", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

featureCounts_f <- args[1]
f_out <- args[2]
# File containing genes to subset.
optional_genes_to_include_f <- args[3]
# coluns number with Geneids to subset. If not specified, use 1
column_number_with_Geneid <- args[4]

if (is.na(column_number_with_Geneid)){
    column_number_with_Geneid <- 1
}


library(tidyverse)
library(edgeR)

mat <- read_tsv(featureCounts_f, comment="#") %>%
    select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
    column_to_rownames("Geneid") %>%
    DGEList()


if (!is.na(args[3])){
    Geneids <- read_tsv(optional_genes_to_include_f) %>% pull(as.numeric(column_number_with_Geneid))
    head(Geneids)
    mat <- mat[Geneids,]
}

mat <- calcNormFactors(mat)

EffectiveLibSize =mat$samples$norm.factors * mat$samples$lib.size
Out <- data.frame(Sample = rownames(mat$samples), EffectiveLibSize = EffectiveLibSize)

write_tsv(Out, f_out)
