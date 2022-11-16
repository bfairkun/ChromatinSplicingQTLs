#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : TempPermuteAndPCA
# @created     : Tuesday May 11, 2021 13:50:56 CDT
#
# @description : Permute QTLtools phenotype table, and calculate number of PCs
# that explain more variance than unpermuted table. Input is a bed.gz file formatted as described in the QTLtools manual.
# Output is a space delimited text file in the covariate file format described by the QTLtools manual
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QTLs/QTLTools/chRNA.IR/OnlyFirstReps.sorted.qqnorm.bed.gz QTLs/QTLTools/chRNA.IR/OnlyFirstReps.sorted.qqnorm.pca", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)

FileIn <- args[1]
FileOut <- args[2]
# FileIn <- "scratch/GSE75220_4su_30_qqnorm.txt.gz"

mat <- read_tsv(FileIn) %>%
    select(-c(1,2,3,5,6)) %>%
    column_to_rownames("pid") %>% as.matrix()

mat.permuted <- matrix(nrow=nrow(mat), ncol=ncol(mat))
for (i in 1:nrow(mat)){
    mat.permuted[i,] <- sample(mat[i,], size=ncol(mat))
}

pca.results <- prcomp(t(mat))

pca <- pca.results %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>%
    rownames_to_column("PC")
pca.permuted <- prcomp(t(mat.permuted)) %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>% rownames_to_column("PC")

merged <- bind_rows(list(pca=pca, pca.permuted=pca.permuted), .id="mat") %>%
    mutate(PC=as.numeric(str_replace(PC, "PC", "")))

#GetNumPCs
NumPCs <- merged %>%
    select(PC, Prop=`Proportion of Variance`, mat) %>%
    spread(key="mat", value="Prop") %>%
    filter(pca > pca.permuted) %>% pull(PC) %>% max()

print("Num PCs for which variance explained is more than in permuted data")
print(NumPCs)

pca.results$x[,1:NumPCs] %>% t() %>%
    round(5) %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    write_delim(FileOut)

# merged %>%
#     ggplot(aes(x=PC, y=`Proportion of Variance`, color=mat)) +
#     geom_line() +
#     theme_bw()
