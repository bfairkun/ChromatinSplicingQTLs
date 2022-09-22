library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)


args <- commandArgs(trailingOnly=TRUE)

IntronCounts_in <- args[1]
output <- args[2]


X <- read_tsv(IntronCounts_in, n_max=Inf) 

dat.matrix <- X %>%
    column_to_rownames("gene_ID") %>%
    select(everything(), -c('Chr', 'start', 'end', 'gene_ID', 'length'))

dat.cpm <- dat.matrix %>% 
    cpm(log=T, prior.count=0.1)

dat.rpkm <- dat.matrix %>% 
    rpkm(gene.length=X$length, prior.count=0.1)

bed <- X %>% mutate(Score=".") %>% select(Chr, start, end, gene_ID, Score, Strand)