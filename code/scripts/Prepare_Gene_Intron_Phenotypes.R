library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)


args <- commandArgs(trailingOnly=TRUE)

IntronCounts_in <- args[1]
output <- args[2]


X <- read_tsv(IntronCounts_in, n_max=Inf) 

dat.matrix <- X %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c('Chr', 'start', 'end', 'pid', 'gid', 'strand', 'length'))

dat.cpm <- dat.matrix %>% 
    cpm(log=T, prior.count=0.1)

dat.rpkm <- dat.matrix %>% 
    rpkm(gene.length=X$length, prior.count=0.1)

#bed <- X %>% mutate(Score=".") %>% select(Chr, start, end, pid, gid, strand)

#RPKM.Out <- bed %>%
#    select(gid, Chr, Start, End, Strand) %>%
#    inner_join(
#               (dat.rpkm %>% as.data.frame() %>% rownames_to_column("gid")),
#               by = "Geneid") %>%
#    mutate(across(where(is.numeric), round, 5)) %>%
#    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=gid, gid=gid, strand=Strand, everything()) %>%
#    arrange(`#Chr`, start)



Genes_bed_f_in <- "ExpressionAnalysis/polyA/ExpressedGeneList.txt" 
gene.list <- read_tsv(Genes_bed_f_in, col_names=c("Chr", "Start", "End", "Geneid", "score", "Strand"))

dat.select = dat.cpm[rownames(dat.cpm) %in% gene.list$Geneid, ]
dat.standardized <- dat.select %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)

qqnorm.Out <- gene.list %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)

rpkm.Out <- gene.list %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.rpkm %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)

write_tsv(qqnorm.Out, "QTLs/QTLTools/GeneIntronCounts/OnlyFirstReps.qqnorm.bed.gz")
write_tsv(rpkm.Out, "QTLs/QTLTools/GeneIntronCounts/OnlyFirstReps.RPKM.bed.gz")
