library(dplyr)
library(tidyverse)
library(data.table)
library(susieR)
set.seed(0)



get.gene.coords <- function(perm, gene){
    chrom <- perm %>% filter(phe_id==gene) %>% pull(phe_chr) %>% as.character()
    start <- perm %>% filter(phe_id==gene) %>% pull(phe_from) %>% as.integer()
    end <- perm %>% filter(phe_id==gene) %>% pull(phe_to) %>% as.integer()
    
    start <- max(0, start - 200000) %>% as.character()
    end <- (end + 200000) %>% as.character()
    
    loc <- paste(start, end, sep='-')
    coords <- paste(chrom, loc, sep=':')
    return(coords)
}

get.nominal.slice <- function(nominal.file, gene, coords){

    nominal.slice <- tabix(coords, nominal.file, check.chr=FALSE, verbose=FALSE) %>%
        filter(phe_id == gene)
    return (nominal.slice)
}

get.genotype <- function(genotype, coords){
    coords_gsub <- gsub('*chr', '', coords)
    head_cmd <- paste("zcat", genotype, "| head")
    header <- fread(cmd = str_glue(head_cmd),
                      data.table = F, header = T) %>% colnames()
    samples <- header[6:length(header)]
    tabix_cmd <- paste("tabix", genotype, coords_gsub)
    GT <- fread(cmd = str_glue(tabix_cmd),
                  data.table = F, header = F)
    colnames(GT) <- header
    GT <- GT[rowSums(GT[,samples]) > 0,]
    return(GT)
    
}

filter.GT <- function(nominal.gene, GT){
    GT <- GT %>% filter(ID %in% nominal.gene$var_id)
    nominal.gene <- nominal.gene %>% filter(var_id %in% GT$ID)
    return (list(GT=GT, nominal.gene=nominal.gene))
}

args = commandArgs(trailingOnly=TRUE)
nominal = args[1]
perm = args[2]
genotype = args[3]
Subset = args[4]
output = args[5]

perm <- read.table(perm, sep=' ', header=TRUE)


