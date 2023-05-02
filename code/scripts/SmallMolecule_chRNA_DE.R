#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : SmallMolecule_chRNA_DE
# @created     : Wednesday Apr 26, 2023 14:24:53 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SmallMolecule/featureCounts/Counts.txt scratch/DE.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(edgeR)

f_in <- args[1]
f_out <- args[2]

IntronAnnotations <- read_tsv("../data/IntronAnnotationsFromYang.tsv.gz") %>%
    filter(!str_detect(SuperAnnotation, "Noncoding"))

ProteinCodingGenes <- IntronAnnotations %>%
    filter(!is.na(symbol)) %>%
    distinct(gene, symbol)

Counts <- read_tsv(f_in, comment="#") %>%
  rename_at(vars(-c(1:6)), ~str_replace(.x, "SmallMolecule/AlignmentsPass2/(.+?)/Aligned.sortedByCoord.out.bam", "\\1")) %>%
  dplyr::select(Geneid, Length, contains("chRNA")) %>%
  filter(Geneid %in% ProteinCodingGenes$gene) %>%
  inner_join(ProteinCodingGenes, by=c("Geneid"="gene")) %>%
  mutate(Geneid = paste(Geneid, symbol, sep="_"))

Counts %>% colnames()

group <- relevel(as.factor(c(3160, 3160, 3160, 100, 100, 100, "DMSO", "DMSO", "DMSO")), "DMSO")

y <- Counts %>%
    dplyr::select(-Length, -symbol) %>%
    column_to_rownames("Geneid") %>%
    DGEList(group = group)

design <- model.matrix(~group)

CONTRASTS <- colnames(design)[-1] %>%
  makeContrasts(contrasts=., levels = design )

keep <- filterByExpr(y)
filtered <- y[keep, , keep.lib.sizes=FALSE] %>%
  calcNormFactors()

fit <- filtered %>%
  estimateDisp(design) %>%
  glmQLFit(design)

N <- ncol(CONTRASTS)
results <- vector("list", N) %>%
  setNames(colnames(CONTRASTS))

for (i in 1:ncol(CONTRASTS)){
  results[[i]] <- glmQLFTest(fit, contrast = CONTRASTS[,i]) %>%
    topTags(n=Inf) %>% as.data.frame() %>%
    rownames_to_column("Geneid")
}


results.df <- bind_rows(results, .id="risdiplam_conc") %>%
    mutate(risdiplam_conc = as.numeric(str_replace_all(risdiplam_conc, "group", "")))

write_tsv(results.df, f_out)
