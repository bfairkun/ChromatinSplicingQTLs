#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PreparePhenotypeTableFromFeatureCounts_SubsetGeneList
# @created     : Tuesday Jun 01, 2021 16:28:47 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/chRNA.Expression/Counts.txt ExpressionAnalysis/polyA/ExpressedGeneList.txt scratch/chRNA.Expression.all.bed.gz scratch/chRNA.Expression.OnlyFirstReps.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

GeneCounts_f_in <- args[1]
Genes_bed_f_in <- args[2]
f_out <- args[3]

### helper "ColumnRenamer" functions to rename the filename column names from featureCounts to the sampleIDs as used in the vcf for qtl calling
rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- "rename_STAR_alignment_samples"

gene.list <- read_tsv(Genes_bed_f_in, col_names=c("Chr", "Start", "End", "Geneid", "score", "Strand"))

dat.genes <- read_tsv(GeneCounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) %>%
    filter(Geneid %in% gene.list$Geneid)

dat.cpm <- dat.genes %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)

#Standardize across individuals (rows),
dat.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)


Out <- gene.list %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)

# Write all samples out
write_tsv(Out, f_out)


