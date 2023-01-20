#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PreparePhenotypeTablesFromFeatureCounts
# @created     : Monday May 10, 2021 10:11:08 CDT
#
# @description : from featureCounts output, prepare standardized and
# qqnormalized phenotype table and also table of covariates (genotype PCs,
# expression PCs and sex) for QTLtools
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/H3K27AC/Counts.txt 1000000 QTLs/QTLTools/H3K27AC/AllReps.qqnorm.bed.gz QTLs/QTLTools/H3K27AC/OnlyFirstReps.qqnorm.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

featureCounts_FileIn <- args[1]
MaxFeatures <- as.numeric(args[2])
f_out_all <- args[3]
f_out_only_first_rep <- args[4]
phenotype <- args[5]


### helper "ColumnRenamer" functions to rename the filename column names from featureCounts to the sampleIDs as used in the vcf for qtl calling

if (phenotype == 'DNaseISensitivity') {
    rename_hisat2_alignment_samples <- function(MyString){
    return(
           str_replace(MyString,"Alignments/Hisat2_Align/.+?/(.+?)\\.merged\\.wasp_filterd\\.markdup\\.sorted\\.bam", "\\1")
        )
}
    } else {
rename_hisat2_alignment_samples <- function(MyString){
    return(
           str_replace(MyString,"Alignments/Hisat2_Align/.+?/(.+?)\\.(\\d+)\\.wasp_filterd\\.markdup\\.sorted\\.bam", "\\1.\\2")
    )
}
}

ColumnRenamerFunction <- "rename_hisat2_alignment_samples"
autosomes <- paste0("chr", 1:22)

dat <- read_tsv(featureCounts_FileIn, comment = "#") %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    # select(1:6, matches("\\.1$")) %>%
    # rename_with(~str_remove(., '\\.1$')) %>%
    filter(Chr %in% autosomes)


dat.cpm <- dat %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)


#Filter for top N autosomal genes based on median expression
MedCpm <- sort(apply(dat.cpm, 1, median), decreasing=T)
GenesToInclude <- MedCpm %>% head(MaxFeatures) %>% names()
dat.cpm.filtered <- dat.cpm[GenesToInclude,] %>%  as.matrix()
print(paste("The ", length(GenesToInclude), " included genes are a cutoff of", 2**(MedCpm %>% head(MaxFeatures) %>% tail(1)), "cpm"))


#Standardize across individuals (rows),
dat.standardized <- dat.cpm.filtered %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
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
    arrange(`#Chr`, start)

# Write all samples out
write_tsv(Out, f_out_all)

if (phenotype == 'DNaseISensitivity'){
    Out %>%
    select(1:6, 8:77) %>%
    write_tsv(f_out_only_first_rep)
    } else {
    Out %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) %>%
    write_tsv(f_out_only_first_rep)
}