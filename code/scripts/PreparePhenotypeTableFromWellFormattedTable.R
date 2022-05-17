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
                 "MiscCountTables/H3K36ME3.bed  QTLs/QTLTools/H3K27AC/OnlyFirstReps.qqnorm.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

featureCounts_FileIn <- args[1]
f_out <- args[2]


### helper "ColumnRenamer" functions to rename the filename column names from featureCounts to the sampleIDs as used in the vcf for qtl calling

dat <- read_tsv(featureCounts_FileIn)


dat.cpm <- dat %>%
    select(-c(1,2,3,5,6)) %>%
    column_to_rownames("pid") %>%
    cpm(log=T, prior.count=0.1)



#Standardize across individuals (rows),
dat.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)


Out <- dat %>%
    select(1:6) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("pid")),
               by = "pid") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`, start, end, pid, gid, strand, everything()) %>%
    arrange(`#Chr`, start)


Out %>%
    write_tsv(f_out)
