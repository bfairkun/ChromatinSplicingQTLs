#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : Collapse_Juncfiles
# @created     : Tuesday Sep 06, 2022 15:01:58 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SmallMolecule/FullSpliceSiteAnnotations/ALL_SAMPLES.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_3160_LCL_chRNA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_3160_LCL_chRNA_2.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_3160_LCL_chRNA_3.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_100_LCL_chRNA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_100_LCL_chRNA_2.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_100_LCL_chRNA_3.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_chRNA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_chRNA_2.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_chRNA_3.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_polyA_2.junc SmallMolecule/leafcutter/juncfiles/autosomes/DMSO_NA_LCL_polyA_3.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_10000_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_3160_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_1000_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_316_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_100_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_31.6_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_10_LCL_polyA_1.junc SmallMolecule/leafcutter/juncfiles/autosomes/Risdiplam_3.16_LCL_polyA_1.junc", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

fout <- args[1]
f_in_list <- args[-1]

library(tidyverse)
library(data.table)

dat <- lapply(f_in_list, fread) %>%
    bind_rows()

dat %>%
    separate(V11, into=c("bSize1", "bSize2"), sep=",", convert=T, remove=F) %>%
    mutate(JuncStart = V2 + bSize1) %>%
    mutate(JuncEnd = V3 - bSize2) %>%
    group_by(V1, JuncStart, JuncEnd, V6) %>%
    mutate(SumCount = sum(V5)) %>%
    ungroup() %>%
    distinct(V1, JuncStart, JuncEnd, .keep_all=T) %>%
    dplyr::select(-c("bSize1", "bSize2", "JuncStart", "JuncEnd", "V5")) %>%
    dplyr::select(V1:V4, SumCount, everything()) %>%
    arrange(V1, V2, V3) %>%
    write_tsv(fout, col_names=F)

