#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : LeafcutterToPSIBed
# @created     : Thursday Jan 20, 2022 09:48:54 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.Expression.Splicing.gz QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.PSI.bed", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_in <- args[1]
f_out <- args[2]

library(tidyverse)
library(data.table)

counts <- fread(f_in, sep=' ')

counts %>%
    # head(1000) %>% as_tibble() %>%
    gather("Sample", "Count", -chrom) %>%
    mutate(Count = as.numeric(str_replace(Count, "(.+?)/.+", "\\1"))) %>%
    mutate(Cluster = str_replace(chrom, ".+:(.+)$", "\\1")) %>%
    group_by(Sample, Cluster) %>%
    mutate(ClusterSum = sum(Count)) %>%
    mutate(PSI = signif(Count/ClusterSum, digits=5)) %>%
    ungroup() %>%
    select(-Count, -Cluster, -ClusterSum) %>%
    spread("Sample", "PSI") %>%
    separate(chrom, into=c("#Chr", "start", "end", "cluster"), sep=":", remove=F) %>%
    mutate(pid = str_remove(chrom, "^chr")) %>%
    mutate(gid = paste(`#Chr`, cluster, sep="_")) %>%
    mutate(strand = str_extract(cluster, "[+-]")) %>%
    select(c("#Chr", "start", "end", "pid", "gid", "strand"), everything()) %>%
    select(-cluster, -chrom) %>%
    write_tsv(f_out)
