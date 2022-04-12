#!/usr/bin/env Rscript

if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI. SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


PrefixOut <- args[1]
leafcutter_count_numers <- args[2]

library(tidyverse)

## Move this to different script
# Make PSI tables
Count.Table.mat <- read.table(leafcutter_count_numers, sep = ' ', nrows = Inf) %>%
  as.matrix()

ClusterMax.mat <- Count.Table.mat %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
    group_by(cluster) %>%
    mutate(across(where(is.numeric), max)) %>%
    ungroup() %>%
    select(junc, everything(), -cluster) %>%
    column_to_rownames("junc") %>%
    as.matrix()


PSI.df <- (Count.Table.mat / as.numeric(ClusterMax.mat) * 100) %>%
    signif() %>%
    as.data.frame()

CountTablePhenotypes <- colnames(Count.Table.mat)[-1] %>%
    str_replace("^(.+?)_.+?_.+$", "\\1") %>% unique()

for (p in CountTablePhenotypes){
    print(p)
    PSI.df %>%
        rownames_to_column("junc") %>%
        select(junc, starts_with(p) & ends_with("_1")) %>%
        rename_with(~ str_replace(.x, "^.+?_(.+?)_.+$", "\\1"), starts_with(p)) %>%
        separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
        mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
        mutate(strand = str_extract(cluster, "[+-]")) %>%
        select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
        arrange(`#Chrom`, start, end) %>%
        write_tsv(paste0(PrefixOut, p, ".bed"))
}
