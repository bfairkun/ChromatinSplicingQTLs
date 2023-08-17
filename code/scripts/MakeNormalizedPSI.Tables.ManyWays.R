#!/usr/bin/env Rscript

if(interactive()){
    args <- scan(text=
                 "scratch/SplicingTablesManyWays. SplicingAnalysis/CombinedJuncTables/All.tsv.gz ../data/IntronAnnotationsFromYang.Updated.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


PrefixOut <- args[1]
junc_counts_all <- args[2]
IntAnnotations <- args[3]

library(tidyverse)
library(data.table)

dat <- fread(junc_counts_all)

IntronAnnotations <- read_tsv(IntAnnotations) %>%
  dplyr::rename(stop=end)

dat.wide <- dat %>%
  split(.$Dataset) %>%
  lapply(pivot_wider, names_from=c("Dataset", "IndID", "RepNumber"), values_from="Count", values_fill=0) %>%
  map(~mutate(.x, stop=stop+1)) %>%
  map(~inner_join(.x, IntronAnnotations))

ClusteredIntrons <- fread("../code/SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz", select=1, col.names=c("Intron")) %>%
  separate(Intron, into=c("chrom", "start", "stop", "cluster"), convert=T, remove=F, sep=":") %>%
  mutate(strand = str_extract(cluster, "[+-]$")) %>%
  mutate(cluster = str_remove(cluster, "_[+-]$"))

IntronNames <- IntronAnnotations %>%
  dplyr::select(chrom:strand) %>%
  left_join(ClusteredIntrons) %>%
  mutate(NewChrom = str_remove(chrom, "^chr")) %>%
  mutate(name = case_when(
    is.na(cluster) ~ str_glue("{NewChrom}:{start}:{stop}:clu_NA_{strand}"),
    TRUE ~ str_glue("{NewChrom}:{start}:{stop}:{cluster}_{strand}")
  )) %>%
  dplyr::select(chrom, start, stop, strand, name, cluster)



test.dat <- dat.wide$chRNA.Expression.Splicing %>% head(10000)

# functions to transform a wide-format counts df

TransformTo_RatioFromMedianProductive <- function(df){
  dat.medians <- df %>%
    filter(SuperAnnotation %in% c("AnnotatedJunc_ProductiveCodingGene")) %>%
    group_by(gene) %>%
    mutate(across(contains(".junc"), median)) %>%
    ungroup() %>%
    dplyr::select(gene, everything())
  dat.medians.mat <- dat.medians %>%
    distinct(gene, .keep_all=T) %>%
    dplyr::select(gene, contains(".junc")) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  dat.AlmostMat <- df %>%
    filter(gene %in% rownames(dat.medians.mat)) %>%
    dplyr::select(chrom, start, stop, strand, gene, contains(".junc")) %>%
    unite(Intron, chrom, start, stop, strand)
  dat.AlmostMat.genes <- dat.AlmostMat$gene
  dat.Mat <- dat.AlmostMat %>%
    dplyr::select(-gene) %>%
    column_to_rownames("Intron") %>%
    as.matrix()
  df.out <- (dat.Mat / dat.medians.mat[dat.AlmostMat.genes,]) %>%
    as.data.frame() %>%
    rownames_to_column("Intron") %>%
    separate(Intron, into=c("chrom", "start", "stop", "strand"), sep='_', convert=T) %>%
    inner_join(IntronNames) %>%
    mutate(score = ".") %>%
    dplyr::select(chrom, start, stop, name, score, strand, contains(".junc"))
  df.out[sapply(df.out, is.infinite)] <- NaN
  
    return(df.out)
}


TransformTo_RatioFromMaxInCluster <- function(df){
    Count.Table.mat <- df %>%
    dplyr::select(chrom:strand, contains(".junc")) %>%
    left_join(IntronNames) %>%
    filter(!is.na(cluster)) %>%
    unite(Intron, chrom, start, stop, strand, cluster, sep=":") %>%
    dplyr::select(Intron, contains(".junc")) %>%
    column_to_rownames("Intron") %>%
    as.matrix()
  ClusterMax.mat <- Count.Table.mat %>%
      as.data.frame() %>%
      rownames_to_column("junc") %>%
      mutate(cluster=str_replace(junc, "^.+?:.+?:.+?:.+?:(.+)$", "\\1")) %>%
      group_by(cluster) %>%
      mutate(across(where(is.numeric), max)) %>%
      ungroup() %>%
      select(junc, everything(), -cluster) %>%
      column_to_rownames("junc") %>%
      as.matrix()
  df.out <- ((Count.Table.mat / ClusterMax.mat)) %>%
    as.data.frame() %>%
    rownames_to_column("Intron") %>%
    separate(Intron, into=c("chrom", "start", "stop", "strand", "cluster"), sep=":", convert=T) %>%
    dplyr::select(-cluster) %>%
    inner_join(IntronNames) %>%
    mutate(score = ".") %>%
    dplyr::select(chrom, start, stop, name, score, strand, contains(".junc"))
  return(df.out)
}

TransformTo_RatioFromClusterTotal <- function(df){
    Count.Table.mat <- df %>%
    dplyr::select(chrom:strand, contains(".junc")) %>%
    left_join(IntronNames) %>%
    filter(!is.na(cluster)) %>%
    unite(Intron, chrom, start, stop, strand, cluster, sep=":") %>%
    dplyr::select(Intron, contains(".junc")) %>%
    column_to_rownames("Intron") %>%
    as.matrix()
  ClusterMax.mat <- Count.Table.mat %>%
      as.data.frame() %>%
      rownames_to_column("junc") %>%
      mutate(cluster=str_replace(junc, "^.+?:.+?:.+?:.+?:(.+)$", "\\1")) %>%
      group_by(cluster) %>%
      mutate(across(where(is.numeric), sum)) %>%
      ungroup() %>%
      select(junc, everything(), -cluster) %>%
      column_to_rownames("junc") %>%
      as.matrix()
  df.out <- ((Count.Table.mat / ClusterMax.mat)) %>%
    as.data.frame() %>%
    rownames_to_column("Intron") %>%
    separate(Intron, into=c("chrom", "start", "stop", "strand", "cluster"), sep=":", convert=T) %>%
    dplyr::select(-cluster) %>%
    inner_join(IntronNames) %>%
    mutate(score = ".") %>%
    dplyr::select(chrom, start, stop, name, score, strand, contains(".junc"))
  return(df.out)
}


rename_func <- function(df){
    rename_with(df, ~str_replace(.x, "^.+?_(.+?)_1.junc", "\\1"),.cols=contains("_1.junc")) %>%
        return()
}

cutoff_columns <- function(df, cutoff=0.5){
    mutate_at(df, ~pmin(.x, cutoff), .vars=vars(contains(".junc"))) %>%
        return()
}

MultiplyColumns <- function(df, n=100){
    mutate_at(df, ~.x*n, .vars=vars(contains(".junc"))) %>%
        return()
}


dat.Out <- dat.wide  %>%
    # lapply(head, 10000) %>%
    lapply(TransformTo_RatioFromClusterTotal) %>%
    map(~dplyr::select(.x, "#Chrom" = chrom, start, stop, name, score, strand, contains("_1.junc"))) %>%
    map(MultiplyColumns) %>%
    map(~mutate_at(.x, .funs=signif, .vars=vars(contains(".junc")))) %>%
    map(rename_func) %>%
    lapply(arrange, `#Chrom`, start, stop) %>%
    mapply(write_tsv, ., file=paste0(PrefixOut, "CountsOverClusterSum_", names(.), '.bed'))

dat.Out <- dat.wide  %>%
    # lapply(head, 10000) %>%
    lapply(TransformTo_RatioFromMaxInCluster) %>%
    map(~dplyr::select(.x, "#Chrom" = chrom, start, stop, name, score, strand, contains("_1.junc"))) %>%
    map(MultiplyColumns) %>%
    map(~mutate_at(.x, .funs=signif, .vars=vars(contains(".junc")))) %>%
    map(rename_func) %>%
    lapply(arrange, `#Chrom`, start, stop) %>%
    mapply(write_tsv, ., file=paste0(PrefixOut, "CountsOverClusterMax_", names(.), '.bed'))

dat.Out <- dat.wide  %>%
    # lapply(head, 10000) %>%
    lapply(TransformTo_RatioFromMedianProductive) %>%
    map(~dplyr::select(.x, "#Chrom" = chrom, start, stop, name, score, strand, contains("_1.junc"))) %>%
    map(MultiplyColumns) %>%
    map(~mutate_at(.x, .funs=signif, .vars=vars(contains(".junc")))) %>%
    map(rename_func) %>%
    lapply(arrange, `#Chrom`, start, stop) %>%
    mapply(write_tsv, ., file=paste0(PrefixOut, "CountsOverMedianProductive_", names(.), '.bed'))

dat.Out <- dat.Out  %>%
    map(cutoff_columns, 100) %>%
    mapply(write_tsv, ., file=paste0(PrefixOut, "CountsOverMedianProductiveCoercedRange_", names(.), '.bed'))
