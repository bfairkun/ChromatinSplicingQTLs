#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PlotJuncsByCategoryForYang
# @created     : Wednesday Jan 18, 2023 12:46:49 CST
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "/project2/yangili1/yangili/chRNA/annotation_leaf.txt.gz scratch/P1", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(viridis)

f_in <- args[1]
P1_out_prefix <- args[2]

print(args)

annotations <- read_tsv(f_in, col_names = c("chrom", "start", "stop", "strand", "OriginalAnnotation", "NewAnnotation"))

JuncSumsAcrossDataTypes <- read_tsv("/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/SplicingAnalysis/regtools_annotate_combined/Comprehensive.YRI.Sample.Sum.Counts.tsv.gz")

ToPlot <- JuncSumsAcrossDataTypes %>%
    separate(junc, into=c("chrom", "start", "stop", "strand"), sep=":", convert=T) %>%
    mutate(strand = str_extract(strand, "[+-]")) %>%
    gather(key="Dataset", value="Counts", -c("chrom", "start", "stop", "strand")) %>%
    inner_join(annotations) %>%
    # group_by(Dataset, NewAnnotation) %>%
    group_by(Dataset, NewAnnotation, OriginalAnnotation) %>%
    summarise(Sum = sum(Counts)) %>%
    ungroup() %>%
    group_by(Dataset) %>%
    mutate(Percent = Sum/sum(Sum)) %>%
    ungroup()

ToPlot2 <- JuncSumsAcrossDataTypes %>%
  separate(junc, into=c("chrom", "start", "stop", "strand"), sep=":", convert=T) %>%
  mutate(strand = str_extract(strand, "[+-]")) %>%
  gather(key="Dataset", value="Counts", -c("chrom", "start", "stop", "strand")) %>%
  filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
  inner_join(annotations) %>%
  group_by(Dataset) %>%
  mutate(Sum = sum(Counts)) %>%
  ungroup() %>%
  mutate(Counts = Counts + 0.1) %>%
  mutate(Percent = Counts/Sum * 100) %>%
  dplyr::select(chrom:strand, Dataset, OriginalAnnotation, NewAnnotation, Percent) %>%
  add_count(chrom, start, stop, strand, OriginalAnnotation, NewAnnotation, Dataset) %>%
  filter(n==1) %>%
  pivot_wider(names_from = "Dataset", values_from = "Percent")

OriginalAnnotationSet <- unique(ToPlot$OriginalAnnotation)



for (annotation in OriginalAnnotationSet){
  print(annotation)
  P <- ToPlot %>%
    filter(OriginalAnnotation == annotation) %>%
    mutate(Dataset = recode(Dataset, !!!c("Expression.Splicing"="polyA RNA", "chRNA.Expression.Splicing"="chRNA", "MetabolicLabelled.30min"="30min 4sU RNA", "MetabolicLabelled.60min"="60min 4sU RNA"))) %>%
    mutate(Dataset = factor(Dataset, levels=c("chRNA", "30min 4sU RNA", "60min 4sU RNA", "polyA RNA"))) %>%
    ggplot(aes(x=Dataset, y=Percent)) +
    geom_col() +
    facet_wrap(~NewAnnotation, scales="free", labeller = label_wrap_gen(width=14)) +
    scale_colour_brewer(type="qual", palette="Dark2") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
    labs(y=str_wrap("Percent of splice junction reads", 20), x=NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(strip.text.x = element_text(size = 12))
  ggsave(paste0(P1_out, "_", annotation, ".pdf"), P, height=20, width=30)

  P2 <- ToPlot2 %>%
    filter(OriginalAnnotation == annotation) %>%
    ggplot(aes(x=chRNA.Expression.Splicing, y=Expression.Splicing)) +
    geom_abline(color='red', slope=1, intercept=0) +
    geom_point(alpha=0.2) +
    scale_fill_continuous(trans="log10") +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10") +
    facet_wrap(~NewAnnotation, labeller = label_wrap_gen(width=14)) +
    theme_bw() +
    labs(x="chRNA, %juncs", y="polyA, %juncs")
  ggsave(paste0(P1_out, "_Scatter_", annotation, ".pdf"), P2, height=20, width=30)

}


