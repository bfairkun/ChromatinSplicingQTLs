#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : ReadLengthMapExperiment.ProcessJuncFiles
# @created     : Wednesday Feb 01, 2023 16:00:57 CST
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

out_f_summed_by_dataset <- args[1]
out_f_longtable <- args[2]
out_f_summed_by_annotation <- args[3]

library(tidyverse)
library(data.table)

dat <- Sys.glob("ReadLengthMapExperimentSpliceCounts/juncfiles/8/*/*.1.junc.gz") %>%
    setNames(str_replace(., "^ReadLengthMapExperimentSpliceCounts/juncfiles/8/(.+?)/(.+?).1.junc.gz$", "\\1_\\2")) %>%
    lapply(fread, sep='\t') %>%
    bind_rows(.id="Sample") %>%
    separate(Sample, into=c("Dataset","IndID"), sep="_")

Annotations <- read_tsv("/project2/yangili1/yangili/chRNA/annotation_leaf_JAN28.txt.gz", col_names=c("chrom", "start","end", "strand","OldAnnotation", "NewAnnotation", "gene", "symbol"))

Annotations %>%
    filter(OldAnnotation == "Unannotated") %>%
    count(NewAnnotation) %>% 
    mutate(percent=n/sum(n)*100) %>%
    arrange(n) %>%
    print(n=Inf)

LongTable <- dat %>%
    dplyr::select(Dataset:strand) %>%
    # head(1000000) %>%
    inner_join(Annotations) %>%
    dplyr::select(-symbol)

SummedBySampleAndAnnotation <- LongTable %>%
    group_by(Dataset, IndID, NewAnnotation) %>%
    summarise(Count = sum(score) + 0.00) %>%
    ungroup()

SummedByDataset <- LongTable %>%
    group_by(Dataset, chrom, start, end, strand, NewAnnotation, gene) %>%
    summarise(Count = sum(score) + 0.00) %>%
    ungroup()

write_tsv(LongTable, out_f_longtable)
write_tsv(SummedByDataset, out_f_summed_by_dataset)
write_tsv(SummedBySampleAndAnnotation, out_f_summed_by_annotation)

#DatasetTotals <- SummedByDataset %>%
#    group_by(Dataset) %>%
#    summarise(TotalJuncsInDataset = sum(Count)) %>%
#    ungroup()

#ToPlot <- SummedByDataset %>%
#    inner_join(DatasetTotals) %>%
#    mutate(FractionTotalJuncs = Count / TotalJuncsInDataset)

#P <- ToPlot %>%
#    filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
#    dplyr::select(-TotalJuncsInDataset, -Count) %>%
#    pivot_wider(names_from="Dataset", values_from = "FractionTotalJuncs") %>%
#    add_count(NewAnnotation) %>%
#    filter(n>1000 ) %>%
#    group_by(NewAnnotation) %>%
#    sample_n(1000) %>%
#    ungroup() %>%
#    # head() %>% print(width=Inf)
#    ggplot(aes(x=chRNA.Expression.Splicing, y=Expression.Splicing)) +
#    geom_point(alpha=0.2) +
#    geom_abline(color='red') +
#    scale_x_continuous(trans='log10', lim=c(1E-11,1E-3)) +
#    scale_y_continuous(trans='log10', lim=c(1E-11,1E-3) ) +
#    facet_wrap(~NewAnnotation) +
#    theme_bw() +
#    labs(title="New mappings")

##ggsave("scratch/ReadLengthMapExperiment.ScatterByAnnotation.png", P, height=15, width=15)

#P2 <- ToPlot %>%
#    filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
#    dplyr::select(-TotalJuncsInDataset, -Count) %>%
#    pivot_wider(names_from="Dataset", values_from = "FractionTotalJuncs") %>%
#    mutate(Ratio = log10(Expression.Splicing/chRNA.Expression.Splicing)) %>%
#    add_count(NewAnnotation) %>%
#    filter(n>1000 ) %>%
#    # distinct(NewAnnotation, .keep_all=T) %>% dplyr::select(NewAnnotation, n) %>%
#    # arrange(n) %>% print(n=100)
#    ggplot(aes(x=Ratio, color=NewAnnotation)) +
#    geom_vline(xintercept=0) +
#    stat_ecdf() +
#    coord_cartesian(xlim=c(-2,2)) +
#    # facet_wrap(~NewAnnotation) +
#    theme_bw() +
#    labs(x="log10(polyA/chRNA)", y="ecdf")

#YRI.Samples <- dat$Sample %>%
#    unique() %>%
#    str_replace(".+?_(.+?)$", "\\1") %>% unique()

##ggsave("scratch/ReadLengthMapExperiment.EcdfRatioByAnnotation.pdf", P2, height=7, width=7)

#SummedByDatasetIndID <- dat %>%
#    dplyr::select(Dataset:strand) %>%
#    # head(1000000) %>%
#    inner_join(Annotations) %>%
#    group_by(Dataset,IndID, NewAnnotation) %>%
#    summarise(Count = sum(score) + 0.00) %>%
#    ungroup()

## write_tsv(SummedByDatasetIndID, out_f_longtable)

#DatasetIndIDTotals <- SummedByDatasetIndID %>%
#    group_by(IndID,Dataset) %>%
#    summarise(TotalJuncsInDataset = sum(Count)) %>%
#    ungroup()

#P.boxplots.dat <- SummedByDatasetIndID %>%
#    inner_join(DatasetIndIDTotals) %>%
#    mutate(FractionTotalJuncs = Count / TotalJuncsInDataset) %>%
#    mutate(Dataset = recode(Dataset, "chRNA.Expression.Splicing"="chRNA", "Expression.Splicing"="polyA", "MetabolicLabelled.30min"="30min", "MetabolicLabelled.60min"="60min")) %>%
#    mutate(Dataset = factor(Dataset, levels=c("chRNA","30min", "60min", "polyA")))
#P.boxplots <- ggplot(P.boxplots.dat, aes(x=Dataset, y=FractionTotalJuncs)) +
#      geom_jitter(alpha=0.2, size=0.5) +
#      geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
#      facet_wrap(~NewAnnotation, scales="free_y", labeller = label_wrap_gen(width=14)) +
#      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
#      labs(y=str_wrap("Fraction of splice junction reads", 20), x=NULL) +
#      guides(colour = guide_legend(override.aes4= list(alpha = 1))) +
#      theme(strip.text.x = element_text(size = 4))
##ggsave("scratch/ReadLengthMapExperiment.NewBoxplots.png", P.boxplots, height=15, width=15)

#P.boxplots.relative <- 
#P.boxplots.dat %>%
#    filter(Dataset == "polyA") %>%
#    group_by(NewAnnotation) %>%
#    summarise(MedianPolyA = median(FractionTotalJuncs)) %>%
#    ungroup() %>%
#    inner_join(
#               P.boxplots.dat
#    ) %>%
#    mutate(AbundanceRelativetoMedianPolyA=log2(FractionTotalJuncs/MedianPolyA)) %>%
#    ggplot(aes(x=Dataset, y=AbundanceRelativetoMedianPolyA)) +
#      geom_jitter(alpha=0.2, size=0.5) +
#      geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
#      facet_wrap(~NewAnnotation, scales="free_y", labeller = label_wrap_gen(width=14)) +
#      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
#      labs(y=str_wrap("Relative fraction junctions", 20), x=NULL) +
#      guides(colour = guide_legend(override.aes4= list(alpha = 1))) +
#      theme(strip.text.x = element_text(size = 4))
##ggsave("scratch/ReadLengthMapExperiment.NewBoxplotsNormalized.png", P.boxplots.relative, height=15, width=15)
 
###################

#Old.Data <- fread("SplicingAnalysis/regtools_annotate_combined/Comprehensive.YRI.Sample.LongTable.Counts.tsv.gz")

#Old.SummedByDataset <- Old.Data %>%
#    filter((Sample %in% YRI.Samples) & Dataset %in% c("Expression.Splicing", "chRNA.Expression.Splicing"))
#    dplyr::select(Dataset, IndID=Sample, start, end, strand, chrom=`#Chrom`, score=Count) %>%
#    inner_join(Annotations) %>%
#    group_by(Dataset, chrom, start, end, strand, NewAnnotation) %>%
#    summarise(Count = sum(score) + 0.00) %>%
#    ungroup()

#Old.DatasetTotals <- Old.SummedByDataset %>%
#    group_by(Dataset) %>%
#    summarise(TotalJuncsInDataset = sum(Count)) %>%
#    ungroup()

#Old.ToPlot <- Old.SummedByDataset %>%
#    inner_join(DatasetTotals) %>%
#    mutate(FractionTotalJuncs = Count / TotalJuncsInDataset)

#Old.P <- Old.ToPlot %>%
#    filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
#    dplyr::select(-TotalJuncsInDataset, -Count) %>%
#    pivot_wider(names_from="Dataset", values_from = "FractionTotalJuncs") %>%
#    add_count(NewAnnotation) %>%
#    filter(n>1000 ) %>%
#    group_by(NewAnnotation) %>%
#    sample_n(1000) %>%
#    ungroup() %>%
#    # head() %>% print(width=Inf)
#    ggplot(aes(x=chRNA.Expression.Splicing, y=Expression.Splicing)) +
#    geom_point(alpha=0.2) +
#    geom_abline(color='red') +
#    scale_x_continuous(trans='log10', lim=c(1E-11,1E-3)) +
#    scale_y_continuous(trans='log10', lim=c(1E-11,1E-3) ) +
#    facet_wrap(~NewAnnotation) +
#    theme_bw() +
#    labs(title="Old  mappings")

##ggsave("scratch/ReadLengthMapExperiment.Old.ScatterByAnnotation.png", Old.P, height=15, width=15)


#MappingRatioComparison <- inner_join(
#ToPlot %>%
#    filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
#    dplyr::select(-TotalJuncsInDataset, -Count) %>%
#    pivot_wider(names_from="Dataset", values_from = "FractionTotalJuncs") %>%
#    mutate(Ratio = log10(Expression.Splicing/chRNA.Expression.Splicing)) %>%
#    add_count(NewAnnotation) %>%
#    filter(n>1000 ) %>%
#    dplyr::select(NewAnnotation, Ratio, chrom, start, end, strand),
#Old.ToPlot %>%
#    filter(Dataset %in% c("chRNA.Expression.Splicing", "Expression.Splicing")) %>%
#    dplyr::select(-TotalJuncsInDataset, -Count) %>%
#    pivot_wider(names_from="Dataset", values_from = "FractionTotalJuncs") %>%
#    mutate(Ratio = log10(Expression.Splicing/chRNA.Expression.Splicing)) %>%
#    add_count(NewAnnotation) %>%
#    filter(n>1000 ) %>%
#    dplyr::select(NewAnnotation, Ratio, chrom, start, end, strand),
#by = c("chrom", "start", "end","strand","NewAnnotation"), suffix=c(".New", ".Old")) %>%
#group_by(NewAnnotation) %>%
#sample_n(1000) %>%
#ungroup() %>%
#ggplot(aes(x=Ratio.Old, y=Ratio.New)) +
#geom_abline(color='red') +
#geom_point() +
#theme_bw() +
#facet_wrap(~NewAnnotation) +
#labs(x="log10(polyA/chRNA), Old mapping", y="log10(polyA/chRNA), New mapping")

##ggsave("scratch/ReadLengthMapExperiment.Old.ScatterRatio_MapA_vs_MapB.png", MappingRatioComparison, height=15, width=15)



#Old.SummedByDatasetIndID <- Old.Data %>%
#    dplyr::select(Dataset, IndID=Sample, start, end, strand, chrom=`#Chrom`, score=Count) %>%
#    # head(1000000) %>%
#    inner_join(Annotations) %>%
#    group_by(Dataset,IndID, NewAnnotation) %>%
#    summarise(Count = sum(score) + 0.00) %>%
#    ungroup()

#Old.DatasetIndIDTotals <- Old.SummedByDatasetIndID %>%
#    group_by(IndID,Dataset) %>%
#    summarise(TotalJuncsInDataset = sum(Count)) %>%
#    ungroup()

#Old.P.boxplots.dat <- Old.SummedByDatasetIndID %>%
#    inner_join(Old.DatasetIndIDTotals) %>%
#    mutate(FractionTotalJuncs = Count / TotalJuncsInDataset) %>%
#    mutate(Dataset = recode(Dataset, "chRNA.Expression.Splicing"="chRNA", "Expression.Splicing"="polyA", "MetabolicLabelled.30min"="30min", "MetabolicLabelled.60min"="60min")) %>%
#    mutate(Dataset = factor(Dataset, levels=c("chRNA","30min", "60min", "polyA")))
#ggplot(Old.P.boxplots.dat, aes(x=Dataset, y=FractionTotalJuncs)) +
#      geom_jitter(alpha=0.2, size=0.5) +
#      geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
#      facet_wrap(~NewAnnotation, scales="free_y", labeller = label_wrap_gen(width=14)) +
#      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
#      labs(y=str_wrap("Fraction of splice junction reads", 20), x=NULL) +
#      guides(colour = guide_legend(override.aes4= list(alpha = 1))) +
#      theme(strip.text.x = element_text(size = 4))
##ggsave("scratch/ReadLengthMapExperiment.OldBoxplots.png", Old.P.boxplots, height=15, width=15)
