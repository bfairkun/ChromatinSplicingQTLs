#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : ProcessIRFeatureCounts
# @created     : Sunday May 16, 2021 14:16:30 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "SplicingAnalysis/IR/chRNA.Expression.Splicing/Counts.txt SplicingAnalysis/IR/chRNA.Expression.Splicing/Counts.txt QTLs/QTLTools/chRNA.IR/OnlyFirstReps.qqnorm.bed.gz SplicingAnalysis/IR/chRNA.IR/IR.Ratio.plot.pdf 100000", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


IntronCounts_f_in <- args[1]
GeneCounts_f_in <- args[2]
f_out <- args[3]
plot_out <- args[4]
SampleIRCountsMinimum <- as.numeric(args[5])

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

### helper "ColumnRenamer" functions to rename the filename column names from featureCounts to the sampleIDs as used in the vcf for qtl calling
rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- "rename_STAR_alignment_samples"

dat.introns <- read_tsv(IntronCounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))

dat.genes <- read_tsv(GeneCounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))



dat.introns.norm <- dat.introns %>%
    separate(Geneid, into=c("Geneid", "IntronID"), sep="_") %>%
    gather(key="sample", value="counts", -c("Geneid", "IntronID", "Chr", "Start", "End", "Strand", "Length")) %>%
    group_by(sample) %>%
    mutate(sum_counts = sum(counts)) %>%
    # summarize(sum_counts=sum(counts)) %>% as.data.frame()
    ungroup() %>%
    filter(sum_counts > SampleIRCountsMinimum) %>%
    mutate(cpk=counts/Length*1000) %>%
    select(-sum_counts)

All.Samples <- colnames(dat.introns %>% select(-(1:6)))
Samples <- dat.introns.norm$sample %>% unique()
length(Samples)

print("Samples filtered out:")
setdiff(All.Samples, Samples)

dat.introns.norm.byhostgene <-
    dat.genes %>%
    select(Geneid, Length, all_of(Samples)) %>%
    gather(key="sample", value="counts", -c("Geneid", "Length")) %>%
    mutate(cpk=counts/Length*1000) %>%
    select(-Length) %>%
    right_join(dat.introns.norm, by=c("sample", "Geneid"), suffix=c(".hostgene", ".intron")) %>%
    mutate(IR.Ratio = cpk.intron/cpk.hostgene)

head(dat.introns.norm.byhostgene)

Plot <- dat.introns.norm.byhostgene %>%
    group_by(IntronID) %>%
    summarize(
           med_ratio = median(IR.Ratio, na.rm=T)
    ) %>%
ggplot(aes(x=med_ratio)) +
    geom_histogram() +
    scale_x_continuous(trans="log10") +
    xlab("IntronicCoverage/HostGeneCoverage")+
    theme_bw()

ggsave(plot_out, Plot, height=3, width=3)
# ggsave("scratch/Plot.Introns.IR.plyA.pdf", Plot, height=3, width=3)


#Filter where all samples < 1E-2, or all samples > 10
Filtered.out <- dat.introns.norm.byhostgene %>%
    group_by(IntronID) %>%
    mutate(
           min_ratio = min(IR.Ratio),
           max_ratio = max(IR.Ratio),
           med_ratio = median(IR.Ratio)
    ) %>%
    ungroup() %>%
    filter(min_ratio <= 10) %>%
    filter(max_ratio >= 1E-3) %>%
    mutate(IR.Ratio.Fixed = case_when(
                                      counts.intron == 0 ~ med_ratio,
                                      is.na(IR.Ratio) ~ med_ratio,
                                      # med_ratio == 0 ~ NA_complex_,
                                      TRUE ~ IR.Ratio,
                                      )) %>%
    mutate(IR.Ratio.Fixed = log(IR.Ratio.Fixed + 1E-3)) %>%
    select(Chr, Start, End, IntronID, Geneid, Strand, sample, IR.Ratio.Fixed) %>%
    spread(key="sample", value="IR.Ratio.Fixed") %>%
    drop_na()

BedInfo <- Filtered.out %>%
    select(`#Chr`=Chr, Start, End, Geneid, IntronID, Strand)

m <- Filtered.out %>%
    select(-c("Chr", "Start", "End", "Geneid", "Strand")) %>%
    column_to_rownames("IntronID") %>%
    as.matrix() %>%
    #standardized across rows
    t() %>% scale() %>% t()
m.qqnorm <- m[rowSums(is.na(m)) != ncol(m), ] %>%
    #inverse normalize across columns
    apply(2, RankNorm)

print("How many introns out:")
print(dim(m.qqnorm))

m.qqnorm %>% as.data.frame() %>%
    rownames_to_column("IntronID") %>%
    as_tibble() %>%
    left_join(BedInfo, by="IntronID") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    select(`#Chr`, start=Start, end=End, pid=IntronID, gid=Geneid, strand=Strand, everything()) %>%
    write_tsv(f_out)

