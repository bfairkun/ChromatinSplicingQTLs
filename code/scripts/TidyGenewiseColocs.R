#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : TidyGenewiseColocs2
# @created     : Saturday Jan 08, 2022 13:30:05 CST
#
# @description : tidy the hyprcoloc results for colocalized traits and merge-in summary stats (eg, beta values) for the top SNP for each cluster. the output text file does not contain any lines for  non-colocalized traits
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../output/hyprcoloc_results/ForColoc/hyprcoloc.results.txt.gz hyprcoloc/LociWiseSummaryStatsInput/ForColoc/ scratch/test.hyprcoloc.tidy.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


library(tidyverse)
library(data.table)

f_in <- args[1]
summary_stats_dir <- args[2]
f_out <- args[3]

sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by_(.dots = grps)
}

GetSummaryStats <- function(df){
    SummaryStats <- fread(paste0(summary_stats_dir, df$Locus[1], ".txt.gz")) %>%
        filter(snp %in% unique(df$topSNP)) %>%
        mutate(phenotype_class = str_replace(source_file, "QTLs/QTLTools/(.+?)/.+", "\\1")) %>%
        mutate(phenotype_full = paste0(phenotype_class, ";", phenotype)) %>%
        select(snp, beta, beta_se, p, Locus=gwas_locus, phenotype_full) %>%
        right_join(df, by=c("phenotype_full"="ColocalizedTraits", "snp"="topSNP", "Locus")) %>%
        return()
}

dat <- read_tsv(f_in, col_names = c("Locus", "iteration", 'ColocalizedTraits', 'ColocPr', 'RegionalPr', "topSNP", "TopSNPFinemapPr", "DroppedTrait"), skip=1) %>%
  filter(!ColocalizedTraits == "None") %>%
  separate_rows(ColocalizedTraits, sep = ', ') %>%
  group_by(Locus)
  # sample_n_groups(100)

Split.list <- setNames(group_split(dat), deframe(group_keys(dat)))

Split.list.w.summarystats <- lapply(Split.list, GetSummaryStats)

bind_rows(Split.list.w.summarystats, .id="Locus") %>%
    select(-DroppedTrait) %>%
    write_tsv(f_out)

