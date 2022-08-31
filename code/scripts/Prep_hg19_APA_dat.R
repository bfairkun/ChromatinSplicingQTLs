#!/usr/bin/env Rscript

## Script to take Briana's nuclear and total APA site raw leafcutter count
## ratios and her permutation results (both of which are in hg19) and create a
## standardized/normalized phenotype table, which I will liftover later with a
## snakemake rule

library(tidyverse)
library(data.table)
library(magrittr)
library(RNOmni)

permutation.pass.dat <-
Sys.glob("../data/APA_dat/APApeak_BrianasPermutationTestResults.*.tsv.gz") %>%
  set_names(str_replace(., "../data/APA_dat/APApeak_BrianasPermutationTestResults.(.+?).tsv.gz", "\\1")) %>%
  lapply(read_delim, delim=' ') %>%
  bind_rows(.id="Fraction")

count.ratios.dat <-
  Sys.glob("../data/APA_dat/APApeak_Phenotype_GeneLocAnno.*.tsv.gz") %>%
  set_names(str_replace(., "../data/APA_dat/APApeak_Phenotype_GeneLocAnno.(.+?).tsv.gz", "\\1")) %>%
  lapply(read_delim, delim=' ') %>%
  bind_rows(.id="Fraction") %>%
  #Only keep the numerator of the fraction
  mutate(across(-c(1:2), ~ as.numeric(str_remove(.x, "/.+")))) %>%
  mutate(chrom = str_replace(chrom, "chr(.+?):(.+?):(.+?):(.+?)_(.+?)_([-+])_(peak.+)$", "\\1;\\2;\\3;\\4;\\5;\\6;\\7")) %>%
  separate(chrom, into=c("Chrom", "start", "stop", "Gene", "feature", "strand", "peak"), sep=";", convert=T) %>%
  mutate(pid = str_glue("{Chrom}:{start}:{stop}:{Gene}_{feature}_{strand}_{peak}")) %>%
  filter(!feature=="008559_end")

Filtered.ratios <- count.ratios.dat %>%
  as_tibble() %>%
  inner_join(
    permutation.pass.dat %>%
      select(pid, Fraction),
    by=c("Fraction", "pid")) %>%
  select(-pid) %>%
  gather("Sample", "PeakCounts", -c(1:8)) %>%
  group_by(Fraction, Gene, Sample) %>%
  mutate(Ratio = PeakCounts/sum(PeakCounts)) %>%
  mutate(pid = paste(Gene, feature, peak, sep = "_")) %>%
  ungroup() %>%
  group_by(Fraction, Gene) %>%
  mutate(sd = sd(Ratio)) %>%
  mutate(Median = median(Ratio, na.rm = T)) %>%
  mutate(Ratio = case_when(
    !is.na(Ratio) ~ Ratio,
    TRUE ~ Median
  )) %>%
  ungroup()

#Filter, standardize, normalize, and write out nuclear data
Filtered.ratios %>%
  filter(Fraction == "Nuclear") %>%
  filter(sd > 0) %>%
  select(pid, Sample, Ratio) %>%
  spread("Sample", "Ratio") %>%
  column_to_rownames("pid") %>%
  t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix() %>%
  apply(2, RankNorm) %>%
  as.data.frame() %>%
  rownames_to_column("pid") %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  inner_join(
    Filtered.ratios %>%
      distinct(pid, .keep_all=T) %>%
      mutate(Chrom=paste0("chr", Chrom)) %>%
      select(`#Chr`=Chrom, start, end=stop, pid, gid=Gene, strand),
    by="pid"
  ) %>%
  select(`#Chr`, start, end, pid, gid, strand, everything()) %>%
  arrange(`#Chr`, start) %>%
  write_tsv("APA_Processing/PhenotypeTables/Nuclear.hg19.bed")

#Filter, standardize, normalize, and write out Total data
Filtered.ratios %>%
  filter(Fraction == "Total") %>%
  filter(sd > 0) %>%
  select(pid, Sample, Ratio) %>%
  spread("Sample", "Ratio") %>%
  column_to_rownames("pid") %>%
  t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix() %>%
  apply(2, RankNorm) %>%
  as.data.frame() %>%
  rownames_to_column("pid") %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  inner_join(
    Filtered.ratios %>%
      distinct(pid, .keep_all=T) %>%
      mutate(Chrom=paste0("chr", Chrom)) %>%
      select(`#Chr`=Chrom, start, end=stop, pid, gid=Gene, strand),
    by="pid"
  ) %>%
  select(`#Chr`, start, end, pid, gid, strand, everything()) %>%
  arrange(`#Chr`, start) %>%
  write_tsv("APA_Processing/PhenotypeTables/Total.hg19.bed")

