#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

A <- fread("SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz") %>%
  column_to_rownames("V1") %>% as.matrix()

sum <- data.frame(SumJuncCount = apply(A, MARGIN = 1, FUN = sum),
                  Intron = rownames(A)) %>%
  as_tibble() %>%
  separate(Intron, into=c("Chrom", "start", "end", "cluster"), sep = ":", convert = T) %>%
  mutate(strand = str_extract(cluster, "[+-]")) %>%
  mutate(Acceptor = case_when(
    strand == "+" ~ paste(Chrom, end),
    strand == "-" ~ paste(Chrom, start)
  )) %>%
  mutate(Donor = case_when(
    strand == "+" ~ paste(Chrom, start),
    strand == "-" ~ paste(Chrom, end)
  )) %>%
  add_count(cluster, Acceptor, name="NumTimesAcceptorInCluster") %>%
  add_count(cluster, Donor, name="NumTimesDonorInCluster")

PotentialSkippingAcceptors <- sum %>%
  filter(NumTimesAcceptorInCluster > 1) %>% pull(Acceptor)
PotentialSkippingDonors <- sum %>%
  filter(NumTimesDonorInCluster > 1) %>% pull(Donor)
PotentialSkippingEvents <- sum %>%
  filter((Acceptor %in% PotentialSkippingAcceptors) & (Donor %in% PotentialSkippingDonors)) %>%
  mutate(Int = paste(Chrom, start, end, strand))

sum %>%
  filter(cluster %in% PotentialSkippingEvents$cluster) %>%
  mutate(Int = paste(Chrom, start, end, strand)) %>%
  mutate(SkippedIntron = Int %in% PotentialSkippingEvents$Int) %>%
  group_by(cluster) %>%
  top_n(3, SumJuncCount) %>%
  head(30) %>%
  group_by(cluster) %>%
  mutate(n = min_rank(SumJuncCount)) %>%
  ggplot(aes(y=n)) +
  geom_segment(aes(x = start, xend=end, yend=n, color=SkippedIntron)) +
  # geom_text(aes(x=mean(start, end), label=SumJuncCount)) %>%
  facet_wrap(~cluster, scales = "free") %>%
  mutate()

