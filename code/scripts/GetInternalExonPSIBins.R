#!/usr/bin/env Rscript


ExonExpressionFileIn <- "featureCounts/AtInternalExons/chRNA.Expression/Counts.txt"
GeneExpressionFileIn <- "featureCounts/chRNA.Expression/Counts.txt"

library(tidyverse)
library(edgeR)

exons <- read_tsv(ExonExpressionFileIn, comment = "#") %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "Alignments/STAR_Align/chRNA.Expression.Splicing/(.+?)/(.+?)/Filtered.bam", "\\1.\\2")) %>%
  select(1:6, ends_with('.1')) %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "(.+?)\\.1$", "\\1"))

ExpressedGenesBed <- read_tsv("ExpressionAnalysis/polyA/ExpressedGeneList.txt", col_names = c("Chr", "start", "stop", "gene", "score", "strand"))
ExpressedGenes <- ExpressedGenesBed$gene

genes <- read_tsv(GeneExpressionFileIn, comment = "#") %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "Alignments/STAR_Align/chRNA.Expression.Splicing/(.+?)/(.+?)/Filtered.bam", "\\1.\\2")) %>%
  select(1:6, ends_with('.1')) %>%
  rename_at(.vars=vars(-(1:6)), .funs=~str_replace(.x, "(.+?)\\.1$", "\\1")) %>%
  filter(Geneid %in% ExpressedGenes)

genes.rpkm <- genes %>%
  select(-(2:6)) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  rpkm(prior.count = 0.1, gene.length=genes$Length)

# Write genes by chRNA expression RPKM quartile
data.frame(MedRPKM = (genes.rpkm %>% apply(1, median))) %>%
  rownames_to_column("gene") %>%
  inner_join(ExpressedGenesBed, by="gene") %>%
  mutate(ExpressionQuartile = paste0("Quartile", ntile(MedRPKM*-1, 4))) %>%
  group_by(ExpressionQuartile) %>%
  arrange(Chr, start, stop) %>%
  select(Chr, start, stop, gene, MedRPKM, strand) %>%
  group_walk(~ write_tsv(.x, paste0("ExpressionAnalysis/chRNAExpression_", .y$ExpressionQuartile, ".bed"), col_names = F))

# MedExonCov <- data.frame(MedExonCoverage =
#              (exons.rpkm %>% apply(1, median))) %>%
#   rownames_to_column("ExonID") %>%
#   separate(ExonID, into=c("Geneid", "ExonNum"), sep="_") %>%
#   filter(Geneid %in% ExpressedGenes) %>%
#   mutate(ExonNum = as.numeric(ExonNum))

SkippedReadCountData <- list.files("SplicingAnalysis/Annotations/InternalExons.noOverlapping.SkippedReads", pattern="*_1.bed.gz", full.names = T) %>%
  setNames(str_replace(., ".+_(NA.+?)_1.bed.gz", "\\1")) %>%
  lapply(read_tsv, col_names = c("Chrom", "start", "stop", "name", "score", "strand", "SkippedCount")) %>%
  bind_rows(.id = "sample") %>%
  select(sample, name, SkippedCount)

Merged.dat <- exons %>%
  select(-Length) %>%
  gather(key = "sample", value="IncludedReadCount", -Geneid, -Chr, -Start, -End, -Strand) %>%
  inner_join(SkippedReadCountData, by=c("sample", "Geneid"="name")) %>%
  separate(Geneid, into=c("HostGene", "ExonNum"), sep="_", remove = F) %>%
  filter(HostGene %in% ExpressedGenes)


Out <- Merged.dat %>%
  mutate(PSI = IncludedReadCount/(IncludedReadCount+SkippedCount)) %>%
  group_by(Geneid) %>%
  mutate(MedPSI = median(PSI)) %>%
  ungroup() %>%
  distinct(Geneid, .keep_all = T) %>%
  drop_na() %>%
  mutate(bin = cut(MedPSI, breaks = c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.1), labels=c("0-20", "20-40", "40-60", "60-80", "80-100"))) %>%
  mutate(name = paste(Geneid, MedPSI, sep="_"))




Out %>%
  # count(bin)
  group_by(bin) %>%
  sample_n(pmin(n(), 1000)) %>%
  ungroup() %>%
  select(Chr, Start, End, name, MedPSI, Strand, bin) %>%
  group_by(bin) %>%
  arrange(Chr, Start, End) %>%
  group_walk(~ write_tsv(.x, paste0("SplicingAnalysis/Annotations/InternalExons.PSI.", .y$bin, ".bed"), col_names = F))

  # ggplot(aes(x=MedPSI)) +
  # geom_histogram() +
  # coord_cartesian(ylim=c(0,1000)) +
  # theme_bw()


