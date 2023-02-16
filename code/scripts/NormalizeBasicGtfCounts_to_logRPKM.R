args <- commandArgs(trailingOnly=TRUE)

counts.fn <- args[1]
out.fn <- args[2]
standardized.fn <- 'QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.sorted.qqnorm.bed.gz'

print('loading libraries')
library(tidyverse)
library(edgeR)

print('finished loading libraries')

# Define funciton to read RPKM table and rename

read_and_rename <- function(fn){
  read_tsv(fn, comment="#") %>%
    dplyr::select(Geneid, Length, contains("/1/Filtered.bam")) %>%
    rename_at(vars(contains("/1/Filtered.bam")), ~str_replace(., 
                "Alignments/STAR_Align/(.+?)/(.+?)/1/Filtered.bam", "\\2")) %>%
    return()
}

standardized <- read_tsv(standardized.fn)

# Read featureCounts table

counts <- read_and_rename(counts.fn) %>%
  dplyr::select(Geneid, Length, !matches(".+?_\\d+")) %>%
  filter(Geneid %in% standardized$pid)

RPKM <- counts %>%
  dplyr::select(-Length) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  rpkm(prior.count=0.1, log=T, gene.length=counts$Length) %>%
  as.data.frame() %>%
  rownames_to_column("pid")

# Merge with BED for top 14000 genes

standardized %>%
  dplyr::select(1:6) %>%
  inner_join(
    RPKM
  ) %>%
  write_tsv(out.fn)