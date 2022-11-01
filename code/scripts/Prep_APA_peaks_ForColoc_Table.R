library(tidyverse)

dat <- paste0("QTLs/QTLTools/", c("APA_Nuclear", "APA_Total"),"/OnlyFirstReps.sorted.qqnorm.bed.gz") %>%
  set_names(str_replace(., "QTLs/QTLTools/(.+?)/OnlyFirstReps.sorted.qqnorm.bed.gz", "\\1")) %>%
  lapply(read_tsv) %>%
  bind_rows(.id = "Dataset")

genes <- read_tsv("ExpressionAnalysis/polyA/ExpressedGeneList.txt", col_names = c("chrom", ""))
