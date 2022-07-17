library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)


dat.ncRNA <- read_tsv("QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.RPKM.filtered.bed.gz")


dat.cpm <- dat.ncRNA %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length"))


ncRNA.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
ncRNA.qqnormed <- apply(ncRNA.standardized, 2, RankNorm)

ncRNA_bed <- dat.ncRNA %>% mutate(Score=".") %>% select(Chr, Start, End, Geneid, Score, Strand)
ncRNA_bed <- ncRNA_bed[ncRNA_bed$Geneid %in% rownames(ncRNA.qqnormed), ]


ncRNA.Out <- ncRNA_bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (ncRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)


write_tsv(ncRNA.Out, "NonCodingRNA/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz")
