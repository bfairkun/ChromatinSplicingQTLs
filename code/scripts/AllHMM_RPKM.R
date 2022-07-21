library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)



GeneCounts_f_in <- "NonCodingRNA/Expression_HMM/chRNA.Expression_featureCounts/Counts.txt"

rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- "rename_STAR_alignment_samples"


dat <- read_tsv(GeneCounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 


dat.cpm <- dat %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)


dat.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
dat.qqnormed <- apply(dat.standardized, 2, RankNorm)




bed <- dat %>% mutate(Score=".") %>% select(Chr, Start, End, Geneid, Score, Strand)
bed <- bed[bed$Geneid %in% rownames(dat.qqnormed), ]



Out <- bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)


write_tsv(Out, "NonCodingRNA/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz")



RPKM.Out <- bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.cpm %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)


write_tsv(RPKM.Out, "NonCodingRNA/Expression_HMM/OnlyFirstReps.RPKM.bed.gz")
