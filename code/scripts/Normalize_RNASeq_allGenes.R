library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

chRNACounts_f_in <- "featureCounts/chRNA.Expression/Counts.txt"
polyACounts_f_in <- "featureCounts/polyA.Expression/Counts.txt"
ml30Counts_f_in <- "featureCounts/MetabolicLabelled.30min/Counts.txt"
ml60Counts_f_in <- "featureCounts/MetabolicLabelled.60min/Counts.txt"


rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- "rename_STAR_alignment_samples"

chRNA.genes <- read_tsv(chRNACounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))

polyA.genes <- read_tsv(polyACounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))

ml30.genes <- read_tsv(ml30Counts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))

ml60.genes <- read_tsv(ml60Counts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$'))

subset_YRI_in <- "QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"

subset_YRI <- read_tsv(subset_YRI_in, n_max=Inf) %>%
            select(everything(), -c("#Chr", "start", "end", "pid", "gid", "strand")) %>% colnames()

chRNA.rpkm <- chRNA.genes %>% 
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length", "NA18855")) %>%
    rpkm(log=F, prior.count=0.1, gene.length=chRNA.genes$Length) %>% log1p()


polyA.rpkm <- polyA.genes %>% 
    column_to_rownames("Geneid") %>%
    select(subset_YRI) %>%
    rpkm(log=F, prior.count=0.1, gene.length=polyA.genes$Length) %>% log1p()


ml30.rpkm <- ml30.genes %>% 
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    rpkm(log=F, prior.count=0.1, gene.length=ml30.genes$Length) %>% log1p()

ml60.rpkm <- ml60.genes %>% 
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    rpkm(log=F, prior.count=0.1, gene.length=ml60.genes$Length) %>% log1p()

chRNA.rpkm %>% as.data.frame() %>% rownames_to_column("Geneid") %>% write_tsv('QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstRepsUnstandardized.AllGenes.qqnorm.bed.gz')

polyA.rpkm %>% as.data.frame() %>% rownames_to_column("Geneid") %>% write_tsv('QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstRepsUnstandardized.AllGenes.qqnorm.bed.gz')

ml30.rpkm %>% as.data.frame() %>% rownames_to_column("Geneid") %>% write_tsv('QTLs/QTLTools/MetabolicLabelled.30min/OnlyFirstRepsUnstandardized.AllGenes.qqnorm.bed.gz')

ml60.rpkm %>% as.data.frame() %>% rownames_to_column("Geneid") %>% write_tsv('QTLs/QTLTools/MetabolicLabelled.60min/OnlyFirstRepsUnstandardized.AllGenes.qqnorm.bed.gz')

