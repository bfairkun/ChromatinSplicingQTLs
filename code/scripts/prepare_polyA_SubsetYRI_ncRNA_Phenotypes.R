library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

subset_YRI_in <- "QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"
ncRNA_RPKM_in <- "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.RPKM.bed.gz"
Genes_bed_f_in <- "ExpressionAnalysis/polyA/ExpressedGeneList.txt" 

subset_YRI <- read_tsv(subset_YRI_in, n_max=Inf) %>%
            select(everything(), -c("#Chr", "start", "end", "pid", "gid", "strand")) %>% colnames()


bed_rpkm <- read_tsv(ncRNA_RPKM_in, n_max=Inf) #%>%
    #rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    #select(1:6, matches("\\.1$")) #%>%
    #rename_with(~str_remove(., '\\.1$')) #%>%

bed <- bed_rpkm %>% select(1:6)

rownames(bed) <- bed$pid


dat.cpm <- bed_rpkm %>% select(subset_YRI)

rownames(dat.cpm) <- bed$pid


ncRNA.standardized <- dat.cpm %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
ncRNA.qqnormed <- apply(ncRNA.standardized, 2, RankNorm)

bed<- bed[rownames(ncRNA.qqnormed), ]
rownames(bed) <- rownames(ncRNA.qqnormed)
ncRNA.Out <- bed %>%
    inner_join(
               (ncRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("gid")),
               by = "gid") %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    arrange(`#Chr`, start)


                
write_tsv(ncRNA.Out, "QTLs/QTLTools/polyA.Expression_ncRNA.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz")
