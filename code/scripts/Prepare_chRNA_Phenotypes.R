library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

GeneCounts_f_in <- "featureCounts/chRNA.Expression/Counts.txt"
GeneCounts_eRNA_in <- "featureCounts/chRNA.Expression_eRNA/Counts.txt"
GeneCounts_cheRNA_in <- "featureCounts/chRNA.Expression_cheRNA/Counts.txt"
GeneCounts_lncRNA_in <- "featureCounts/chRNA.Expression_lncRNA/Counts.txt"
GeneCounts_snoRNA_in <- "featureCounts/chRNA.Expression_snoRNA/Counts.txt"
Genes_bed_f_in <- "ExpressionAnalysis/polyA/ExpressedGeneList.txt" 

rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

ColumnRenamerFunction <- "rename_STAR_alignment_samples"

gene.list <- read_tsv(Genes_bed_f_in, col_names=c("Chr", "Start", "End", "Geneid", "score", "Strand"))

dat.genes <- read_tsv(GeneCounts_f_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) #%>%
#     filter(Geneid %in% gene.list$Geneid)


dat.genes.eRNA <- read_tsv(GeneCounts_eRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) #%>%
#     filter(Geneid %in% gene.list$Geneid)

dat.genes.cheRNA <- read_tsv(GeneCounts_cheRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 

dat.genes.lncRNA <- read_tsv(GeneCounts_lncRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 

dat.genes.snoRNA <- read_tsv(GeneCounts_snoRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 


X <- rbind(dat.genes, dat.genes.eRNA, dat.genes.cheRNA)


dat.cpm <- X %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)


protein_coding = dat.cpm[rownames(dat.cpm) %in% gene.list$Geneid, ]
eRNA_ <- dat.cpm[rownames(dat.cpm) %in% dat.genes.eRNA$Geneid, ]
lncRNA_ <- dat.cpm[rownames(dat.cpm) %in% dat.genes.lncRNA$Geneid, ]
cheRNA_ <- dat.cpm[rownames(dat.cpm) %in% dat.genes.cheRNA$Geneid, ]
snoRNA_ <- dat.cpm[rownames(dat.cpm) %in% dat.genes.snoRNA$Geneid, ]


snoRNA <- snoRNA_[apply(snoRNA_, 1, median) > -7,]
cheRNA <- cheRNA_[apply(cheRNA_, 1, median) > -7,]
lncRNA <- lncRNA_[apply(lncRNA_, 1, median) > -7,]
eRNA <- eRNA_[apply(eRNA_, 1, median) > -7,]



eRNA.standardized <- eRNA %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
eRNA.qqnormed <- apply(eRNA.standardized, 2, RankNorm)

snoRNA.standardized <- snoRNA %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
snoRNA.qqnormed <- apply(snoRNA.standardized, 2, RankNorm)

lncRNA.standardized <- lncRNA %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
lncRNA.qqnormed <- apply(lncRNA.standardized, 2, RankNorm)

cheRNA.standardized <- cheRNA %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
cheRNA.qqnormed <- apply(cheRNA.standardized, 2, RankNorm)

protein_coding.standardized <- protein_coding %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
protein_coding.qqnormed <- apply(protein_coding.standardized, 2, RankNorm)


genesBed_FileIn <- "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
genes_bed <- read_tsv(genesBed_FileIn, col_names=c("Chr", "Start", "End", "Strand", "Geneid", "geneName"), col_types='cnnccc') %>% select(-geneName)




lncRNA_bed <- data.frame(Geneid = rownames(lncRNA.qqnormed)) %>%
    inner_join(genes_bed, by="Geneid") %>%
    mutate(Chr=paste0("chr", Chr), Score=".") %>%
    select(Chr, Start, End, Geneid, Score, Strand) %>%
    arrange(Chr, Start)


snoRNA_bed <- data.frame(Geneid = rownames(snoRNA.qqnormed)) %>%
    inner_join(genes_bed, by="Geneid") %>%
    mutate(Chr=paste0("chr", Chr), Score=".") %>%
    select(Chr, Start, End, Geneid, Score, Strand) %>%
    arrange(Chr, Start)


eRNA_bed <- dat.genes.eRNA %>% select(Chr, Start, End, Strand, Geneid)
eRNA_bed <- eRNA_bed[eRNA_bed$Geneid %in% rownames(eRNA.qqnormed), ]



cheRNA_bed <- dat.genes.cheRNA %>% select(Chr, Start, End, Strand, Geneid)
cheRNA_bed <- cheRNA_bed[cheRNA_bed$Geneid %in% rownames(cheRNA.qqnormed), ]


protein_coding.Out <- gene.list %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (protein_coding.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)



cheRNA.Out <- cheRNA_bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (cheRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)



eRNA.Out_ <- eRNA_bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (eRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)

eRNA.Out <- rbind(eRNA.Out_[eRNA.Out_$strand == '-',], eRNA.Out_[eRNA.Out_$strand == '+',])

snoRNA.Out <- snoRNA_bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (snoRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)


lncRNA.Out <- lncRNA_bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (lncRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)




write_tsv(protein_coding.Out, "QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz")
write_tsv(snoRNA.Out, "QTLs/QTLTools/chRNA.Expression_snoRNA/OnlyFirstReps.qqnorm.bed.gz")
write_tsv(lncRNA.Out, "QTLs/QTLTools/chRNA.Expression_lncRNA/OnlyFirstReps.qqnorm.bed.gz")
write_tsv(cheRNA.Out, "QTLs/QTLTools/chRNA.Expression_cheRNA/OnlyFirstReps.qqnorm.bed.gz")
write_tsv(eRNA.Out, "QTLs/QTLTools/chRNA.Expression_eRNA/OnlyFirstReps.qqnorm.bed.gz")

