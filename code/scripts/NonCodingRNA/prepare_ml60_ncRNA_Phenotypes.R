library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

GeneCounts_f_in <- "featureCounts/MetabolicLabelled.60min/Counts.txt"
GeneCounts_ncRNA_in <- "featureCounts/MetabolicLabelled.60min_ncRNA/Counts.txt"
GeneCounts_lncRNA_in <- "featureCounts/MetabolicLabelled.60min_annotated_ncRNA/Counts.txt"
Genes_bed_f_in <- "ExpressionAnalysis/polyA/ExpressedGeneList.txt" 
annotation_f_in <- "NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz"

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

dat.genes.ncRNA <- read_tsv(GeneCounts_ncRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 

dat.genes.lncRNA <- read_tsv(GeneCounts_lncRNA_in, comment = "#", n_max=Inf) %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) 

# dat.genes.snoRNA <- read_tsv(GeneCounts_snoRNA_in, comment = "#", n_max=Inf) %>%
#     rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
#     select(1:6, matches("\\.1$")) %>%
#     rename_with(~str_remove(., '\\.1$')) 


X <- rbind(dat.genes, dat.genes.ncRNA)




annot <- read_tsv(annotation_f_in, )

x <- apply(annot[annot$lncRNA != '.','lncRNA'], 2, function(x) c(strsplit(x, "\\|")))
lncRNA_ <- do.call(c, unlist(x, recursive=FALSE))
           
           
# y <- apply(annot[annot$pseudogene != '.','pseudogene'], 2, function(x) c(strsplit(x, "\\|")))
# pseudogene_ <- do.call(c, unlist(y, recursive=FALSE))

# snoRNA_ <- rownames(annot[annot$snoRNA != '.',])
           
#dat.cpm <- dat.cpm %>% as.data.frame() %>%
#  filter(!rownames(dat.cpm) %in% c(lncRNA_, pseudogene_)) %>% as.matrix()


dat.matrix <- X %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length"))
    

dat.matrix.renamed <- dat.matrix %>%
    as.data.frame() %>%
    filter(!rownames(dat.matrix) %in% lncRNA_) %>% #c(lncRNA_, pseudogene_, snoRNA_)) %>% 
    as.matrix() 


dat.cpm <- dat.matrix.renamed %>% # dat.matrix.expressed %>%
    cpm(log=T, prior.count=0.1)

print('created dat.cpm')
protein_coding = dat.cpm[rownames(dat.cpm) %in% gene.list$Geneid, ]



print('selected protein coding')
ncRNA_names <- c(dat.genes.lncRNA$Geneid, dat.genes.ncRNA$Geneid)#, dat.genes.snoRNA$Geneid)
ncRNA.dat <- dat.cpm[rownames(dat.cpm) %in% ncRNA_names, ]
ncRNA <- ncRNA.dat[apply(exp(ncRNA.dat), 1, quantile, probs=0.9) >= 1e-4,]
print('ncRNA')
ncRNA.standardized <- ncRNA %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
ncRNA.qqnormed <- apply(ncRNA.standardized, 2, RankNorm)
print('qqnormed')
protein_coding.standardized <- protein_coding %>% t() %>% scale() %>% t() %>% as.data.frame() %>% drop_na() %>% as.matrix()
protein_coding.qqnormed <- apply(protein_coding.standardized, 2, RankNorm)

genesBed_FileIn <- "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
genes_bed <- read_tsv(genesBed_FileIn, col_names=c("Chr", "Start", "End", "Strand", "Geneid", "geneName"), col_types='cnnccc') %>% select(-geneName)


lncRNA_bed <- data.frame(Geneid = rownames(ncRNA.qqnormed)[(rownames(ncRNA.qqnormed) %in% dat.genes.lncRNA$Geneid)]) %>%
    inner_join(genes_bed, by="Geneid") %>%
    mutate(Chr=paste0("chr", Chr), Score=".") %>%
    select(Chr, Start, End, Geneid, Score, Strand) %>%
    arrange(Chr, Start)


# snoRNA_bed <- data.frame(Geneid = rownames(ncRNA.qqnormed)[(rownames(ncRNA.qqnormed) %in% dat.genes.snoRNA$Geneid)]) %>%
#     inner_join(genes_bed, by="Geneid") %>%
#     mutate(Chr=paste0("chr", Chr), Score=".") %>%
#     select(Chr, Start, End, Geneid, Score, Strand) %>%
#     arrange(Chr, Start)



ncRNA_bed <- dat.genes.ncRNA %>% mutate(Score=".") %>% select(Chr, Start, End, Geneid, Score, Strand)
ncRNA_bed <- ncRNA_bed[ncRNA_bed$Geneid %in% rownames(ncRNA.qqnormed), ]




protein_coding.Out <- gene.list %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (protein_coding.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)


bed <- rbind(lncRNA_bed, ncRNA_bed)#, snoRNA_bed)


ncRNA.Out <- bed %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (ncRNA.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    # mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start=Start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)  %>% as.data.frame()


write_tsv(ncRNA.Out, "QTLs/QTLTools/MetabolicLabelled.60min_ncRNA/OnlyFirstReps.qqnorm.bed.gz")

gene.length = X[X$Geneid %in% rownames(dat.matrix.renamed),'Length']
dat.rpkm <- rpkm(dat.matrix.renamed, gene.length=gene.length$Length, prior.count=0.1)
           
RPKM.Out <- cbind(GeneID = rownames(dat.rpkm), dat.rpkm) %>% as.data.frame()
rownames(RPKM.Out) <- 1:nrow(RPKM.Out) 


write_tsv(RPKM.Out, "RPKM_tables/MetabolicLabelled.60min.RPKM.bed.gz", )

