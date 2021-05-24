#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : PreparePhenotypeTablesFromFeatureCounts
# @created     : Monday May 10, 2021 10:11:08 CDT
#
# @description : from featureCounts output, prepare standardized and
# qqnormalized phenotype table and also table of covariates (genotype PCs,
# expression PCs and sex) for QTLtools
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "featureCounts/polyA.Expression/Counts.txt rename_STAR_alignment_samples ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed ExpressionAnalysis/polyA/ExpressedGeneList.txt QTLs/QTLTools/Expression.Splicing/AllReps.qqnorm.bed.gz QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

featureCounts_FileIn <- args[1]
ColumnRenamerFunction <- args[2]
Gtf_FileIn <- args[3]
genesBed_FileIn <- args[4]
GeneListOut <- args[5]
PhenotypesBedOut_All <- args[6]
PhenotypesBedOut_OnlyFirstReps <- args[7]



### helper "ColumnRenamer" functions to rename the filename column names from featureCounts to the sampleIDs as used in the vcf for qtl calling
rename_STAR_alignment_samples <- function(MyString){
    return(
           str_replace(MyString, "Alignments/STAR_Align/.+?/(.+?)/(\\d+)/Filtered\\.bam", "\\1.\\2")
    )
}

rename_hisat2_alignment_samples <- function(MyString){
    return(
           str_replace(MyString,"Alignments/Hisat2_Align/.+?/(.+?)\\.(\\d+)\\.wasp_filterd\\.markdup\\.sorted\\.bam", "\\1.\\2")
    )
}

# function to help parse gtf file to identify protein coding genes
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}


#Get protein coding features, using same gtf as used in featureCounts so that gene_ids perfectly match
gtf <- read_tsv(Gtf_FileIn, comment="#", n_max=Inf, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )) %>%
    filter(feature=="gene") %>%
    select(attribute)
gtf$gene_id <- unlist(lapply(gtf$attribute, extract_attributes, "gene_id"))
gtf$gene_type <- unlist(lapply(gtf$attribute, extract_attributes, "gene_type"))
ProteinCodingGenes <- filter(gtf, gene_type=="protein_coding") %>% pull(gene_id)

genes_bed <- read_tsv(genesBed_FileIn, col_names=c("Chr", "Start", "End", "Strand", "Geneid", "geneName"), col_types='cnnccc') %>% select(-geneName)

dat <- read_tsv(featureCounts_FileIn, comment = "#") %>%
    rename_with(get(ColumnRenamerFunction), starts_with("Alignments")) %>%
    # mutate(Strand=str_replace(Strand, "^([+-]);.+", "\\1")) %>%
    # mutate(Chr=str_replace(Chr, "^(.+?);.+", "\\1")) %>%
    # mutate(Start=str_replace(Start, "^(.+?);.+", "\\1")) %>%
    # mutate(End=str_replace(End, "^.+?;(.+?)$", "\\1"))
    select(-c("Strand", "Chr", "Start", "End")) %>%
    inner_join(genes_bed, ., by="Geneid") %>%
    mutate(Chr=paste0("chr", Chr))



dat.cpm <- dat %>%
    #Filter for protein coding genes, or non genes (for chip-seq peaks)
    filter(Geneid %in% ProteinCodingGenes  ) %>%
    filter(Chr %in% paste0("chr", 1:22)) %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length")) %>%
    cpm(log=T, prior.count=0.1)


#Filter for top N autosomal genes based on median expression
MedCpm <- sort(apply(dat.cpm, 1, median), decreasing=T)
NGenesToInclude <- 14000
GenesToInclude <- MedCpm[1:NGenesToInclude] %>% names()
dat.cpm.filtered <- dat.cpm[GenesToInclude,] %>% as.matrix()
print(paste("The ", NGenesToInclude, " included genes are a cutoff of", 2**MedCpm[14000], "cpm"))


#Standardize across individuals (rows),
dat.standardized <- dat.cpm.filtered %>% t() %>% scale() %>% t()
#then qqnorm across genes (columns)
dat.qqnormed <- apply(dat.standardized, 2, rankNorm)


Out <- dat %>%
    select(Geneid, Chr, Start, End, Strand) %>%
    inner_join(
               (dat.qqnormed %>% as.data.frame() %>% rownames_to_column("Geneid")),
               by = "Geneid") %>%
    mutate(start= as.numeric(Start)) %>%
    mutate(across(where(is.numeric), round, 5)) %>%
    dplyr::select(`#Chr`=Chr, start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything()) %>%
    arrange(`#Chr`, start)

# Write all samples out
write_tsv(Out, PhenotypesBedOut_All )

Out %>%
    select(1:6, matches("\\.1$")) %>%
    rename_with(~str_remove(., '\\.1$')) %>%
    write_tsv(PhenotypesBedOut_OnlyFirstReps)

#Write out genes included in analysis
data.frame(Geneid = GenesToInclude) %>%
    inner_join(genes_bed, by="Geneid") %>%
    mutate(Chr=paste0("chr", Chr), Score=".") %>%
    select(Chr, Start, End, Geneid, Score, Strand) %>%
    arrange(Chr, Start) %>%
    write_tsv(GeneListOut, col_names=F)

#library(magrittr)

#permute_mat <- function(mat){
#    mat.permuted <- matrix(nrow=nrow(mat), ncol=ncol(mat))
#    for (i in 1:nrow(mat)){
#    mat.permuted[i,] <- sample(mat[i,], size=ncol(mat))
#    }
#    colnames(mat.permuted) <- colnames(mat)
#    rownames(mat.permuted) <- rownames(mat)
#    return(mat.permuted)
#}

#mat <- dat.qqnormed
#mat.permuted <- permute_mat(dat.qqnormed)

#pca <- prcomp(t(mat)) %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>%
#    rownames_to_column("PC")
#pca.permuted <- prcomp(t(mat.permuted)) %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>% rownames_to_column("PC")
#merged <- bind_rows(list(pca=pca, pca.permuted=pca.permuted), .id="mat") %>%
#    mutate(PC=as.numeric(str_replace(PC, "PC", "")))

##GetNumPCs
#merged %>%
#    select(PC, Prop=`Proportion of Variance`, mat) %>%
#    spread(key="mat", value="Prop") %>%
#    filter(pca > pca.permuted) %>% pull(PC) %>% max()


#merged %>%
#    ggplot(aes(x=PC, y=`Proportion of Variance`, color=mat)) +
#    geom_line() +
#    theme_bw()
