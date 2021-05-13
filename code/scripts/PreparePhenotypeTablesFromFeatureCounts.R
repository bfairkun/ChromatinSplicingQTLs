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
                 "featureCounts/polyA.Expression/Counts.txt rename_STAR_alignment_samples ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf scratch/genes.bed", what='character')
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
PhenotypesBedOut <- args[6]


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

genes_bed <- read_tsv(genesBed_FileIn, col_names=c("Chr", "Start", "End", "Strand", "Geneid", "geneName")) %>% select(-geneName)

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
    dplyr::select(`#Chr`=Chr, start, end=End, pid=Geneid, gid=Geneid, strand=Strand, everything())


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
