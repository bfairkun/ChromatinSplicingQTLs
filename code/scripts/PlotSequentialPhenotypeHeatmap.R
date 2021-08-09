#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : PlotSequentialPhenotypeHeatmap
# @created     : Thursday Aug 05, 2021 13:21:44 CDT
#
# @description :
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "ExpressionAnalysis/polyA/ExpressedGeneList.txt scratch/MultiphenotypeSpearmanMat.pdf featureCounts/chRNA.Expressiong/Counts.txt featureCounts/polyA.Expression/Counts.txt featureCounts/MetabolicLabelled.30min/Counts.txt featureCounts/MetabolicLabelled.60min/Counts.txt featureCounts/AtTSS/H3K27AC/Counts.txt featureCounts/AtTSS/H3K4ME3/Counts.txt", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

geneListBedfile <- args[1]
PdfOut <- args[2]
CountFiles <- args[-c(1:2)]

library(tidyverse)
library(edgeR)
library(gplots)
library(scales)
library(magrittr)
library(viridis)

NumbersToColors <- function(NumberVector){
  N <- length(unique(NumberVector))
  Key <- setNames(hue_pal()(N), 1:N)
  return(list(ColorVector=recode(NumberVector, !!!Key), Key=Key))
}

GetCountMatrix <- function(FileIn, genesToInclude){
    CountsFile.df <- read_tsv(FileIn, comment="#") %>%
        filter(Geneid %in% genesToInclude)
    GeneLengths <- CountsFile.df$Length
    Mat <- CountsFile.df %>%
        select(-Chr, -Start, -End, -Strand, -Length) %>%
        column_to_rownames("Geneid") %>%
        DGEList() %>%
        rpkm(log=T, prior.count=0.1, gene.length=GeneLengths)
    Out <- Mat[,sample(1:ncol(Mat), 15)]
    colnames(Out) <- 1:15
    Out %>% as.data.frame() %>% rownames_to_column("Geneid") %>% return()
}


genesToInclude <- read_tsv(geneListBedfile, col_names=c("chr", "start", "stop", "gene", "score", "strand")) %>% pull(gene)

set.seed(1)
mylist <- lapply(CountFiles, GetCountMatrix, genesToInclude)
names(mylist) <- str_replace(CountFiles, "featureCounts/(.+)/Counts.txt", "\\1")

CombinedMat <- bind_rows(mylist, .id="Source") %>% as.tibble() %>%
    gather("SampleNum", "rpkm", -Source, -Geneid) %>%
    unite("Sample", Source, SampleNum ) %>%
    spread(Sample, rpkm) %>%
    column_to_rownames("Geneid")

MyColors <- colnames(CombinedMat) %>%
    str_replace("_\\d+", "") %>%
    factor() %>%
    as.numeric() %>%
    NumbersToColors() %>%
    extract2(1)


pdf(PdfOut)
CombinedMat %>%
    cor(method="spearman") %>%
    heatmap.2(trace="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, ColSideColors = MyColors, RowSideColors=MyColors, col=viridis(30, option="B", direction = 1), labRow=F)
dev.off()
