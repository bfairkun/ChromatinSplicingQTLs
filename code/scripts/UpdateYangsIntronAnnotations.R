#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : UpdateYangsIntronAnnotations
# @created     : Wednesday Jun 21, 2023 12:02:53 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../data/IntronAnnotationsFromYang.tsv.gz ReferenceGenome/Annotations/gencode.v37.chromasomal.annotation.GeneTypes.tsv.gz ../data/IntronAnnotationsFromYang.Updated.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

IntronAnnotations_f <- args[1]
GeneTypeAnnotations_f <- args[2]
output_f <- args[3]

library(tidyverse)

IntronAnnotations <- read_tsv(IntronAnnotations_f)
GeneTypeAnnotations <- read_tsv(GeneTypeAnnotations_f, col_names=c("gene", "type"))

GeneTypeAnnotations %>% pull(type) %>% unique()

GeneTypeAnnotations %>%
    count(type) %>%
    arrange(desc(n))



NewIntronAnnotations <- IntronAnnotations %>%
    left_join(GeneTypeAnnotations, by=c("gene")) %>%
    separate(SuperAnnotation, into=c("IsAnnotated", "IsProductive"), sep="_", remove=F) %>%
    # distinct(IsAnnotated, IsProductive)
    mutate(SemiSupergroupAnnotations = case_when(
        type == "protein_coding" ~ SemiSupergroupAnnotations,
        TRUE ~ type
                                                 )) %>%
    mutate(SuperAnnotation = case_when(
        type == "protein_coding" ~ SuperAnnotation,
        !type == "protein_coding" & IsAnnotated == "AnnotatedJunc" ~ "AnnotatedJunc_NoncodingGene",
        !type == "protein_coding" & IsAnnotated == "UnannotatedJunc" ~ "UnannotatedJunc_NoncodingJunc"

                                       )) %>%
    dplyr::select(chrom, start, end, strand, NewAnnotation, gene, symbol, SuperAnnotation, SemiSupergroupAnnotations, gene_type=type)

IntronAnnotations %>%
    count(SuperAnnotation)

NewIntronAnnotations %>%
    count(SuperAnnotation)

NewIntronAnnotations %>%
    write_tsv(output_f)
