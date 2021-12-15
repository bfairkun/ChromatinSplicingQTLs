#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : hyprcoloc_scratch
# @created     : Wednesday Aug 04, 2021 15:12:12 CDT
#
# @description : for loop thru loci attempting to colocalize mol QTL
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 " scratch/testfilelist.txt  scratch/test.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

# SummaryStatsInList FileOut

library(data.table)
library(tidyverse)
library(hyprcoloc)

FilelistIn <- args[1]
FileOut <- args[2]

Loci<- read_tsv(FilelistIn, col_names=c("f")) %>%
    filter(file.exists(f)) %>%
    mutate(Loci_names = str_replace(f, "hyprcoloc/LociWiseSummaryStatsInput/ForColoc/(.+?)\\.txt\\.gz", "\\1")) %>%
    select(Loci_names, f) %>%
    deframe()


HyprcolocResults.list <- vector("list", length(Loci))

for (i in seq_along(Loci)){

    TestLocus <- Loci[1]
    TestLocus <- Loci[i]
    print(paste0("running loci", i, " : ", names(TestLocus)))


    SummaryStats.filtered <- fread(TestLocus) %>%
        filter(gwas_locus == names(TestLocus)) %>%
        mutate(phenotype_class = str_replace(source_file, "QTLs/QTLTools/(.+?)/(.+?)\\.txt\\.gz", "\\1")) %>%
        unite(phenotype, phenotype_class, phenotype, sep=";") %>%
        select(snp, phenotype, beta, beta_se, p) %>%
        distinct(.keep_all=T)


    # Make snp x phenotype matrix of betas
    betas <- SummaryStats.filtered %>%
        as_tibble() %>%
        select(phenotype, snp, beta) %>%
        distinct(phenotype, snp, .keep_all=T) %>%
        spread(phenotype, beta) %>%
        drop_na() %>%
        column_to_rownames("snp") %>%
        as.matrix()

    ## Consider computing LD matrix as optional input for hyprcoloc
    # snps <- betas %>%
    #     rownames() %>%
    #     as.data.frame() %>%
    #     separate(".", into=c("chrom", "pos", "ref", "alt"), convert=T)

    # Min <- min(snps$pos)
    # Max <- max(snps$pos)
    # Chrom <- snps$chrom[1]
    # CMD <- paste0("vcftools --hap-r2 --gzvcf Genotypes/1KG_GRCh38/", Chrom, ".vcf.gz --chr chr", Chrom, " --from-bp ", Min, " --to-bp ", Max, " -c")
    # LD <- fread(cmd=CMD)

    # Make snp x phenotype matrix of beta std errors
    ses <- SummaryStats.filtered %>%
        as_tibble() %>%
        select(phenotype, snp, beta_se) %>%
        distinct(phenotype, snp, .keep_all=T) %>%
        spread(phenotype, beta_se) %>%
        drop_na() %>%
        column_to_rownames("snp") %>%
        as.matrix()

    # hyprcoloc does not accept 0 values for ses. Manually Replace 0 values with next smallest value
    ses[ses == 0] <- sort(ses[ses>0])[1]

    traits <- colnames(betas)
    rsid <- rownames(betas)
    if (length(traits)==1){
        print("Only one triat... skipping locus")
        next
    }
    res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)

    HyprcolocResults.list[[i]] <- 
         res[[1]] %>%
         as.data.frame()
}
HyprcolocResults.list <- setNames(HyprcolocResults.list, names(Loci))

HyprcolocResults.df <-
    bind_rows(HyprcolocResults.list, .id="loci_name")

write_tsv(HyprcolocResults.df, FileOut, append=T)
