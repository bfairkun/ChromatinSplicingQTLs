#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : hyprcoloc_scratch
# @created     : Wednesday Aug 04, 2021 15:12:12 CDT
#
# @description : for loop thru gwas loci attempting to colocalize mol QTL
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 'hyprcoloc/LociWiseSummaryStatsInput/ForGWASColoc/GCST004603.txt.gz gwas_summary_stats/StatsForColoc/GCST004603.standardized.txt.gz scratch/test.out.txt.gz ', what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(data.table)
library(tidyverse)
library(hyprcoloc)

FileIn <- args[1]
FileIn_GWAS <- args[2]
FileOut <- args[3]

if (args[4] == '' | is.na(args[4])){
    TraitClassesRestrictions = NA
} else {
    TraitClassesRestrictions <- 
        unlist(strsplit(args[4], ' '))
}

SummaryStats <- fread(FileIn, nrows=Inf) %>%
    separate(snp, into=c("chrom", "pos", "A1_mol", "A2_mol"), sep=":", convert=T, remove=F) %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    mutate(gwas_locus = str_replace(gwas_locus, "(.+?)_(.+?)_.+_(.+?)$", "chr\\1_\\2_\\3"))

Gwas_SummaryStats <- fread(FileIn_GWAS) %>%
    mutate(phenotype = str_replace(loci, "(.+?)_(.+?)_N_N_(.+)$", "\\1_\\2_\\3")) %>%
    select(chrom, pos=start, phenotype, beta, beta_se=SE, A1, A2)


Loci <- SummaryStats$gwas_locus %>% unique()

HyprcolocResults.list <- vector("list", length(Loci))

for (i in seq_along(Loci)){

    # TestLocus <- Loci[1]
    TestLocus <- Loci[i]
    print(paste0("running loci", i, " : ", TestLocus))

    SummaryStats.filtered <- SummaryStats %>%
        filter(gwas_locus == TestLocus) %>%
        mutate(phenotype_class = str_replace(source_file, "QTLs/QTLTools/(.+?)/(.+?)/.+$", "\\1")) %>%
        filter( if (is.na(TraitClassesRestrictions) ) TRUE else phenotype_class %in% TraitClassesRestrictions ) %>%
        unite(phenotype, phenotype_class, phenotype, sep=";") %>%
        select(chrom, pos, snp, phenotype, beta, beta_se, p, A1_mol, A2_mol) %>%
        distinct(.keep_all=T) %>% as_tibble()

    # Make snp x phenotype matrix of betas
    betas <- SummaryStats.filtered %>%
        select(phenotype, snp, chrom, pos, A1_mol, A2_mol, beta) %>%
        distinct(phenotype, chrom, pos, .keep_all=T) %>%
        drop_na() %>%
        spread(phenotype, beta) %>%
        inner_join(
                   Gwas_SummaryStats %>%
                       filter(phenotype == TestLocus) %>%
                       select(beta, A1, A2, chrom, pos) %>%
                       rename(!!TestLocus:=beta),
                   by=c("chrom", "pos")
        ) %>%
        rowwise() %>%
        # Check that molQTL allele set matches gwas, and by permissive if gwas
        # alleles are unknown (listed as "N"). Order of alleles (and sign of
        # betas) doesn't matter for coloc analysis
        filter((setequal(c(A1_mol, A2_mol), c(A1, A2))) | ("N" %in% c(A1, A2))) %>%
        ungroup() %>%
        select(-chrom, -pos, -A1, -A2, -A1_mol, -A2_mol) %>%
        drop_na() %>%
        column_to_rownames("snp") %>%
        as.matrix()

    colnames(betas)


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
        select(phenotype, snp, chrom, pos, A1_mol, A2_mol, beta_se) %>%
        distinct(phenotype, chrom, pos, .keep_all=T) %>%
        drop_na() %>%
        spread(phenotype, beta_se) %>%
        inner_join(
                   Gwas_SummaryStats %>%
                       filter(phenotype == TestLocus) %>%
                       select(beta_se, A1, A2, chrom, pos) %>%
                       rename(!!TestLocus:=beta_se),
                   by=c("chrom", "pos")
        ) %>%
        rowwise() %>%
        # Check that molQTL allele set matches gwas, and by permissive if gwas
        # alleles are unknown (listed as "N"). Order of alleles (and sign of
        # betas) doesn't matter for coloc analysis
        filter((setequal(c(A1_mol, A2_mol), c(A1, A2))) | ("N" %in% c(A1, A2))) %>%
        ungroup() %>%
        select(-chrom, -pos, -A1, -A2, -A1_mol, -A2_mol) %>%
        drop_na() %>%
        column_to_rownames("snp") %>%
        as.matrix()

    # hyprcoloc does not accept 0 values for ses. Manually Replace 0 values with next smallest value
    ses[ses == 0] <- sort(ses[ses>0])[1]

    traits <- colnames(betas)
    rsid <- rownames(betas)
    if (length(traits)==1){
        print("Only one trait... skipping locus")
        next
    }
    res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)

    HyprcolocResults.list[[i]] <- 
         res[[1]] %>%
         as.data.frame()
}
HyprcolocResults.list <- setNames(HyprcolocResults.list, Loci)

HyprcolocResults.df <-
    bind_rows(HyprcolocResults.list, .id="loci_name") %>%
    # skipped gwas loci will still contain a row, but with NAs for hyprcoloc fields
    complete(loci_name = names(HyprcolocResults.list))

write_tsv(HyprcolocResults.df, FileOut, append=T)
