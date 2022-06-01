#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : hyprcoloc_scratch
# @created     : Wednesday Aug 04, 2021 15:12:12 CDT
#
# @description : for loop thru loci attempting to colocalize mol QTL
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
# Usage: <Filelist of loci-wise summary stats files> <Fileout_clusters> <Fileout_finemappingScores> '[Optional list of space delimited phenotype classes to include, surrounded by single quotes]' 
if(interactive()){
    # args <- scan(text=
    #              "hyprcoloc/Results/ForColoc/Chunks/4.list.txt scratch/testshyprcoloc.txt.gz scratch/testhyprcoloc.snpscores.txt.gz 0.01 'MetabolicLabelled.30min MetabolicLabelled.60min Expression.Splicing.Subset_YRI'", what='character')
    args <- scan(text=
                 "hyprcoloc/Results/ForColoc/Chunks/4.list.txt scratch/testshyprcoloc.txt.gz scratch/testhyprcoloc.snpscores.txt.gz 0.01 ''", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

# SummaryStatsInList FileOut

library(data.table)
library(tidyverse)
library(hyprcoloc)

FilelistIn <- args[1]
FileOut <- args[2]
FinemappingOut <- args[3]
Threshold <- as.numeric(args[4])

PermutationResultFiles <- Sys.glob("QTLs/QTLTools/*/PermutationPassForColoc.txt.gz") %>%
    setNames(str_replace(., "QTLs/QTLTools/(.+?)/PermutationPassForColoc.txt.gz", "\\1"))

if (args[5] == '' | is.na(args[5])){
    TraitClassesRestrictions = NA
} else {
    TraitClassesRestrictions <- unlist(strsplit(args[5], ' '))
    PermutationResultFiles <- PermutationResultFiles[TraitClassesRestrictions]
}

# Get list of phenotypes to include for coloc attempts
PhenotypesToColoc <- 
    lapply(PermutationResultFiles, read_delim, delim=' ') %>%
    bind_rows(.id = "phenotype_class") %>%
    select(phenotype_class, phe_id, adj_emp_pval) %>%
    filter(adj_emp_pval < Threshold) %>%
    mutate(phe_id = str_replace(phe_id, "(.+):(.+)", "\\1;\\2")) %>%
    separate(phe_id, into=c("phenotype", "Locus"), extra="merge", sep=";") %>%
    unite(phenotype, phenotype_class, phenotype, sep=";")

Loci<- read_tsv(FilelistIn, col_names=c("f")) %>%
    filter(file.exists(f)) %>%
    mutate(Loci_names = str_replace(f, "hyprcoloc/LociWiseSummaryStatsInput/ForColoc/(.+?)\\.txt\\.gz", "\\1")) %>%
    select(Loci_names, f) %>%
    deframe()


# HyprcolocResults.list <- vector("list", length(Loci))

for (i in seq_along(Loci)){

    # TestLocus <- Loci[875]
    TestLocus <- Loci[i]
    print(paste0("running loci", i, " : ", names(TestLocus)))


    SummaryStats.filtered <- fread(TestLocus) %>%
        filter(gwas_locus == names(TestLocus)) %>%
        mutate(phenotype_class = str_replace(source_file, "QTLs/QTLTools/(.+?)/(.+?)/.+$", "\\1")) %>%
        # filter( if (is.na(TraitClassesRestrictions) ) TRUE else phenotype_class %in% TraitClassesRestrictions )  %>%
        unite(phenotype, phenotype_class, phenotype, sep=";") %>%
        inner_join(
                   PhenotypesToColoc %>%
                       filter(Locus == names(TestLocus)),
                   by = "phenotype"
        ) %>%
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
    if (length(colnames(betas)) >= 2) {
        print(i)
    }

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
    # changed == to <= to account for when we are not colocalizing some traits
    if (length(traits)<=1){
        print("Only one trait... skipping locus")
        next
    }
    res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, snpscores=T)
    print(dim(res[[1]]))
    try(
        res[[2]] %>%
            setNames(1:length(res[[2]])) %>%
            bind_rows() %>%
            mutate(snps=rsid) %>%
            gather(key="ColocalizedCluster", value="FinemappingPr", -snps) %>%
            group_by(ColocalizedCluster) %>%
            arrange(desc(FinemappingPr)) %>%
            mutate(cumPr=cumsum(FinemappingPr)) %>%
            mutate(lagcumPr=lag(cumPr)) %>%
            filter((lagcumPr < 0.95) | is.na(lagcumPr)) %>%
            ungroup() %>%
            select(snps, ColocalizedCluster, FinemappingPr) %>%
            mutate(Locus=names(TestLocus)) %>%
            write_tsv(FinemappingOut, append=T)

    )

    # ResToWrite <- res[[1]] %>%
    #      as.data.frame() %>%
    #      mutate(Loci = names(TestLocus)) %>%
    #      select(Loci, everything())
     # ResToWrite %>%
     #     separate_rows(ColocalizedTraits, sep = ', ') %>%
     #     mutate(Trait = case_when(
     #        ColocalizedTraits == "None" ~ DroppedTrait,
     #        TRUE ~ ColocalizedTraits
     #     ))
     #     separateres[[1]] %>%

     res[[1]] %>%
         as.data.frame() %>%
         mutate(dropped_trait = as.character(dropped_trait)) %>%
         mutate(Traits = case_when(
                                   traits == "none" ~ dropped_trait,
                                   TRUE ~ traits
                                   )) %>%
         select(-traits, -dropped_trait) %>%
         separate_rows(Traits, sep=", ") %>%
         right_join(
                    data.frame(Traits = as.character(colnames(betas))),
                    by = "Traits"
         ) %>%
         mutate(Loci = names(TestLocus)) %>%
         select(Loci, everything()) %>%
         write_tsv(FileOut, append=T)

}
