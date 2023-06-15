library(dplyr)
library(tidyverse)
library(data.table)
library(bedr)
library(readr)
library(susieR)
options(readr.num_columns = 0)
set.seed(0)



get.gene.coords <- function(perm, gene){
    chrom <- perm %>% filter(phe_id==gene) %>% pull(phe_chr) %>% as.character()
    start <- perm %>% filter(phe_id==gene) %>% pull(phe_from) %>% as.integer()
    end <- perm %>% filter(phe_id==gene) %>% pull(phe_to) %>% as.integer()
    
    start <- max(0, start - 500000) %>% as.character()
    end <- (end + 500000) %>% as.character()
    
    loc <- paste(start, end, sep='-')
    coords <- paste(chrom, loc, sep=':')
    return(coords)
}

get.nominal.slice <- function(nominal.file, gene, coords){

    nominal.slice <- tabix(coords, nominal.file, check.chr=FALSE, verbose=FALSE) %>%
        filter(phe_id == gene)
    return (nominal.slice)
}

get.genotype <- function(genotype, coords){
    coords_gsub <- gsub('*chr', '', coords)
    head_cmd <- paste("zcat", genotype, "| head")
    header <- fread(cmd = str_glue(head_cmd),
                      data.table = F, header = T) %>% colnames()
    samples <- header[6:length(header)]
    tabix_cmd <- paste("tabix", genotype, coords_gsub)
    GT <- fread(cmd = str_glue(tabix_cmd),
                  data.table = F, header = F)
    colnames(GT) <- header
    GT <- GT[rowSums(GT[,samples]) > 0,]
    return(GT)
    
}

filter.GT <- function(nominal.gene, GT){
    
    var_id_ng <- nominal.gene %>% pull(var_id) 
    var_id_gt <- GT %>% pull(ID) 
    shared_vars <- var_id_ng[var_id_ng %in% var_id_gt]
    GT_filtered <- GT %>% filter(ID %in% shared_vars)
    row.names(GT_filtered) <- GT_filtered %>% pull(ID)
    row.names(nominal.gene) <- nominal.gene %>% pull(var_id)
    nominal.gene_filtered <- nominal.gene[row.names(GT_filtered),]
    
    return (list(GT=GT_filtered, nominal.gene=nominal.gene_filtered))
}


#filter.GT <- function(nominal.gene, GT){
#    GT <- GT %>% filter(ID %in% nominal.gene$var_id)
#    nominal.gene <- nominal.gene %>% filter(var_id %in% GT$ID)
#    return (list(GT=GT, nominal.gene=nominal.gene))
#}


genotype.cor <- function(GT){
    XCor <- cor(t(GT[,colnames(GT[6:dim(GT)[2]])]))
    return(XCor)
}

prepare.susie.input <- function(nominal.gene, GT){
    slope <- nominal.gene$slope %>% as.numeric()
    slope_se <- nominal.gene$slope_se %>% as.numeric()
    nom.pval <- nominal.gene$nom_pval %>% as.numeric()
    Z <- (slope / slope_se)
    n <- dim(GT)[2] - 5
    XCor <- genotype.cor(GT)
    return(list(Z=Z, XCor=XCor, n=n, slope=slope, slope_se=slope_se, nom.pval=nom.pval))
}


run.susie <- function(susie.input){
    Z <- susie.input$Z
    XCor <- susie.input$XCor
    n <- susie.input$n
    
    tryCatch(
        {
          susie.out <- susie_rss(Z, XCor, L=10, n=n,
                             estimate_prior_variance = TRUE,
                             scaled_prior_variance = 1,
                             min_abs_corr = 0.1,
                             verbose = FALSE)
        }, error = function(e){
          print('No prior variance')
          susie.out <- susie_rss(Z, XCor, L=10, n=n,
                             estimate_prior_variance = FALSE,
                             scaled_prior_variance = 1,
                             min_abs_corr = 0.1,
                             verbose = FALSE)
          return(susie.out)
        })
    
    
    
    return(susie.out)
}

finemap.susie <- function(perm, nominal.file, genotype.file, gene){
    coords <- get.gene.coords(perm, gene)
    nominal.gene <- get.nominal.slice(nominal.file, gene, coords)
    GT <- get.genotype(genotype.file, coords)
    GT.filtered <- filter.GT(nominal.gene, GT)
    nominal.gene <- GT.filtered$nominal.gene
    GT <- GT.filtered$GT
    susie.input <- prepare.susie.input(nominal.gene, GT)
    susie.out <- run.susie(susie.input)
    return(list(nominal.gene=nominal.gene, susie.input=susie.input, susie.out=susie.out))
}

args = commandArgs(trailingOnly=TRUE)
nominal = args[1]
perm = args[2]
genotype = args[3]
Subset = args[4]
output = args[5]

perm <- read.table(perm, sep=' ', header=TRUE)


genes <- perm[perm$q <= 1e-1, ]$phe_id

cs_names <- c()
cs_vars <- c()
cs_gene <- c()
cs_pip <- c()
cs_beta <- c()
cs_beta_se <- c()
cs_z <- c()
cs_pval <- c()

for (gene in genes) {

    print(gene)
    
    sus <- finemap.susie(perm, nominal,
              genotype, gene
             )
    
    print(sus$susie.out$sets$cs)
    
    susie_rdata_file <- paste0('FineMapping/susie_runs_', Subset, '/', gene, '.RData')
    
    save(sus, file = susie_rdata_file)
    
    credible_sets <- names(sus$susie.out$sets$cs) 
    
    for (cs in credible_sets){
        cs_idx <- sus$susie.out$sets$cs[[cs]]
        cs_names <- c(cs_names, rep(cs, length(cs_idx)))
        vars <- sus$nominal.gene[cs_idx,] %>% select('var_id') 
        vars <- vars[[1]] %>% as.vector()
        cs_vars <- c(cs_vars, vars)
        cs_gene <- c(cs_gene, rep(gene, length(cs_idx)))
        cs_pip <- c(cs_pip, sus$susie.out$pip[cs_idx])
        
        beta <- sus$susie.input$slope[cs_idx]#nominal.gene[cs_idx,] %>% select('slope') 
        
        beta_se <- sus$susie.input$slope_se[cs_idx]#beta_se[[1]] %>% as.vector()
        
        z <- sus$susie.input$Z[cs_idx]#beta/beta_se
        
        pval <- sus$susie.input$nom.pval[cs_idx]
        
        cs_beta <- c(cs_beta, beta)
        cs_beta_se <- c(cs_beta_se, beta_se)
        cs_z <- c(cs_z, z)
        cs_pval <- c(cs_pval, pval)
    }
    
    
    #cs_vars %>% length() %>% print()

}

cs_gene %>% length %>% print()
cs_vars %>% length %>% print()
cs_names %>% length %>% print()
cs_pip %>% length %>% print()
cs_beta %>% length %>% print()
cs_beta_se %>% length %>% print()
cs_z %>% length %>% print()
cs_pval %>% length %>% print()

print('Finished running SuSiE. Creating output table...')
out_df <- data.frame(cs_gene, cs_vars, cs_names, cs_pip, cs_beta, cs_beta_se, cs_z, cs_pval)
print("If this message doesn't show up, the error is creating the matrix")
out_df %>% write_delim(output, delim='\t')

