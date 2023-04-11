library(dplyr)
library(tidyverse)
library(data.table)
library(susieR)
set.seed(0)
args = commandArgs(trailingOnly=TRUE)
nominal = args[1]
perm = args[2]
genotype = args[3]
Subset = args[4]
output = args[5]


nominal <- read.table(nominal, sep=' ', header=TRUE)
perm <- read.table(perm, sep=' ', header=TRUE)

run_susie <- function(nominal.gene, gene, header, genotype) {
  
    # Get variants positions
    chrom <- nominal.gene$phe_chr[1] %>% as.character() %>% substr(4, 5)
    var_from <- nominal.gene$var_from[1] %>% as.character()
    var_to <- tail(nominal.gene$var_from, n=1) %>% as.character()
    position = paste0(chrom, ':', var_from, '-', var_to)
    
    # Get genotypes of region
    
    tabix_cmd <- paste("tabix", genotype, position)
        
    GT <- fread(cmd = str_glue(tabix_cmd),
                  data.table = F, header = F)
    
    colnames(GT) <- header
    
    X <- GT %>% filter(ID %in% nominal.gene$var_id)
    XCor <- cor(t(X[,colnames(X[6:dim(GT)[2]])]))
    XCor[is.na(XCor)] = 0
    
    Z <- (nominal.gene$slope / nominal.gene$slope_se)
    
    n <- dim(GT)[2] - 6
    
    tryCatch({
        susie_out <- susie_rss(Z, XCor, L=3, n=n)
    } , error = function(e){
        susie_out <- susie_rss(Z, XCor, L=3, n=n, estimate_prior_variance = FALSE)
    }
    )
    
    return(susie_out)
}

head_cmd <- paste("zcat", genotype, "| head")

header <- fread(cmd = str_glue(head_cmd),
                  data.table = F, header = T) %>% colnames()

genes <- perm[perm$q <= 1e-1, ]$phe_id

cs_names <- c()
cs_vars <- c()
cs_gene <- c()
cs_pip <- c()
cs_beta <- c()
cs_beta_se <- c()
cs_z <- c()

for (gene in genes) {

    print(gene)
    # Get the nominal pass rows that match the gene
    nominal.gene <- nominal %>% filter(phe_id == gene)
    
    susie_out <- run_susie(nominal.gene, gene, header, genotype)
    
    susie_rdata_file <- paste0('FineMapping/susie_runs_', Subset, '/', gene, '.RData')
    
    save(susie_out, file = susie_rdata_file)
    
    credible_sets <- names(susie_out$sets$cs) 
    
    for (cs in credible_sets){
        cs_idx <- susie_out$sets$cs[[cs]]
        cs_names <- c(cs_names, rep(cs, length(cs_idx)))
        vars <- nominal.gene[cs_idx,] %>% select('var_id') 
        vars <- vars[[1]] %>% as.vector()
        cs_vars <- c(cs_vars, vars)
        cs_gene <- c(cs_gene, rep(gene, length(cs_idx)))
        cs_pip <- c(cs_pip, susie_out$pip[cs_idx])
        
        beta <- nominal.gene[cs_idx,] %>% select('slope') 
        beta <- beta[[1]] %>% as.vector()
        
        beta_se <- nominal.gene[cs_idx,] %>% select('slope_se') 
        beta_se <- beta_se[[1]] %>% as.vector()
        
        z <- beta/beta_se
        
        cs_beta <- c(cs_beta, beta)
        cs_beta_se <- c(cs_beta_se, beta_se)
        cs_z <- c(cs_z, z)
    }

}

print('Finished running SuSiE. Creating output table...')
out_df <- data.frame(cs_gene, cs_vars, cs_names, cs_pip, cs_beta, cs_beta_se, cs_z)
print("If this message doesn't show up, the error is creating the matrix")
out_df %>% write_delim(output, delim='\t')
