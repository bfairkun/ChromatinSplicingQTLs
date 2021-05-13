#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : TempPermuteAndPCA
# @created     : Tuesday May 11, 2021 13:50:56 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "InteractiveMode Test Args", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


# FileIn <- "scratch/GSE75220_4su_30_qqnorm.txt.gz"
FileIn <- "scratch/GSE75220_4su_60_qqnorm.txt.gz"

mat <- read_delim(, delim=" ") %>%
    column_to_rownames("gene") %>% as.matrix()

mat.permuted <- matrix(nrow=nrow(mat), ncol=ncol(mat))
for (i in 1:nrow(mat)){
    mat.permuted[i,] <- sample(mat[i,], size=ncol(mat))
}

pca <- prcomp(t(mat)) %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>%
    rownames_to_column("PC")
pca.permuted <- prcomp(t(mat.permuted)) %>% summary() %>% extract2("importance") %>% t() %>% as.data.frame() %>% rownames_to_column("PC")

merged <- bind_rows(list(pca=pca, pca.permuted=pca.permuted), .id="mat") %>%
    mutate(PC=as.numeric(str_replace(PC, "PC", "")))

#GetNumPCs
merged %>%
    select(PC, Prop=`Proportion of Variance`, mat) %>%
    spread(key="mat", value="Prop") %>%
    filter(pca > pca.permuted) %>% pull(PC) %>% max()


merged %>%
    ggplot(aes(x=PC, y=`Proportion of Variance`, color=mat)) +
    geom_line() +
    theme_bw()
