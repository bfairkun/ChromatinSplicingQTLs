#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : FitSmallMoleculeModels
# @created     : Friday Jan 13, 2023 14:26:58 CST
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


library(tidyverse)
library(edgeR)
library(data.table)
library(drc)
library(broom)
library(qvalue)

sample_n_of <- function(data, size, ...) {
  dots <- quos(...)
  
  group_ids <- data %>% 
    group_by(!!! dots) %>% 
    group_indices()
  
  sampled_groups <- sample(unique(group_ids), size)
  
  data %>% 
    filter(group_ids %in% sampled_groups)
}

Samples <- read_tsv("../code/config/SmallMoleculeRNASeq.Samples.tsv")
GeneCounts <- read_tsv("../code/SmallMolecule/featureCounts/Counts.txt", comment="#") %>%
  rename_at(vars(-c(1:6)), ~str_replace(.x, "SmallMolecule/AlignmentsPass2/(.+?)/Aligned.sortedByCoord.out.bam", "\\1"))


#limit to top expressed protein coding genes
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

#Get protein coding features, using same gtf as used in featureCounts so that gene_ids perfectly match
gtf <- read_tsv("../code/ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf", comment="#", n_max=Inf, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )) %>%
    filter(feature=="gene") %>%
    dplyr::select(attribute)
gtf$gene_id <- unlist(lapply(gtf$attribute, extract_attributes, "gene_id"))
gtf$gene_type <- unlist(lapply(gtf$attribute, extract_attributes, "gene_type"))
ProteinCodingGenes <- filter(gtf, gene_type=="protein_coding") %>% pull(gene_id)

GenesToInclude <- GeneCounts %>%
  dplyr::select(Geneid, contains("DMSO_NA_LCL_polyA")) %>%
  filter(Geneid %in% ProteinCodingGenes) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors() %>%
  cpm(prior.count=0.1, log=F) %>%
  log2() %>%
  apply(1, mean) %>%
  as.data.frame() %>%
  filter(`.` > 1) %>%
  rownames_to_column("Geneid")

GeneCounts.filtered <- GeneCounts %>%
  filter(Geneid %in% GenesToInclude$Geneid)

Count.table.StandardNormFactors <- GeneCounts.filtered %>%
  dplyr::select(-c(2:6)) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors()


CPM.StandardNormFactors <- Count.table.StandardNormFactors %>%
  cpm(prior.count=0.1, log=F) 


CPM.StandardNormFactors.tidy <- CPM.StandardNormFactors %>%
  as.data.frame() %>%
  rownames_to_column("Geneid") %>%
  gather("Sample", "CPM", -Geneid) %>%
  separate(Sample, into=c("treatment", "dose.nM", "Cell.type", "LibraryType", "rep"), sep="_", convert=T) %>%
  replace_na(list(dose.nM=0)) %>%
  mutate(ensembl_gene_id = str_replace(Geneid, "^(.+?)\\..+$", "\\1"))


### TIDY SPLICING DATA
leafcutter.counts <- read.table("../code/SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz", header=T, sep=' ') %>% as.matrix()

ClusterMax.mat <- leafcutter.counts %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
    group_by(cluster) %>%
    mutate(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    dplyr::select(junc, everything(), -cluster) %>%
    column_to_rownames("junc") %>%
    as.matrix()

PSI.df <- (leafcutter.counts / as.numeric(ClusterMax.mat) * 100) %>%
  signif() %>%
  as.data.frame() %>%
  rownames_to_column("Leafcutter.ID") %>%
  mutate(Intron = str_replace(Leafcutter.ID, "(.+?):(.+?):(.+?):clu_.+?_([+-])$", "\\1:\\2:\\3:\\4"))

Intron.Donors <- fread("../code/SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz", col.names = c("Intron", "DonorSeq", "DonorScore")) %>%
  mutate(Intron = str_replace(Intron, "(.+?)_(.+?)_(.+?)_(.+?)::.+?$", "\\1:\\2:\\3:\\4"))

PSI.tidy <- PSI.df %>%
  left_join(Intron.Donors) %>%
  gather("Sample", "PSI",contains("_")) %>%
  inner_join(
    leafcutter.counts %>%
      as.data.frame() %>%
      rownames_to_column("Leafcutter.ID") %>%
      mutate(Intron = str_replace(Leafcutter.ID, "(.+?):(.+?):(.+?):clu_.+?_([+-])$", "\\1:\\2:\\3:\\4")) %>%
      gather("Sample", "Counts",contains("_"))
  ) %>%
  separate(Sample, into=c("treatment", "dose.nM", "Cell.type", "LibraryType", "rep"), sep="_", convert=T) %>%
  replace_na(list(dose.nM=0))

# Spearman correlation of dose and gene
Spearman.cors.CPM <- CPM.StandardNormFactors.tidy %>%
  nest(-Geneid) %>% 
  mutate(cor=map(data,~cor.test(.x$dose.nM, .x$CPM, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied, .drop = T)

# Spearman of dose and PSI
chRNADoses <- c(0, 100, 3160)


model.dat.df.AllGAGT <- PSI.tidy %>%
  filter(LibraryType == "polyA") %>%
  drop_na() %>%
  add_count(Intron) %>%
  filter(n==11) %>%
  mutate(Is.GA.GT = substr(DonorSeq, 3,4)) %>%
  filter(Is.GA.GT == "GA") %>%
  nest(-Intron) %>%
  mutate(cor=map(data,~cor.test(.x$dose.nM, .x$PSI, method = "sp", alternative="greater"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  dplyr::select(Intron:data, spearman=estimate, spearman.p = p.value) %>%
  mutate(q = qvalue(spearman.p)$qvalues)


#pre-filter data for model fitting. just do GAGT introns with positive dose response cor (q<0.01)
model.dat.df.FilteredGAGT <- model.dat.df.AllGAGT %>%
  filter(q<0.01) %>%
  dplyr::select(junc=Intron, everything())

# Fit model for splicing data
Results <- list()
for(i in 1:nrow(model.dat.df.FilteredGAGT)) {
# for(i in 1:10) {
  tryCatch(
    expr = {
      junc <- model.dat.df.FilteredGAGT$junc[i]
      data <- model.dat.df.FilteredGAGT$data[i] %>% as.data.frame()
      fit <- drm(formula = PSI ~ dose.nM,
                  data = data,
                  fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50")),
                  lowerl=c(NA,0,NA,NA),
                  upperl=c(NA,NA,100,NA),
                  robust = "mean"
                  )
      df.out <- 
        bind_rows(
          coef(summary(fit)) %>%
            as.data.frame() %>%
            rownames_to_column("param") %>%
            dplyr::select(param, Estimate, SE=`Std. Error`),
          predict(fit, data.frame(dose.nM=chRNADoses), se.fit = T) %>%
            as.data.frame() %>%
            dplyr::rename("Estimate"="Prediction") %>%
            mutate(param=paste("Pred",chRNADoses, sep="_"))
        )
      Results[[junc]] <- df.out
      message("Successfully fitted model.")
    },
    error=function(e){
      if (i < 100){
        cat("ERROR :",conditionMessage(e), junc, "\n")
      }
      })
}

ModelFits.Coefficients.GAGTIntrons <- bind_rows(Results, .id="junc")


### TIDY SPLICING DATA JUST COMPETING 3ss
leafcutter.counts <- read.table("../code/SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz", header=T, sep=' ') %>% as.matrix()

ClusterMax.mat <- leafcutter.counts %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
    group_by(cluster) %>%
    mutate(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    dplyr::select(junc, everything(), -cluster) %>%
    column_to_rownames("junc") %>%
    as.matrix()

PSI.df <- (leafcutter.counts / as.numeric(ClusterMax.mat) * 100) %>%
  signif() %>%
  as.data.frame() %>%
  rownames_to_column("Leafcutter.ID") %>%
  mutate(Intron = str_replace(Leafcutter.ID, "(.+?):(.+?):(.+?):clu_.+?_([+-])$", "\\1:\\2:\\3:\\4"))

Intron.Donors <- fread("../code/SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz", col.names = c("Intron", "DonorSeq", "DonorScore")) %>%
  mutate(Intron = str_replace(Intron, "(.+?)_(.+?)_(.+?)_(.+?)::.+?$", "\\1:\\2:\\3:\\4"))

PSI.tidy <- PSI.df %>%
  left_join(Intron.Donors) %>%
  gather("Sample", "PSI",contains("_")) %>%
  inner_join(
    leafcutter.counts %>%
      as.data.frame() %>%
      rownames_to_column("Leafcutter.ID") %>%
      mutate(Intron = str_replace(Leafcutter.ID, "(.+?):(.+?):(.+?):clu_.+?_([+-])$", "\\1:\\2:\\3:\\4")) %>%
      gather("Sample", "Counts",contains("_"))
  ) %>%
  separate(Sample, into=c("treatment", "dose.nM", "Cell.type", "LibraryType", "rep"), sep="_", convert=T) %>%
  replace_na(list(dose.nM=0))

# Spearman correlation of dose and gene
Spearman.cors.CPM <- CPM.StandardNormFactors.tidy %>%
  nest(-Geneid) %>% 
  mutate(cor=map(data,~cor.test(.x$dose.nM, .x$CPM, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied, .drop = T)

# Spearman of dose and PSI
chRNADoses <- c(0, 100, 3160)


model.dat.df.AllGAGT <- PSI.tidy %>%
  filter(LibraryType == "polyA") %>%
  drop_na() %>%
  add_count(Intron) %>%
  filter(n==11) %>%
  mutate(Is.GA.GT = substr(DonorSeq, 3,4)) %>%
  filter(Is.GA.GT == "GA") %>%
  nest(-Intron) %>%
  mutate(cor=map(data,~cor.test(.x$dose.nM, .x$PSI, method = "sp", alternative="greater"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  dplyr::select(Intron:data, spearman=estimate, spearman.p = p.value) %>%
  mutate(q = qvalue(spearman.p)$qvalues)


#pre-filter data for model fitting. just do GAGT introns with positive dose response cor (q<0.01)
model.dat.df.FilteredGAGT <- model.dat.df.AllGAGT %>%
  filter(q<0.01) %>%
  dplyr::select(junc=Intron, everything())

# Fit model for splicing data
Results <- list()
for(i in 1:nrow(model.dat.df.FilteredGAGT)) {
# for(i in 1:10) {
  tryCatch(
    expr = {
      junc <- model.dat.df.FilteredGAGT$junc[i]
      data <- model.dat.df.FilteredGAGT$data[i] %>% as.data.frame()
      fit <- drm(formula = PSI ~ dose.nM,
                  data = data,
                  fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50")),
                  lowerl=c(NA,0,NA,NA),
                  upperl=c(NA,NA,100,NA),
                  robust = "mean"
                  )
      df.out <- 
        bind_rows(
          coef(summary(fit)) %>%
            as.data.frame() %>%
            rownames_to_column("param") %>%
            dplyr::select(param, Estimate, SE=`Std. Error`),
          predict(fit, data.frame(dose.nM=chRNADoses), se.fit = T) %>%
            as.data.frame() %>%
            dplyr::rename("Estimate"="Prediction") %>%
            mutate(param=paste("Pred",chRNADoses, sep="_"))
        )
      Results[[junc]] <- df.out
      message("Successfully fitted model.")
    },
    error=function(e){
      if (i < 100){
        cat("ERROR :",conditionMessage(e), junc, "\n")
      }
      })
}

ModelFits.Coefficients.GAGTIntrons <- bind_rows(Results, .id="junc")



#Model genes

model.dat.df.CPM.StandardNormFactors <- CPM.StandardNormFactors.tidy %>%
  filter(LibraryType == "polyA") %>%
  mutate(Log2CPM = log2(CPM)) %>%
  nest(data=-Geneid)


Results <- list()
for(i in 1:nrow(model.dat.df.CPM.StandardNormFactors)) {
# for(i in 1:100) {
  tryCatch(
    expr = {
      Geneid <- model.dat.df.CPM.StandardNormFactors$Geneid[i]
      data <- model.dat.df.CPM.StandardNormFactors$data[i] %>% as.data.frame()
      fit <- drm(formula = Log2CPM ~ dose.nM,
                  data = data,
                  fct = LL.4(names=c("Steepness", "LowerLimit", "UpperLimit", "ED50")),
                  robust = "mean"
                  )
      df.out <- 
        bind_rows(
          coef(summary(fit)) %>%
            as.data.frame() %>%
            rownames_to_column("param") %>%
            dplyr::select(param, Estimate, SE=`Std. Error`),
          predict(fit, data.frame(dose.nM=chRNADoses), se.fit = T) %>%
            as.data.frame() %>%
            dplyr::rename("Estimate"="Prediction") %>%
            mutate(param=paste("Pred",chRNADoses, sep="_"))
        )
      Results[[Geneid]] <- df.out
      message("Successfully fitted model.")
    },
    error=function(e){
      if (i < 100){
        cat("ERROR :",conditionMessage(e), Geneid, "\n")
      }
      })
}


ModelFits.Coefficients.Genes <- bind_rows(Results, .id="Geneid")

write_tsv(ModelFits.Coefficients.Genes, "SmallMolecule/FitModels/polyA_genes.tsv.gz")
write_tsv(ModelFits.Coefficients.GAGTIntrons, "SmallMolecule/FitModels/polyA_GAGTIntrons.tsv.gz")
