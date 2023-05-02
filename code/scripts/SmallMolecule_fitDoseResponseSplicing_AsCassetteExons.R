#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : SmallMolecule_fitDoseResponseSplicing_AsCassetteExons
# @created     : Thursday Apr 27, 2023 13:10:52 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "../code/SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz ../output/SmallMoleculeGAGT_CassetteExonclusters.bed SmallMolecule/CassetteExons/TidyPSI.tsv.gz SmallMolecule/FitModels/polyA_GAGTIntrons_asPSI.tsv.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_leaf_in <- args[1]
f_exons_in <- args[2]
f_tidy_out <- args[3]
f_fit_out <- args[4]

library(tidyverse)
library(edgeR)
library(data.table)
library(drc)
library(broom)


leafcutter.counts <- read.table(f_leaf_in, header=T, sep=' ') %>%
    rownames_to_column("rowname") %>%
    as_tibble()

cassette.exons <- read_tsv(f_exons_in, col_names=c("chrom","start", "stop","name", "score", "strand", "thickStart", "thickEnd", "color" )) %>%
    separate(name, into=c("type", "dummy", "clusterNum", "GAGTjunc"), sep="_", convert=T) %>%
    dplyr::select(-dummy, -score, -thickStart, -thickEnd, -color) %>%
    mutate(rowname = str_glue("{chrom}:{start}:{stop}:clu_{clusterNum}_{strand}"))

counts.tidy <- leafcutter.counts %>%
    inner_join(cassette.exons) %>%
    pivot_longer(contains("LCL"),names_to="Sample", values_to="Counts") %>%
    dplyr::select(type, GAGTjunc, Counts, Sample) %>%
    pivot_wider(names_from="type", values_from="Counts") %>%
    mutate(IncludedCounts = (junc.GAGT + junc.UpstreamIntron)/2) %>%
    mutate(PSI = IncludedCounts / (IncludedCounts + junc.skipping)) %>%
    separate(Sample, into=c("treatment", "dose.nM", "cell.type", "libType", "RepNumber"), convert=T, sep="_") %>%
    replace_na(list(dose.nM=0))

counts.tidy %>%
    write_tsv(f_tidy_out)

model.dat.df.FilteredGAGT <- counts.tidy %>%
    dplyr::rename(junc=GAGTjunc) %>%
    filter(libType == "polyA") %>%
    nest(-junc)

chRNADoses <- c(0, 100, 3160)

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

ModelFits.Coefficients.GAGTIntrons %>% write_tsv(f_fit_out)
