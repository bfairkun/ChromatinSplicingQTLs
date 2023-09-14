#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : CalculatePi1_PlotHeatmap
# @created     : Thursday Aug 10, 2023 15:13:08 CDT
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
library(data.table)
library(qvalue)
library(viridis)


Out_dat <- args[1]
Out_P <- args[2]

theme_set(
  theme_classic() +
  theme(text=element_text(size=16,  family="Helvetica")))


PeaksToTSS <- Sys.glob("../code/Misc/PeaksClosestToTSS/*_assigned.tsv.gz") %>%
  setNames(str_replace(., "../code/Misc/PeaksClosestToTSS/(.+?)_assigned.tsv.gz", "\\1")) %>%
  lapply(read_tsv) %>%
  bind_rows(.id="ChromatinMark") %>%
  mutate(GenePeakPair = paste(gene, peak, sep = ";")) %>%
  distinct(ChromatinMark, peak, gene, .keep_all=T)

PC1.filter <- c("Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1", "Expression.Splicing.Subset_YRI", "chRNA.Expression.Splicing")

dat <- fread("../code/pi1/PairwisePi1Traits.P.all.txt.gz") %>%
  filter((PC1 %in% PC1.filter) & (PC2 %in% PC1.filter))

dat$PC1 %>% unique()

PhenotypeAliases <- read_tsv("../data/Phenotypes_recode_for_Plotting.txt")


# Approximate the null
MaxSampleSizeToCreateANull <- 250
NSamplesToEstimateDistribution <- 10000
NullSimulatedTestStats <- matrix(nrow=MaxSampleSizeToCreateANull, ncol=NSamplesToEstimateDistribution)
rownames(NullSimulatedTestStats) <- paste0("runif_samplesize", 1:MaxSampleSizeToCreateANull)
colnames(NullSimulatedTestStats) <- paste0("Sample_", 1:NSamplesToEstimateDistribution)


for (i in 1:MaxSampleSizeToCreateANull){
  SampleSizeFromUniform <- i
  SampledDat <- matrix(runif(SampleSizeFromUniform*NSamplesToEstimateDistribution), nrow=NSamplesToEstimateDistribution)
  SampleDatNullTestStatistics <- -log10(apply(SampledDat, 1, min))
  NullSimulatedTestStats[i,] <- SampleDatNullTestStatistics
}


ecdf.functions <- apply(NullSimulatedTestStats, 1, ecdf)
ecdf.functions[[1]](1)

# pi0est function breaks sometimes. this wrapper function returns 1 when that happens
CalculatePi1 <- function (dat.in, ...) {
  return(tryCatch(1-pi0est(dat.in$Pvals.For.Pi1, ...)$pi0, error=function(e) 1))
}

# And test it

Recode.PCs = c("ActivatingMark_HostGenePromoter"="Promoter marks*", "H3K36ME3"="H3K36ME3", "chRNA.Expression.Splicing"="chRNA", "MetabolicLabelled.30min"="4sU 30min", "MetabolicLabelled.60min"="4sU 60min", "Expression.Splicing.Subset_YRI"="polyA YRI", "Expression.Splicing"="polyA")

ActivatingMarksToConsider <- c("H3K27AC", "H3K4ME3")

ProcessPi1HeatmapDat <- function(FDR.cutoff=0.1){
  dat.split <- dat %>%
  filter(FDR.x < FDR.cutoff) %>%
  mutate(PC1 = case_when(
    paste(GeneLocus, P1, sep=";") %in% PeaksToTSS$GenePeakPair & PC1 %in% ActivatingMarksToConsider ~ "ActivatingMark_HostGenePromoter",
    P1 %in% PeaksToTSS$peak ~ "ActivatingMark_OtherGenePromoter",
    PC1 %in% c("H3K4ME1", "H3K4ME3", "H3K27AC") ~ "ActivatingMark_Enhancer",
    TRUE ~ PC1
  )) %>%
  mutate(PC2 = case_when(
    paste(GeneLocus, P2, sep=";") %in% PeaksToTSS$GenePeakPair & PC2 %in% ActivatingMarksToConsider ~ "ActivatingMark_HostGenePromoter",
    P2 %in% PeaksToTSS$peak ~ "ActivatingMark_OtherGenePromoter",
    PC2 %in% c("H3K4ME1", "H3K4ME3", "H3K27AC") ~ "ActivatingMark_Enhancer",
    TRUE ~ PC2
  )) %>%
  group_by(PC1, P1, PC2) %>%
  mutate(test.stat.obs = -log10(min(trait.x.p.in.y))) %>%
  ungroup() %>%
  add_count(PC1, P1, PC2) %>%
  filter(n<=250) %>%
  group_by(PC1, PC2) %>%
  rowwise() %>%
  mutate(Pvals.For.Pi1 = 1-ecdf.functions[[n]](test.stat.obs)) %>%
  ungroup() %>%
  dplyr::select(PC1, PC2, Pvals.For.Pi1) %>%
  filter(!PC1==PC2) %>%
  split(paste(.$PC1, .$PC2, sep = ";"))
dat.out <- lapply(dat.split, CalculatePi1, pi0.method="bootstrap") %>%
  unlist() %>%
  data.frame(pi1=.) %>%
  rownames_to_column("PC1_PC2") %>%
  separate(PC1_PC2, into=c("PC1", "PC2"), sep=';') %>%
  filter((PC1%in% c(PC1.filter, "ActivatingMark_HostGenePromoter") & (PC2%in% c(PC1.filter, "ActivatingMark_HostGenePromoter") ))) %>%
  mutate(TextColor = case_when(
    pi1 >= 0.7 ~ "white",
    TRUE ~ "black"
  )) %>%
mutate(PC1 = recode(PC1, !!!Recode.PCs)) %>%
mutate(PC2 = recode(PC2, !!!Recode.PCs)) %>%
mutate(PC1 = factor(PC1, levels=Recode.PCs)) %>%
mutate(PC2 = factor(PC2, levels=Recode.PCs))
P <- ggplot(dat.out, aes(x=PC1, y=PC2, fill=pi1, color=TextColor)) +
  geom_raster() +
  geom_text(aes(label=signif(pi1*100, 2))) +
  scale_color_identity() +
  scale_fill_viridis_c(option="A", direction = -1, limits=c(0,1)) +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Discovery QTL phenotype", y="Phenotype assessed for overlap", caption=paste("*Promoter marks include", paste(ActivatingMarksToConsider, collapse=",")))
return(list(dat=dat.out, P=P))
}


ActivatingMarksToConsider <- c("H3K27AC", "H3K4ME3")

pi.dat <- bind_rows(
    FDR10 = ProcessPi1HeatmapDat(FDR.cutoff=0.1)$dat,
    FDR05 = ProcessPi1HeatmapDat(FDR.cutoff = 0.05)$dat,
    FDR01 = ProcessPi1HeatmapDat(FDR.cutoff = 0.01)$dat,
    .id = "FDR"
)

P <- pi.dat %>%
  mutate(PC1 = recode(PC1, "chRNA"="naRNA", "Promoter marks*"="H3K4Me3|H3K27Ac\natPromoter")) %>%
  mutate(PC2 = recode(PC2, "chRNA"="naRNA", "Promoter marks*"="H3K4Me3|H3K27Ac\nat promoter")) %>%
  ggplot(aes(x=PC1, y=PC2, fill=pi1, color=TextColor)) +
  geom_raster() +
  geom_text(aes(label=signif(pi1*100, 2))) +
  scale_color_identity() +
  scale_fill_viridis_c(option="A", direction = -1, limits=c(0,1)) +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  facet_wrap(~FDR) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Discovery QTL phenotype", y="Phenotype assessed for overlap")

write_tsv(pi.dat, Out_dat)
ggsave(Out_P, P, height=5, width=12)
