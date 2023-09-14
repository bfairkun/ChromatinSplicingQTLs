#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Plot_eQTL_QQ
# @created     : Thursday Aug 10, 2023 16:45:33 CDT
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


library(data.table)
library(tidyverse)
library(scattermore)

Out_dat <- args[1]
Out_P <- args[2]

# PeaksToTSS <- Sys.glob("../code/Misc/PeaksClosestToTSS/*_assigned.tsv.gz") %>%
#   setNames(str_replace(., "../code/Misc/PeaksClosestToTSS/(.+?)_assigned.tsv.gz", "\\1")) %>%
#   lapply(read_tsv) %>%
#   bind_rows(.id="ChromatinMark") %>%
#   mutate(GenePeakPair = paste(gene, peak, sep = ";")) %>%
#   distinct(ChromatinMark, peak, gene, .keep_all=T)

sQTLs <- fread("SplicingAnalysis/sQTLs_p_and_u.tsv.gz")

PC1.filter <- c("polyA.Splicing", "H3K27AC")
PC2.filter <- c("Expression.Splicing" )
molQTLs <- fread("pi1/PairwisePi1Traits.P.all.txt.gz") %>%
    filter(PC1 %in% PC1.filter & PC2 %in% PC2.filter)

Test.SNPs <- fread("QTLs/QTLTools/Expression.Splicing/NominalPassForColoc.RandomSamplePvals.txt.gz")

QQ.dat <- bind_rows(
  molQTLs %>%
        left_join(
                  sQTLs %>%
                      dplyr::select(PC1, P1, sQTL.type, trait.x.p.in.y)
        ) %>%
      filter(PC1 == "H3K27AC" | !is.na(sQTL.type)) %>%
      mutate(SNP_group = case_when(
                                   PC1 == "H3K27AC" ~ "hQTL",
                                   !is.na(sQTL.type) ~ sQTL.type
                                   )) %>%
      dplyr::select(P=trait.x.p.in.y, SNP_group),
  Test.SNPs %>%
      mutate(SNP_group = "Test SNPs") %>%
      dplyr::select(SNP_group, P=V1)) %>%
  group_by(SNP_group) %>%
  mutate(MyRank = rank(P, ties.method='random')) %>%
  mutate(ExpectedP = MyRank/(max(MyRank) + 1)) %>%
  ungroup() %>%
  dplyr::select(-MyRank) %>%
  mutate(ExpectedP = signif(ExpectedP))

Factor_Orders <- count(QQ.dat, SNP_group) %>%
    arrange(desc(n)) %>%
    pull(SNP_group)

QQ.eQTL<- 
  QQ.dat %>%
  mutate(Y = -log10(as.numeric(P))) %>%
  # mutate(Y = if_else(Y>20, 20, Y)) %>%
  arrange(SNP_group) %>%
  ggplot(aes(x=-log10(as.numeric(ExpectedP)), y=Y, color=SNP_group)) +
  geom_abline(slope=1, intercept=0) +
  # geom_point() +
  geom_scattermore(pixels=c(2E3, 2E3), pointsize=20, alpha=1) +
  # scale_color_brewer(palette = "Set3") +
  scale_color_manual(values=
                       c("Test SNPs"="#000000", "u-sQTL"="#e31a1c", "p-sQTL"="#1f78b4", "hQTL"="#6a3d9a")) +
  labs(y=NULL,x=NULL, fill="SNP category") +
  theme_classic() +
  theme(legend.position='none')

ggsave(Out_P, QQ.eQTL, height=2, width=2)

QQ.dat %>%
    write_tsv(Out_dat)
