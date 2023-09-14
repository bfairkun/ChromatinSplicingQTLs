#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : MolQTL_GWAS_QQ
# @created     : Thursday Aug 10, 2023 10:59:20 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/test.qq.pdf scratch/test.qq.dat.tsv.gz SplicingAnalysis/sQTLs_p_and_u.tsv.gz gwas_summary_stats/MolQTLIntersections/IMSGC2019.bed.gz gwas_summary_stats/MolQTLIntersections_ControlSNPs/IMSGC2019/ALL.txt.gz gwas_summary_stats/sorted_index_summarystat_hg38beds/IMSGC2019.bed.gz IMSGC2019", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

library(data.table)
library(tidyverse)
library(scattermore)

Out_pdf <- args[1]
Out_dat <- args[2]
In_sQTLs <- args[3]
In_MolQTLs <- args[4]
In_OtherControlSNPs <- args[5]
In_TestSNPs <- args[6]
gwas_accession <- args[7]

gwas_name <- read_tsv("config/gwas_table.tsv") %>%
    filter(gwas == gwas_accession) %>% distinct() %>% pull(trait)
print(gwas_name)

sQTLs <- fread(In_sQTLs)

molQTL.gwas.P <- fread(In_MolQTLs, col.names = c("chrom", "varPos", "GWAS.P", "molQTL.name", "molQTL.P", "strand", "molQTL.beta", "molQTL.se", "molQTL.q", "varID", "overlap")) %>%
  separate(molQTL.name, into=c("PhenotypeClass", "MolPhenotypeName"), sep = ";")

gwas.P <- fread(In_TestSNPs, select = 4, col.names = "GWAS.P") %>%
  mutate(PhenotypeClass = "All SNPs") %>%
  mutate(GWAS.P = as.numeric(GWAS.P))

Other.Control.SNPs <- ControlSNPs <- fread(In_OtherControlSNPs, col.names=c("GWAS.P", "PhenotypeClassControl"))

PhenotypeClass.filter <- c("Expression.Splicing", "polyA.Splicing", "H3K27AC")


QQ.gwas.dat <- bind_rows(
  molQTL.gwas.P %>%
      filter(PhenotypeClass %in% PhenotypeClass.filter) %>%
      filter(molQTL.q < 0.1) %>%
      left_join(
        sQTLs %>%
          dplyr::select(MolPhenotypeName=P1, sQTL.type, PhenotypeClass = PC1)
      ) %>%
    mutate(PhenotypeClass = case_when(
      !is.na(sQTL.type) ~ sQTL.type,
      TRUE ~ PhenotypeClass
    )),
  Other.Control.SNPs %>%
    filter(PhenotypeClassControl %in% c("H3K27AC", "Expression.Splicing")) %>%
    mutate(PhenotypeClass = recode(PhenotypeClassControl, "H3K27AC"="H3K27AC peak SNPs", "Expression.Splicing"="genic SNPs")) %>%
    dplyr::select(-PhenotypeClassControl),
  gwas.P
  # gwas.P %>%
  #   sample_n(1E5)
  #   sample_frac(1)
) %>%
mutate(GWAS.P = as.numeric(GWAS.P)) %>%
dplyr::select(SNP_group = PhenotypeClass, GWAS.P) %>%
filter(!SNP_group=="polyA.Splicing") %>%
  # pull(SNP_group) %>% unique()
group_by(SNP_group) %>%
  mutate(MyRank = rank(GWAS.P, ties.method='random')) %>%
  mutate(ExpectedP = MyRank/(max(MyRank) + 1)) %>%
  ungroup() %>%
  dplyr::select(-MyRank) %>%
  mutate(ExpectedP = signif(ExpectedP))

# QQ.gwas.dat %>%
#     # arrange(desc(ExpectedP))
#     mutate(Y=-log10(ExpectedP)) %>%
#     arrange(Y)

Factor_Orders <- count(QQ.gwas.dat, SNP_group) %>%
    arrange(desc(n)) %>%
    pull(SNP_group)

QQ.gwas <- 
  QQ.gwas.dat %>%
  # mutate(SNP_group = factor(SNP_group, levels=Factor_Orders, ordered=T)) %>%
  mutate(SNP_group = factor(SNP_group, levels=c("All SNPs", "H3K27AC peak SNPs", "genic SNPs", "Expression.Splicing", "H3K27AC", "u-sQTL", "p-sQTL"))) %>%
  filter(!SNP_group %in% c("H3K27AC peak SNPs", "Expression.Splicing", "genic SNPs", "H3K27AC")) %>%
  mutate(Y = -log10(as.numeric(GWAS.P))) %>%
  mutate(Y = if_else(Y>20, 20, Y)) %>%
  arrange(SNP_group) %>%
  ggplot(aes(x=-log10(as.numeric(ExpectedP)), y=Y, color=SNP_group)) +
  geom_abline(slope=1, intercept=0) +
  # geom_point() +
  geom_scattermore(pixels=c(2E3, 2E3), pointsize=20, alpha=1) +
  # scale_color_brewer(palette = "Set3") +
  scale_color_manual(values=
                       c("genic SNPs"="#969696", "All SNPs"="#000000", "u-sQTL"="#e31a1c", "p-sQTL"="#1f78b4", "H3K27AC"="#6a3d9a", "H3K27AC peak SNPs"="#cab2d6", "Expression.Splicing"="#ff7f00"),
                     labels=c("genic SNPs"="genic SNPs", "All SNPs"="All SNPs", "u-sQTL"="u-sQTL", "p-sQTL"="p-sQTL", "H3K27AC"="H3K27AC QTL", "H3K27AC peak SNPs"="H3K27AC peak SNPs", "Expression.Splicing"="eQTL")) +
  scale_y_continuous(breaks=c(0,5,10,15,20), labels=c(0,5,10,15,">20")) +
  labs(y=NULL,x=NULL, fill="SNP category") +
  theme_classic() +
  theme(legend.position='none') +
  labs(title=str_wrap(gwas_name, 20))

# QQ.gwas <- 
#   QQ.gwas.dat %>%
#   filter(!SNP_group=="polyA.Splicing") %>%
#   ggplot(aes(x=GWAS.P )) +
#   geom_histogram(bins=100) +
#   facet_wrap(~SNP_group, scales="free")
#   labs(y="-log10(ObservedP)", title="MS GWAS QQ plot", caption="GWAS SNPs sub-sampled to 100K for plotting speed", fill="SNP category") +
#   theme_classic()

ZoomIn <- QQ.gwas +
  coord_cartesian(ylim=c(0,20))
  # ylim(c(0,20))
ggsave("/project2/yangili1/carlos_and_ben_shared/rough_figs/OriginalSubplots/Final_GWAS_MS_QQ_V3.pdf", ZoomIn, height=2, width=2)
# ggsave(Out_pdf, ZoomIn, height=2, width=2)

QQ.gwas.dat %>%
    write_tsv(Out_dat)



