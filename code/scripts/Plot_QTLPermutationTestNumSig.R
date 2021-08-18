#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : Plot_QTLPermutationTestNumSig
# @created     : Tuesday Jul 06, 2021 10:37:29 CDT
#
# @description : Plot number gene level hits at varying FDR thresholds
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "QC/NumQTLsPerPhenotype.pdf
                     QC/AllQTLPhenotypes.PermutationTest.bed.gz QTLs/QTLTools/chRNA.IR/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/Expression.Splicing/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/chRNA.Expression.Splicing/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/H3K27AC/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/CTCF/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/H3K4ME3/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/chRNA.Splicing/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/polyA.Splicing/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/MetabolicLabelled.30min/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/MetabolicLabelled.60min/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/Expression.Splicing.Subset_YRI/PermutationPass.FDR_Added.txt.gz QTLs/QTLTools/ProCap/PermutationPass.FDR_Added.txt.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

# PermutationPassFDR.Added.Files <- Sys.glob("QTLs/QTLTools/*/PermutationPass.FDR_Added.txt.gz")
PlotOut_f <- args[1]
AllPhenotypesOutBed_f <- args[2]
PermutationPassFDR.Added.Files <- args[-c(1:2)]

library(tidyverse)
library(scales)

dat <- sapply(PermutationPassFDR.Added.Files, read_delim, simplify=FALSE, delim=' ') %>%
    bind_rows(.id = "id") %>%
    mutate(id = str_replace(id, "QTLs/QTLTools/(.+?)/PermutationPass.FDR_Added.txt.gz", "\\1")) %>%
    select(id,  phe_chr, phe_from, phe_to, phe_id, q, phe_strd, var_id, nom_pval, adj_beta_pval, slope )

Pheno_classes <- dat %>% pull(id) %>% unique()
Cols <- data.frame(id=Pheno_classes, hexcol=hue_pal()(length(Pheno_classes)))
Cols$rgb <- col2rgb(Cols$hexcol) %>% apply(2, paste, collapse=",")

feats.tested <- dat %>%
    filter(!is.na(q)) %>%
    count(id) %>%
    mutate(NewLabel = paste0(id, "\n(", n, " feats tested)"))

FDR_HitsTable <-
    dat %>%
    select(id, q) %>%
    mutate(
        `0.1`=q<0.1,
        `0.05`=q<0.05,
        `0.01`=q<0.01) %>%
    select(-q) %>%
    gather("Threshold", "IsSig", -id) %>%
    group_by(Threshold, id) %>%
    summarize(NumSignificant = sum(IsSig, na.rm=T)) %>%
    ungroup() %>%
    left_join(feats.tested, by="id") %>%
    left_join(Cols, by="id")

PlotOut <- ggplot(FDR_HitsTable, aes(x=Threshold, y=NumSignificant, fill=NewLabel)) +
    geom_col() +
    # scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    facet_wrap(~NewLabel, scales="free_y") +
    scale_fill_manual(values=FDR_HitsTable %>% select(NewLabel, hexcol) %>% deframe()) +
    xlab("FDR cutoff") +
    theme_classic() +
    theme(legend.position="none") +
    theme(strip.text = element_text(size=4.5))

dat %>%
    mutate(name = paste(id, phe_id, var_id, adj_beta_pval, q, sep=";")) %>%
    mutate(thickStart = phe_from, thickEnd = phe_to) %>%
    left_join(Cols, by="id") %>%
    select(phe_chr, phe_from, phe_to, name, q, phe_strd, thickStart, thickEnd, rgb) %>%
    arrange(phe_chr, phe_from, phe_to) %>%
    write_tsv(AllPhenotypesOutBed_f, col_names=F)

ggsave(PlotOut_f, PlotOut, height=6, width=6)
