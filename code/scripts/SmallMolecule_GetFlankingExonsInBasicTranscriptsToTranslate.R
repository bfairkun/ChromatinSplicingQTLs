#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : SmallMolecule_GetFlankingExonsInBasicTranscriptsToTranslate
# @created     : Wednesday May 03, 2023 10:37:56 CDT
#
# @description : 
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "scratch/test.tsv.gz scratch/test.tsv.long.gz", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

f_out <- args[1]
f_long_out <- args[2]

library(tidyverse)

salmon = "SmallMolecule/salmon.DMSO.merged.txt"
OverlappingExons = "SmallMolecule/CassetteExons/FlankingExons.tsv.gz"
CassetteExons = "../output/SmallMoleculeGAGT_CassetteExonclusters.bed"

IntersectingExons <- read_tsv(OverlappingExons, col_names=c("SkipJuncChrom", "SkipJuncStart", "SkipJuncStop", "SkipJuncName", "SkipJuncScore", "Strand", "dummy_thickStart", "dummy_thickEnd", "dummy_color", "ExonChrom","ExonStart", "ExonStop", "gene_transcript", "dummy", "dummy_2", "overlappingBases")) %>%
  dplyr::select(-contains("dummy")) %>%
  separate(gene_transcript, into=c("gene", "transcript"), sep="\\_") %>%
  mutate(SkipJuncStart = SkipJuncStart+1, SkipJuncStop = SkipJuncStop -1)

TPM.transcripts.DMSO <- read_tsv(salmon) %>%
  column_to_rownames("Name") %>%
  apply(1, median) %>%
  as.data.frame() %>%
  rownames_to_column("transcript") %>%
  dplyr::rename("Median.TPM"=".")

cassette.exons <- read_tsv(CassetteExons,col_names=c("chrom","start", "end","name", "score", "strand", "thickStart", "thickEnd", "color" )) %>%
    mutate(end = end-1) %>%
    dplyr::select(name, chrom, start, end, strand) %>%
    separate(name, into=c("type", "dummy","cluster", "GAGTInt"), sep="\\_") %>%
    group_by(GAGTInt) %>%
    filter(n()==3) %>%
    ungroup() %>%
    pivot_wider(names_from="type", values_from=c("start", "end")) %>%
    mutate(CassetteExon_ExonStart = if_else(strand=="+", end_junc.UpstreamIntron, end_junc.GAGT)) %>%
    mutate(CassetteExon_ExonStop = if_else(strand=="+", start_junc.GAGT, start_junc.UpstreamIntron))

ExonTrios <- IntersectingExons %>%
  inner_join(TPM.transcripts.DMSO) %>%
  group_by(SkipJuncName, gene) %>%
  filter(Median.TPM == max(Median.TPM)) %>%
  ungroup() %>%
  filter(overlappingBases==1) %>%
  add_count(SkipJuncName, transcript) %>%
  filter(n==2) %>%
  # print(width=Inf) %>%
  mutate(ExonRelativePos = case_when(
    Strand == "+" & SkipJuncStart == ExonStop ~ "Upstream",
    Strand == "+" & SkipJuncStop == ExonStart ~ "Downstream",
    Strand == "-" & SkipJuncStart == ExonStop ~ "Downstream",
    Strand == "-" & SkipJuncStop == ExonStart ~ "Upstream",
    TRUE ~ "Other"
                                     )) %>%
  dplyr::select(ExonStart, ExonStop, gene, transcript, ExonRelativePos, SkipJuncName, Strand) %>%
  pivot_wider(names_from="ExonRelativePos", values_from=c("ExonStart", "ExonStop")) %>%
  separate(SkipJuncName, into=c("type", "dummy","cluster", "GAGTInt"), sep="\\_", remove=F) %>%
  dplyr::select(-type, -dummy, -cluster) %>%
  inner_join(
             cassette.exons %>%
                 dplyr::select(chrom, GAGTInt, ExonStart_Cassette=CassetteExon_ExonStart,ExonStop_Cassette=CassetteExon_ExonStop ),
             by="GAGTInt"
  ) %>%
  dplyr::select(chrom, GAGTInt, Strand, ExonStart_Upstream, ExonStop_Upstream, ExonStart_Cassette, ExonStop_Cassette, ExonStart_Downstream, ExonStop_Downstream, everything())

ExonTrios %>%
    pivot_longer(contains("Exon"), names_sep="\\_", names_to=c("StartOrStop", "Type")) %>%
    pivot_wider(names_from="StartOrStop", values_from="value") %>%
    unite(name, gene, transcript, Type, GAGTInt) %>%
    mutate(score = ".") %>%
    dplyr::select(chrom, ExonStart, ExonStop, name, score, Strand) %>%
    arrange(chrom, ExonStart, ExonStop) %>%
    write_tsv(f_long_out, col_names=F)

ExonTrios %>%
  write_tsv(f_out)

