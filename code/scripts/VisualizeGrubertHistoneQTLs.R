library(tidyverse)

# setwd("~/CurrentProjects/ChromatinSplicingQTLs/code/")

args = commandArgs(trailingOnly=TRUE)
Pheno = args[1]


# Pheno = "H3K4ME1"



Grubert.QTLs <- read.delim(paste0("PlotGruberQTLs/Data/localQTLs_", Pheno, ".FDR0.1.hg38.bedpe"), col.names=c("SNP_chr", "SNP_pos", "SNP_stop", "Peak_Chr", "Peak_start", "Peak_stop", "Name", "Score", "strand1", "stand2"), stringsAsFactors = F) %>%
  separate(col = "Name", into=c("PEAKid","SNPrsid","beta","p.value","FDR_TH","pvalTH","pass.pvalTH","mod","peak.state"), sep = ";", convert = T) %>%
  # filter(peak.state=="TSS") %>%
  # filter(pass.pvalTH=="pass") %>%
  mutate(SNP_pos = as.numeric(SNP_pos)) %>%
  mutate(SNP=paste(SNP_chr, SNP_pos))

GEUVADIS.eQTLs <- read.delim("../data/QTLBase.GEUVADIS.eQTLs.hg38.txt.gz", stringsAsFactors = F) %>%
  mutate(SNP_chr=paste0("chr",SNP_chr)) %>%
  filter(!SNP_chr == "chr6") %>%
  mutate(SNP_pos = as.numeric(SNP_pos)) %>%
  mutate(SNP=paste(SNP_chr, SNP_pos))


  
# GEUVADIS.eQTLs %>% 
#   inner_join(Grubert.QTLs, by="SNP")  %>%
#   filter(pass.pvalTH=="pass") %>% dim()
  


# H3K4me3_QQ <- Grubert.H3K4me3.QTLs %>%
#   dplyr::select(p.value, SNP_chr, SNP_chr, SNP_pos) %>%
#   filter(!SNP_chr == "chr6") %>%
#   mutate(SNP=paste(SNP_chr, SNP_pos)) %>%
#   mutate(SNPIsEqtl= SNP %in% GEUVADIS.eQTLs$SNP) %>%
#   group_by(SNPIsEqtl) %>%
#   mutate(ExpectedP=percent_rank(p.value)) %>%
#   ggplot(aes(x=-log10(ExpectedP), y=-log10(p.value), color=SNPIsEqtl)) +
#   geom_point() +
#   geom_abline() +
#   theme_bw()

QQ.Plot <- Grubert.QTLs %>%
  dplyr::select(p.value, SNP_chr, SNP_chr, SNP_pos) %>%
  filter(!SNP_chr == "chr6") %>%
  mutate(SNP=paste(SNP_chr, SNP_pos)) %>%
  mutate(SNPIsEqtl= SNP %in% GEUVADIS.eQTLs$SNP) %>%
  group_by(SNPIsEqtl) %>%
  mutate(ExpectedP=percent_rank(p.value)) %>%
  sample_n(400) %>%
  ggplot(aes(x=-log10(ExpectedP), y=-log10(p.value), color=SNPIsEqtl)) +
  geom_point() +
  geom_abline() +
  theme_bw()
QQ.Plot
ggsave(paste0("PlotGruberQTLs/Plots/", Pheno, "_QQPlot_400TestsSampledPerGroup.png"))

Grubert.QTLs %>%
  dplyr::select(SNP_chr, SNP_pos) %>% dplyr::distinct() %>% dim()
  
ShellScript.A <- paste0("PlotGruberQTLs/Plots/Plot", Pheno, "QTL_eQTLs.sh")
cat("set +xe",file=ShellScript.A,sep="\n")
Grubert.QTLs %>%
  # filter(peak.state=="TSS") %>%
  inner_join(GEUVADIS.eQTLs, by=c("SNP_chr","SNP_pos")) %>%
  filter(pass.pvalTH=="pass") %>%
  group_by(PEAKid) %>%
  dplyr::select(p.value, SNP_chr, SNP_pos, PEAKid, Peak_start, Peak_stop, beta, Mapped_gene) %>%
  slice(which.min(p.value)) %>%
  ungroup() %>% head(10) %>%
  rowwise() %>%
  mutate(LeftBorder=min(SNP_pos, Peak_start)-1000, RightBorder=max(SNP_pos, Peak_stop)+1000) %>%
  mutate(NewSnpPos=SNP_pos+1) %>%
  mutate(Pheno=Pheno) %>%
  mutate(ExpressionMakeTracks=glue::glue('python scripts/GenometracksByGenotype/NormalizedBigwigsByGenotype.py Genotypes/GEUVADIS_Lappalainnen.vcf.gz {SNP_chr}:{NewSnpPos} {SNP_chr}:{LeftBorder}-{RightBorder} "Bigwigs/Grubert_ChIPSeq/{Pheno}.GM1*.bw" --Normalization WholeGenome --BigwigListType GlobPattern --OutputPrefix PlotGruberQTLs/Data/Tracks/{Pheno}QTL_eQTL_Peak_{PEAKid}. --TracksTemplate scripts/GenometracksByGenotype/tracks_templates/tracks.ini.template3.txt --OutputNormalizedBigwigsPerSample')) %>%
  mutate(ExpressionPlotTracks=glue::glue('pyGenomeTracks --tracks <(cat PlotGruberQTLs/Data/Tracks/{Pheno}QTL_eQTL_Peak_{PEAKid}.output_tracks.ini ../data/pygenometracks/gtf_tracks.ini PlotGruberQTLs/Data/Tracks/{Pheno}.tracks.ini) --out PlotGruberQTLs/Plots/{Pheno}QTL_eQTL_Peak_{PEAKid}.pdf --region {SNP_chr}:{LeftBorder}-{RightBorder}')) %>%
  dplyr::select(ExpressionMakeTracks, ExpressionPlotTracks) %>%
  gather(key="ExpressionType", value="Expression") %>%
  arrange(ExpressionType) %>%
  dplyr::select(Expression) %>%
  write.table(ShellScript.A, quote = F, row.names = F,col.names = F, append=T)
  
ShellScript.B <- paste0("PlotGruberQTLs/Plots/Plot", Pheno, "TopQTL.sh")
cat("set +xe",file=ShellScript.B,sep="\n")
Grubert.QTLs %>%
  filter(pass.pvalTH=="pass") %>%
  group_by(PEAKid) %>%
  dplyr::select(p.value, SNP_chr, SNP_pos, PEAKid, Peak_start, Peak_stop, beta) %>%
  slice(which.min(p.value)) %>%
  ungroup() %>%
  arrange(p.value) %>% head(10) %>%
  rowwise() %>%
  mutate(LeftBorder=min(SNP_pos, Peak_start)-1000, RightBorder=max(SNP_pos, Peak_stop)+1000) %>%
  mutate(NewSnpPos=SNP_pos+1) %>%
  mutate(Pheno=Pheno) %>%
  mutate(ExpressionMakeTracks=glue::glue('python scripts/NormalizedBigwigsByGenotype.py Genotypes/GEUVADIS_Lappalainnen.vcf.gz {SNP_chr}:{NewSnpPos} {SNP_chr}:{LeftBorder}-{RightBorder} "Bigwigs/Grubert_ChIPSeq/{Pheno}.GM1*.bw" --Normalization WholeGenome --BigwigListType GlobPattern --OutputPrefix PlotGruberQTLs/Data/Tracks/{Pheno}TopQTL_Peak_{PEAKid}. --TracksTemplate scripts/GenometracksByGenotype/tracks_templates/tracks.ini.template3.txt --OutputNormalizedBigwigsPerSample')) %>%
  mutate(ExpressionPlotTracks=glue::glue('pyGenomeTracks --tracks <(cat PlotGruberQTLs/Data/Tracks/{Pheno}TopQTL_Peak_{PEAKid}.output_tracks.ini ../data/pygenometracks/gtf_tracks.ini PlotGruberQTLs/Data/Tracks/{Pheno}.tracks.ini) --out PlotGruberQTLs/Plots/{Pheno}TopQTL_Peak_{PEAKid}.pdf --region {SNP_chr}:{LeftBorder}-{RightBorder}')) %>%
  dplyr::select(ExpressionMakeTracks, ExpressionPlotTracks) %>%
  gather(key="ExpressionType", value="Expression") %>%
  arrange(ExpressionType) %>%
  dplyr::select(Expression) %>%
  write.table(ShellScript.B, quote = F, row.names = F,col.names = F, append=T)
