library(tidyverse)

#setwd("CurrentProjects/ChromatinSplicingQTLs/code/")



GEUVADIS.eQTLs %>% 
  dplyr::select(SNP_chr, SNP_pos) %>% dplyr::distinct() %>% dim()

Grubert.H3K4me3.QTLs <- read.delim("PlotGruberQTLs/Data/localQTLs_H3K4ME3.FDR0.1.hg38.bedpe", col.names=c("SNP_chr", "SNP_pos", "SNP_stop", "Peak_Chr", "Peak_start", "Peak_stop", "Name", "Score", "strand1", "stand2"), stringsAsFactors = F) %>%
  separate(col = "Name", into=c("PEAKid","SNPrsid","beta","p.value","FDR_TH","pvalTH","pass.pvalTH","mod","peak.state"), sep = ";", convert = T) %>%
  # filter(peak.state=="TSS") %>%
  # filter(pass.pvalTH=="pass") %>%
  mutate(SNP_pos = as.numeric(SNP_pos))



GEUVADIS.eQTLs <- read.delim("../data/QTLBase.GEUVADIS.eQTLs.hg38.txt.gz", stringsAsFactors = F) %>%
  mutate(SNP_chr=paste0("chr",SNP_chr)) %>%
  filter(!SNP_chr == "chr6") %>%
  mutate(SNP_pos = as.numeric(SNP_pos))  %>%
  mutate(SNP=paste(SNP_chr, SNP_pos))

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

H3K4me3_QQ <- Grubert.H3K4me3.QTLs %>%
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
H3K4me3_QQ
ggsave("PlotGruberQTLs/Plots/H3K4Me3_QQPlot_400TestsSampledPerGroup.png")

Grubert.H3K4me3.QTLs %>%
  dplyr::select(SNP_chr, SNP_pos) %>% dplyr::distinct() %>% dim()
  
ShellScript.A <- "PlotGruberQTLs/Plots/PlotH3K4Me3QTL_eQTLs.sh"
cat("set +xe",file=ShellScript.A,sep="\n")
Grubert.H3K4me3.QTLs %>%
  filter(peak.state=="TSS") %>%
  inner_join(GEUVADIS.eQTLs, by=c("SNP_chr","SNP_pos")) %>%
  filter(pass.pvalTH=="pass") %>%
  group_by(PEAKid) %>%
  dplyr::select(p.value, SNP_chr, SNP_pos, PEAKid, Peak_start, Peak_stop, beta, Mapped_gene) %>%
  slice(which.min(p.value)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(LeftBorder=min(SNP_pos, Peak_start)-1000, RightBorder=max(SNP_pos, Peak_stop)+1000) %>%
  mutate(NewSnpPos=SNP_pos+1) %>%
  mutate(ExpressionMakeTracks=glue::glue('python scripts/NormalizedBigwigsByGenotype.py Genotypes/GEUVADIS_Lappalainnen.vcf.gz {SNP_chr}:{NewSnpPos} {SNP_chr}:{LeftBorder}-{RightBorder} "Bigwigs/Grubert_ChIPSeq/H3K4ME3.GM1*.bw" --Normalization WholeGenome --BigwigListType GlobPattern --OutputPrefix PlotGruberQTLs/Data/Test/H3K4me3QTL_eQTL_Peak_{PEAKid}. --TracksTemplate scratch/tracks.ini.template3.txt')) %>%
  mutate(ExpressionPlotTracks=glue::glue('pyGenomeTracks --tracks <(cat PlotGruberQTLs/Data/Test/H3K4me3QTL_eQTL_Peak_{PEAKid}.output_tracks.ini ../data/pygenometracks/gtf_tracks.ini ../data/pygenometracks/GruberH3K4me3.tracks.ini) --out PlotGruberQTLs/Plots/H3K4me3QTL_eQTL_Peak_{PEAKid}.pdf --region {SNP_chr}:{LeftBorder}-{RightBorder}')) %>%
  dplyr::select(ExpressionMakeTracks, ExpressionPlotTracks) %>%
  gather(key="ExpressionType", value="Expression") %>%
  arrange(ExpressionType) %>%
  dplyr::select(Expression) %>%
  write.table(ShellScript.A, quote = F, row.names = F,col.names = F, append=T)
  
ShellScript.B <- "PlotGruberQTLs/Plots/PlotH3K4Me3TopQTL.sh"
cat("set +xe",file=ShellScript.B,sep="\n")
Grubert.H3K4me3.QTLs %>%
  filter(pass.pvalTH=="pass") %>%
  group_by(PEAKid) %>%
  dplyr::select(p.value, SNP_chr, SNP_pos, PEAKid, Peak_start, Peak_stop, beta) %>%
  slice(which.min(p.value)) %>%
  ungroup() %>%
  arrange(p.value) %>% head(10) %>%
  rowwise() %>%
  mutate(LeftBorder=min(SNP_pos, Peak_start)-1000, RightBorder=max(SNP_pos, Peak_stop)+1000) %>%
  mutate(NewSnpPos=SNP_pos+1) %>%
  mutate(ExpressionMakeTracks=glue::glue('python scripts/NormalizedBigwigsByGenotype.py Genotypes/GEUVADIS_Lappalainnen.vcf.gz {SNP_chr}:{NewSnpPos} {SNP_chr}:{LeftBorder}-{RightBorder} "Bigwigs/Grubert_ChIPSeq/H3K4ME3.GM1*.bw" --Normalization WholeGenome --BigwigListType GlobPattern --OutputPrefix PlotGruberQTLs/Data/Test/H3K4me3QTL_eQTL_Peak_{PEAKid}. --TracksTemplate scratch/tracks.ini.template3.txt')) %>%
  mutate(ExpressionPlotTracks=glue::glue('pyGenomeTracks --tracks <(cat PlotGruberQTLs/Data/Test/H3K4me3QTL_eQTL_Peak_{PEAKid}.output_tracks.ini ../data/pygenometracks/gtf_tracks.ini ../data/pygenometracks/GruberH3K4me3.tracks.ini) --out PlotGruberQTLs/Plots/H3K4me3TopQTL_Peak_{PEAKid}.pdf --region {SNP_chr}:{LeftBorder}-{RightBorder}')) %>%
  dplyr::select(ExpressionMakeTracks, ExpressionPlotTracks) %>%
  gather(key="ExpressionType", value="Expression") %>%
  arrange(ExpressionType) %>%
  dplyr::select(Expression) %>%
  write.table(ShellScript.B, quote = F, row.names = F,col.names = F, append=T)
