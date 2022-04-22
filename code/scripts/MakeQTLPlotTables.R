#!/usr/bin/env Rscript

if(interactive()){
    args <- scan(text=
                 "scratch/ ../output/hyprcoloc_results/ForColoc/MolColocStandard/hyprcoloc.results.txt.gz  'H3K27AC H3K4ME3 Expression.Splicing chRNA.Expression.Splicing'", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


PrefixOut <- args[1]
hyprcoloc_results_fn <- args[2]
Phenotypes <- unlist(strsplit(args[3], ' '))

library(tidyverse)
library(RColorBrewer)
library(scales)
library(data.table)
library(readxl)


ColocBed_fout <- paste0(PrefixOut, "ColocFeatures.bed")
TestedFeaturesBed_fout <- paste0(PrefixOut, "ColocTestFeatures.bed")
bwList_fout <-paste0(PrefixOut, "bwList.tsv")
bwGroupsFile_fout <- paste0(PrefixOut, "bwList.Groups.tsv")

# Read pre defined color key
# Colors for phenotypes, each of which is associated with a parent assay...
ColorsForPhenotypes <- read_excel("../data/ColorsForPhenotypes.xlsx", sheet = 1)
PhenotypePallette <-
  ColorsForPhenotypes %>%
  select(Phenotype, Hex) %>%
  deframe()
ColorsForPhenotypes %>%
  arrange(Phenotype) %>%
  mutate(n=row_number()) %>%
  ggplot(aes(x=1, y=n, fill=Phenotype, label=paste(Phenotype, Hex))) +
  geom_tile() +
  geom_text() +
  scale_fill_manual(values=PhenotypePallette) +
  theme_void()
PhenotypePallette_rgb <- col2rgb(PhenotypePallette) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Phenotype") %>%
  unite(rgb, red, green, blue, sep=",") %>% deframe()

# Colors for assays
ColorsForAssays <- read_excel("../data/ColorsForPhenotypes.xlsx", sheet = 2)
AssayPallette <-
  ColorsForAssays %>%
  select(ParentAssay, Hex) %>%
  deframe()
ColorsForAssays %>%
  arrange(ParentAssay) %>%
  mutate(n=row_number()) %>%
  ggplot(aes(x=1, y=n, fill=ParentAssay, label=paste(ParentAssay, Hex))) +
  geom_tile() +
  geom_text() +
  scale_fill_manual(values=AssayPallette) +
  theme_void()
AssayPallette_rgb <- col2rgb(AssayPallette) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Phenotype") %>%
  unite(rgb, red, green, blue, sep=",") %>% deframe()

# Colors for coloc clusters
ClusterBrewerPal <- "Accent"
ClusterBrewerPalN <- 8
ClusterPallette <- setNames(brewer.pal(ClusterBrewerPalN, ClusterBrewerPal), 1:ClusterBrewerPalN)
ClusterPallette_rgb <- col2rgb(ClusterPallette) %>% t() %>% as.data.frame() %>%
  rownames_to_column("ClusterColorNum") %>%
  unite(rgb, red, green, blue, sep=",") %>% deframe()

## Define phenotypes and their colors for plotting
coloc.results <- read_tsv(hyprcoloc_results_fn, col_names = c("Locus", "iteration", 'ColocalizedTraits', 'ColocPr', 'RegionalPr', "topSNP", "TopSNPFinemapPr", "DroppedTrait"), skip=1) %>%
  select(Locus, ColocalizedTraits, iteration, DroppedTrait, topSNP) %>%
  separate_rows(ColocalizedTraits, sep = ', ') %>%
  mutate(full_pid = case_when(
    ColocalizedTraits == "None" ~ DroppedTrait,
    TRUE ~ ColocalizedTraits
  )) %>%
  mutate(Phenotype=str_replace(full_pid, "(.+?);(.+)", "\\1")) %>%
  mutate(pid=str_replace(full_pid, "(.+?);(.+)", "\\2")) %>%
  select(-ColocalizedTraits, -DroppedTrait)

AllSamples <- read_tsv("config/samples.tsv")

Phenotypes <- coloc.results$Phenotype %>% unique()

## Read in bed of all features
beds <- setNames(paste0("QTLs/QTLTools/", Phenotypes ,"/OnlyFirstReps.sorted.qqnorm.bed.gz"), Phenotypes)

bed_list <- lapply(beds, fread, select=1:6)
combined_bed <- bind_rows(bed_list, .id = "Phenotype")


BedWithColors <- combined_bed %>%
  inner_join(coloc.results, by=c("Phenotype", "pid")) %>%
  mutate(PhenotypeColor = recode(Phenotype, !!!PhenotypePallette_rgb)) %>%
  arrange(Locus, iteration) %>%
  add_count(Locus, iteration) %>%
  mutate(ClusterColorNum = case_when(
    n>=2 ~ as.integer(iteration %% ClusterBrewerPalN + 1),
    TRUE ~ NA_integer_
  )) %>%
  mutate(ClusterColor = recode(ClusterColorNum, !!!ClusterPallette_rgb))

write_lines(c('#track name="ColocalizedFeatures" itemRgb="On"'), ColocBed_fout )
BedWithColors %>%
  filter(!is.na(ClusterColorNum)) %>%
  mutate(score=paste(Locus,topSNP, sep="_"), strand=".", thickStart=start, thickEnd=end) %>%
  select(`#Chr`, start, end, full_pid, score, strand, thickStart, thickEnd, ClusterColor) %>%
  arrange(`#Chr`, start, end) %>%
  write_tsv(ColocBed_fout, append = T)

write_lines(c('#track name="TestFeatures" itemRgb="On"'), TestedFeaturesBed_fout )
BedWithColors %>%
  mutate(full_pid_with_cluster = paste(Locus, topSNP, full_pid, sep="_")) %>%
  mutate(score="0", strand=".", thickStart=start, thickEnd=end) %>%
  select(`#Chr`, start, end, full_pid, score, strand, thickStart, thickEnd, PhenotypeColor) %>%
  distinct(full_pid, .keep_all = T) %>%
  arrange(`#Chr`, start, end) %>%
  write_tsv(TestedFeaturesBed_fout, append = T)


# Read geuvadis YRI
YRI <- read_tsv("QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.sorted.qqnorm.bed.gz") %>%
  select(-(1:6)) %>% colnames()

# make file of bigwigs
AllSamples$Phenotype %>% unique()
BwTable.df <- AllSamples %>%
  select(IndID, Assay, Phenotype, RepNumber) %>%
  filter(RepNumber==1) %>%
  filter(Phenotype %in% c(Phenotypes, "Expression.Splicing" )) %>%
  distinct(.keep_all = T) %>%
  mutate(Strand = case_when(
    Phenotype == "chRNA.Expression.Splicing" ~ "+,-",
    TRUE ~ ".")) %>%
  separate_rows(Strand, sep=",") %>%
  mutate(bw= case_when(
    Strand == "+" ~ str_glue("bigwigs/{Phenotype}_stranded/{IndID}.{RepNumber}.minus.bw"),
    Strand == "-" ~ str_glue("bigwigs/{Phenotype}_stranded/{IndID}.{RepNumber}.plus.bw"),
    Strand == "." ~ str_glue("bigwigs/{Phenotype}/{IndID}.{RepNumber}.bw")
  )) %>%
  filter((!Phenotype=="Expression.Splicing") | (IndID %in% YRI)) %>%
  select(IndID, bw, Phenotype, Strand) %>%
  mutate(Phenotype = recode(Phenotype, Expression.Splicing="polyA.RNA", chRNA.Expression.Splicing="chRNA")) %>%
  select(SampleID=IndID, BigwigFilepath=bw, Group_label=Phenotype, Strand)

write_tsv(BwTable.df, bwList_fout)

# Make group-specific settings
# 1: Group_label (must match a Group_label in the KeyFile) 2: Group_color (optional). Hex or rgb colors to plot potentially plot each group as defined in the output ini 3: BedgzFile (optional).
BwTable.df %>%
  distinct(Group_label, .keep_all = T) %>%
  mutate(Phenotype = str_replace(BigwigFilepath, ".+?/(.+?)[/_].+", "\\1")) %>%
  mutate(Group_color = recode(Phenotype, !!!AssayPallette)) %>%
  mutate(BedgzFile = case_when(
    Phenotype %in% c("Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing") ~ paste0("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.", Phenotype, ".bed.gz"),
    TRUE ~ ""
  )) %>%
  select(Group_label, Group_color, BedgzFile) %>%
  write_tsv(bwGroupsFile_fout)

