#!/usr/bin/env Rscript

if(interactive()){
    args <- scan(text=
                 "scratch/ scratch/hyprcoloc.results.txt.gz SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz 'H3K27AC H3K4ME3 Expression.Splicing chRNA.Expression.Splicing'", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

PrefixOut <- args[1]
hyprcoloc_results_fn <- args[2]
leafcutter_count_numers <- args[3]
Phenotypes <- unlist(strsplit(args[4], ' '))
## Define phenotypes and their colors for plotting
Phenotypes <- Sys.glob("QTLs/QTLTools/*/PermutationPassForColoc.txt.gz") %>%
  str_replace("QTLs/QTLTools/(.+?)/.+", "\\1")

library(tidyverse)
library(RColorBrewer)
library(scales)
library(data.table)


ColocBed_fout <- paste0(PrefixOut, "ColocFeatures.bed")
TestedFeaturesBed_fout <- paste0(PrefixOut, "ColocTestFeatures.bed")
bwList_fout <-paste0(PrefixOut, "bwList.tsv")
PSI_fout_prefix <-paste0(PrefixOut, "PSI.")
bwGroupsFile_fout <- paste0(PrefixOut, "bwList.Groups.tsv")


# display.brewer.all()
# display.brewer.pal(12, "Set3")


ColorBrewerPal <- "Paired"
PhenotypePallette <- colorRampPalette(brewer.pal(12, ColorBrewerPal))(length(Phenotypes))
PhenotypePallette <- setNames(PhenotypePallette, Phenotypes)
# show_col(PhenotypePallette)
PhenotypePallette_rgb <- col2rgb(PhenotypePallette) %>% t() %>% as.data.frame() %>%
  rownames_to_column("Phenotype") %>%
  unite(rgb, red, green, blue, sep=",") %>% deframe()

ColorKeyPlot <- data.frame(PhenotypePallette) %>%
  rownames_to_column("Phenotype") %>%
  mutate(n=row_number()) %>%
  ggplot(aes(x=1, y=n, fill=Phenotype, label=Phenotype)) +
  geom_tile() +
  geom_text() +
  scale_fill_manual(values=PhenotypePallette) +
  theme_void()
ggsave(paste0(PrefixOut, "PhenotypeColorLegend.pdf"), ColorKeyPlot)

ClusterBrewerPal <- "Accent"
ClusterBrewerPalN <- 8
ClusterPallette <- setNames(brewer.pal(ClusterBrewerPalN, ClusterBrewerPal), 1:ClusterBrewerPalN)
ClusterPallette_rgb <- col2rgb(ClusterPallette) %>% t() %>% as.data.frame() %>%
  rownames_to_column("ClusterColorNum") %>%
  unite(rgb, red, green, blue, sep=",") %>% deframe()



## Read in bed of all features
beds <- setNames(paste0("QTLs/QTLTools/", Phenotypes ,"/OnlyFirstReps.sorted.qqnorm.bed.gz"), Phenotypes)

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


# OnlyColocalized <- read_tsv("../output/hyprcoloc_results/ForColoc/hyprcoloc.results.OnlyColocalized.Stats.txt.gz")


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
  mutate(score=".", strand=".", thickStart=start, thickEnd=end) %>%
  select(`#Chr`, start, end, full_pid, score, strand, thickStart, thickEnd, PhenotypeColor) %>%
  distinct(full_pid, .keep_all = T) %>%
  arrange(`#Chr`, start, end) %>%
  write_tsv(TestedFeaturesBed_fout, append = T)

# Make PSI tables
Count.Table.mat <- read.table(leafcutter_count_numers, sep = ' ', nrows = Inf) %>%
  as.matrix()

ClusterMax.mat <- Count.Table.mat %>%
    as.data.frame() %>%
    rownames_to_column("junc") %>%
    mutate(cluster=str_replace(junc, "^(.+?):.+?:.+?:(.+)$", "\\1_\\2")) %>%
    group_by(cluster) %>%
    mutate(across(where(is.numeric), max)) %>%
    ungroup() %>%
    select(junc, everything(), -cluster) %>%
    column_to_rownames("junc") %>%
    as.matrix()


PSI.df <- (Count.Table.mat / as.numeric(ClusterMax.mat) * 100) %>%
    as.data.frame()

CountTablePhenotypes <- colnames(Count.Table.mat)[-1] %>%
    str_replace("^(.+?)_.+?_.+$", "\\1") %>% unique()

for (p in CountTablePhenotypes){
    PSI.df %>%
        rownames_to_column("junc") %>%
        select(junc, starts_with(p) & ends_with("_1")) %>%
        rename_with(~ str_replace(.x, "^.+?_(.+?)_.+$", "\\1"), starts_with(p)) %>%
        separate(junc, into=c("#Chrom", "start", "end", "cluster"), convert=T, remove=F, sep=":") %>%
        mutate(gid = paste(`#Chrom`, cluster, sep="_" )) %>%
        mutate(strand = str_extract(cluster, "[+-]")) %>%
        select(`#Chrom`, start, end, junc, gid, strand, everything(), -cluster) %>%
        arrange(`#Chrom`, start, end) %>%
        write_tsv(paste0(PSI_fout_prefix, p, ".bed"))
}


# Read geuvadis YRI
YRI <- read_tsv("QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.sorted.qqnorm.bed.gz") %>%
  select(-(1:6)) %>% colnames()

# make file of bigwigs
AllSamples <- read_tsv("config/samples.tsv")
AllSamples$Phenotype %>% unique()
BwTable.df <- AllSamples %>%
  select(IndID, Assay, Phenotype, RepNumber) %>%
  filter(RepNumber==1) %>%
  filter(Phenotype %in% Phenotypes) %>%
  distinct(.keep_all = T) %>%
  mutate(Strand = case_when(
    # Phenotype == "chRNA.Expression.Splicing" ~ "+,-",
    TRUE ~ ".")) %>%
  separate_rows(Strand, sep=",") %>%
  mutate(bw= case_when(
    Strand == "+" ~ str_glue("bigwigs/{Phenotype}_stranded/{IndID}.{RepNumber}.minus.bw"),
    Strand == "-" ~ str_glue("bigwigs/{Phenotype}_stranded/{IndID}.{RepNumber}.plus.bw"),
    Strand == "." ~ str_glue("bigwigs/{Phenotype}/{IndID}.{RepNumber}.bw")
  )) %>%
  filter((!Phenotype=="Expression.Splicing") | (IndID %in% YRI)) %>%
  select(IndID, bw, Phenotype, Strand) %>%
  mutate(Phenotype = recode(Phenotype, Expression.Splicing="polyA RNA", chRNA.Expression.Splicing="chRNA")) %>%
  select(SampleID=IndID, BigwigFilepath=bw, Group_label=Phenotype, Strand)

write_tsv(BwTable.df, bwList_fout)

# Make group-specific settings
# 1: Group_label (must match a Group_label in the KeyFile) 2: Group_color (optional). Hex or rgb colors to plot potentially plot each group as defined in the output ini 3: BedgzFile (optional).
BwTable.df %>%
  distinct(Group_label, .keep_all = T) %>%
  mutate(Phenotype = str_replace(BigwigFilepath, ".+?/(.+?)[/_].+", "\\1")) %>%
  mutate(Group_color = recode(Phenotype, !!!PhenotypePallette)) %>%
  mutate(BedgzFile = case_when(
    Phenotype %in% CountTablePhenotypes ~ paste0(PSI_fout_prefix, Phenotype, ".bed.gz"),
    TRUE ~ ""
  )) %>%
  select(Group_label, Group_color, BedgzFile) %>%
  write_tsv(bwGroupsFile_fout)

