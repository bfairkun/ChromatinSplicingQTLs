library(plyr)
library(tidyverse)
library(data.table)
library(gtools)
library(broom)
library(MASS)
library(magrittr)
library(knitr)

options(warn=-1)
# equalLengthDir = '/project2/yangili1/cfbuenabadn/ChromatinSplicing/code/IntronWindowCounts/'
# intronCoverageFiles = list.files(equalLengthDir)

# intronCoverageFiles <- intronCoverageFiles[grepl('_equalLength', intronCoverageFiles)]
       
# for (intronFile in intronCoverageFiles) {
    
#     tryCatch({

args = commandArgs(trailingOnly=TRUE)
intronFile = args[1]

minIntronCounts = as.integer(args[2])
minCoverageCounts = as.integer(args[3])
minCoverage = as.numeric(args[4])
WinLen = as.integer(args[5])

print(args[5])

print(WinLen)



# intronFile = 'IntronWindowCounts/NA19209.IntronWindows_equalLength.bed.gz' 

# minIntronCounts = 10
# minCoverageCounts = 1
# minCoverage = 0.25
# WinLen = 100

print(args[5])


    
Counts.EqualSizeBins <- fread(intronFile, #fread(paste(equalLengthDir, intronFile, sep=''), 
                              sep = '\t', col.names = c("chr", "start", "stop", "name", "score", "strand", "Readcount"), header=F) %>%

  separate(name, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand", "Window"),convert=T, sep = "_") %>%
  mutate(IntronLength = IntStop - IntStart) %>%
  unite(IntronName, gene, IntChr, IntStart, IntStop, IntStrand, remove=F) %>%
  group_by(IntronName) %>%
  mutate(WindowCount=max(Window)) %>%
  ungroup() %>%
  mutate(Window = case_when(
    IntStrand == "+" ~ Window,
    IntStrand == "-" ~ WindowCount + as.integer(1) - Window
  )) %>%
  dplyr::select(chr, start, stop, Window, IntronName, Readcount, IntronLength)


IntronCountSums <- Counts.EqualSizeBins %>%
  group_by(IntronName) %>%
  summarise(IntronSum=sum(Readcount), IntronObsSum=sum(Readcount>0), IntronCoverage=mean(Readcount>=minCoverageCounts))


#This plot should look similar if not identical...
ggplot(IntronCountSums, aes(x=IntronSum)) +
  geom_histogram() +
  scale_x_continuous(trans="log10") +
  theme_bw()

# Hard-coded IntronObsSum minimum of 2, despite min Intron Coverage, to ensure rlm can run
IntronsToGetSlope <- IntronCountSums %>% filter((IntronSum>=minIntronCounts), (IntronObsSum>=2), (IntronCoverage>=minCoverage)) %>% pull(IntronName)


sample_n_of <- function(data, size, ...) {
  dots <- quos(...)

  group_ids <- data %>% 
    group_by(!!! dots) %>% 
    group_indices()

  sampled_groups <- sample(unique(group_ids), size)

  data %>% 
    filter(group_ids %in% sampled_groups)
}

set.seed(0)
# Counts.EqualSizeBins %>%
#   filter(IntronName %in% IntronsToGetSlope) %>%
#   mutate(RelativeIntronPosInBp = Window*50) %>%
#   sample_n_of(9, IntronName) %>%
#   ggplot(aes(x=RelativeIntronPosInBp, y=Readcount)) +
#   geom_point() +
#   geom_smooth(method='rlm') +
#   # geom_smooth(method='glm.nb') +
#   facet_wrap(~IntronName, nrow=3, scales="free") +
#   theme_bw()

# Counts.EqualSizeBins %>%
#   filter(IntronName %in% IntronsToGetSlope) %>%
#   #Add a pseudocount so we don't filter out windows with 0 read counts if we try to plot on a log scale
#   mutate(Readcount = Readcount + 0.01)  %>%
#   ggplot(aes(x=Readcount)) +
#   geom_histogram() +
#   scale_x_continuous(limits = c(0,250)) +
#   theme_bw()

# Counts.EqualSizeBins %>%
#   filter(IntronName %in% IntronsToGetSlope) %>%
#   sample_n_of(9, IntronName) %>%
#   #Add a pseudocount so we don't filter out windows with 0 read counts if we try to plot on a log scale
#   mutate(Readcount = Readcount + 0.01)  %>%
#   ggplot(aes(x=Readcount)) +
#   geom_histogram() +
#   scale_x_continuous() +
#   facet_wrap(~IntronName, nrow=3, scales="free") +
#   theme_bw()

CoverageFits.df <- Counts.EqualSizeBins %>%
  filter(IntronName %in% IntronsToGetSlope) %>% filter(Readcount >= 1) %>%
   ## See if we can make length independent
  mutate(RelativeIntronPosInBp = Window*WinLen) %>%
  # sample_n_of(9, IntronName) %>%
  group_by(IntronName) %>%
# mutate(RelativeIntronPosInBp = seq(1/n(), 1, by=1/n())) %>% #### THIS IS NEW
  do(CoverageFit = rlm(Readcount ~ RelativeIntronPosInBp, data = ., maxit=40)) %>%
  tidy(CoverageFit)

#Let's tidy the data so there is one row per intron
CoverageFits.df.tidy <- CoverageFits.df %>%
  ungroup() %>%
  dplyr::select(-c("std.error", "statistic")) %>%
  spread(key="term", value="estimate") %>%
  separate(IntronName, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand"),convert=T, sep = "_", remove = F) %>%
  mutate(IntronLength = IntStop - IntStart)

print(table(CoverageFits.df.tidy$RelativeIntronPosInBp <= 0))

SampleName <- strsplit(intronFile, ".IntronWindows")[[1]][1]
SampleName <- strsplit(SampleName, "/")[[1]][2]

CoverageFits.error.tidy <- CoverageFits.df %>%
  ungroup() %>%
  dplyr::select(-c("estimate", "statistic")) %>%
  spread(key="term", value="std.error") %>%
  separate(IntronName, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand"),convert=T, sep = "_", remove = F) %>%
  mutate(IntronLength = IntStop - IntStart)

CoverageFits.df.tidy['std.error'] <- CoverageFits.error.tidy$RelativeIntronPosInBp

colnames(CoverageFits.df.tidy)[8] <- 'slope'

CoverageFits.df.tidy['counts'] <- IntronCountSums[rownames(CoverageFits.df.tidy),]$IntronSum
CoverageFits.df.tidy['coverageMean'] <- IntronCountSums[rownames(CoverageFits.df.tidy),]$IntronCoverage

CoverageFits.df.tidy <- CoverageFits.df.tidy %>%
  filter(IntronName %in% IntronsToGetSlope)

print(paste('slopes/', SampleName, '.tab.gz', sep=''))

CoverageFits.df.tidy %>%
  mutate(IsSlopeNegative = slope <= 0 ) %>%
  write_delim(paste('slopes/', SampleName, '.tab.gz', sep=''), delim = '\t')


             
