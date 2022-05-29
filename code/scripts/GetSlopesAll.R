library(dplyr)
library(tidyverse)
library(data.table)
library(gtools)
library(broom)
library(MASS)
library(magrittr)
library(knitr)

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
intronFile = args[1]

SampleName <- strsplit(intronFile, "[.]")[[1]][1]
SampleName <- strsplit(SampleName, "/")[[1]][3]

windowStyle = strsplit(intronFile, "[.]")[[1]][2]

#dir.create(file.path('IntronSlopes', 'slopes', showWarnings = FALSE)

print(SampleName)
print(windowStyle)

minIntronCounts = as.integer(args[2])
minCoverageCounts = as.integer(args[3])
minCoverage = as.numeric(args[4])
WinLen = as.integer(args[5])
minIntronLen = as.integer(args[6])

print('loaded params')


print('maybe it crashes when loading data')
print(intronFile)
if (windowStyle == 'IntronWindows_equalLength')   {     
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
} else if  (windowStyle == 'IntronWindows'){
    Counts.EqualSizeBins <- fread(intronFile, #fread(paste(equalLengthDir, intronFile, sep=''), 
                              sep = '\t', col.names = c("chr", "start", "stop", "name", "score", "strand", "Readcount"), header=F) %>%
    separate(name, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand", "Window"),convert=T, sep = "_") %>%
    mutate(IntronLength = IntStop - IntStart) %>%
    unite(IntronName, gene, IntChr, IntStart, IntStop, IntStrand, remove=F) %>%
    group_by(IntronName) %>%
    mutate(WindowCount=max(Window)) %>%
    ungroup() %>%
    dplyr::select(chr, start, stop, Window, IntronName, Readcount, IntronLength)
}

print('or when processing')
            
IntronCountSums <- Counts.EqualSizeBins %>%
  group_by(IntronName) %>%
  summarise(IntronSum=sum(Readcount), IntronObsSum=sum(Readcount>0), IntronCovBins = sum(Readcount>=minCoverageCounts),
            IntronLength = mean(IntronLength), IntronCoverage=mean(Readcount>=minCoverageCounts))

print('or when filtering')
# Hard-coded IntronObsSum minimum of 2, despite min Intron Coverage, to ensure rlm can run
IntronsToGetSlope <- IntronCountSums %>% filter((IntronSum>=minIntronCounts), (IntronCovBins>=2), (IntronLength>=minIntronLen), ((IntronObsSum>=20) | (IntronCoverage>=minCoverage))) %>% pull(IntronName)

IntronCountSums %>% write_delim(paste('IntronSlopes/AllSlopes/', SampleName, '.', windowStyle, '.IntronCountSums.tab.gz', sep=''), delim = '\t')
            
print('or when whatever this is')

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

print('or getting coverage fits')
            
if (windowStyle == 'IntronWindows_equalLength')   {  
            
CoverageFits.df <- Counts.EqualSizeBins %>%
  filter(IntronName %in% IntronsToGetSlope) %>% #filter(Readcount >= 1) %>%
   ## See if we can make length independent
  #mutate(RelativeIntronPosInBp = Window*WinLen) %>% ##### This is what we were doing with equalLength
  mutate(RelativeIntronPosInBp =  (Window*WinLen)/IntronLength) %>%
  mutate(logCounts =  log1p(Readcount)) %>%
  #mutate(RelativeIntronPosInBp = Window/100) %>%
  # sample_n_of(9, IntronName) %>%
  group_by(IntronName) %>%
  mutate(NormCounts = scale(logCounts)) %>%
  do(CoverageFit = rlm(NormCounts ~ RelativeIntronPosInBp, data = ., maxit=40)) %>%
  summarise(IntronName, tidy(CoverageFit))
    
  } else if (windowStyle == 'IntronWindows')   {  
    CoverageFits.df <- Counts.EqualSizeBins %>%
    filter(IntronName %in% IntronsToGetSlope) %>% #filter(Readcount >= 1) %>%
    ## See if we can make length independent
    #mutate(RelativeIntronPosInBp = Window*WinLen) %>% ##### This is what we were doing with equalLength
    mutate(RelativeIntronPosInBp = Window/100) %>%
    #mutate(RelativeIntronPosInBp = Window/100) %>%
    # sample_n_of(9, IntronName) %>%
    group_by(IntronName) %>%
    mutate(NormCounts = scale(Readcount)) %>%
    do(CoverageFit = rlm(NormCounts ~ RelativeIntronPosInBp, data = ., maxit=40)) %>%
    summarise(IntronName, tidy(CoverageFit))
  }
            
print('or when tidying up')

#Let's tidy the data so there is one row per intron
CoverageFits.df.tidy <- CoverageFits.df %>%
  ungroup() %>%
  dplyr::select(-c("std.error", "statistic")) %>%
  spread(key="term", value="estimate") %>%
  separate(IntronName, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand"),convert=T, sep = "_", remove = F) %>%
  mutate(IntronLength = IntStop - IntStart)

print(table(CoverageFits.df.tidy$RelativeIntronPosInBp <= 0))

CoverageFits.error.tidy <- CoverageFits.df %>%
  ungroup() %>%
  dplyr::select(-c("estimate", "statistic")) %>%
  spread(key="term", value="std.error") %>%
  separate(IntronName, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand"),convert=T, sep = "_", remove = F) %>%
  mutate(IntronLength = IntStop - IntStart)

CoverageFits.df.tidy['std.error'] <- CoverageFits.error.tidy$RelativeIntronPosInBp

colnames(CoverageFits.df.tidy)[1] <- 'IntronName'
colnames(CoverageFits.df.tidy)[8] <- 'Slope'
            
IntronCountSums_filtered <- IntronCountSums[IntronCountSums$IntronName %in% CoverageFits.df.tidy$IntronName,]
IntronCountSums_filtered <- IntronCountSums_filtered %>% column_to_rownames(., var = "IntronName")
            
CoverageFits.df.tidy['counts'] <- IntronCountSums_filtered[CoverageFits.df.tidy$IntronName,]$IntronSum
CoverageFits.df.tidy['coverageMean'] <- IntronCountSums_filtered[CoverageFits.df.tidy$IntronName,]$IntronCoverage

CoverageFits.df.tidy <- CoverageFits.df.tidy %>%
  filter(IntronName %in% IntronsToGetSlope)
            
            
print(paste('IntronSlopes/AllSlopes/', SampleName, '.tab.gz', sep=''))

CoverageFits.df.tidy %>%
  mutate(IsSlopeNegative = Slope < 0 ) %>%
  write_delim(paste('IntronSlopes/AllSlopes/', SampleName, '.', windowStyle, '.tab.gz', sep=''), delim = '\t')
            





#mutate(RelativeIntronPosInBp =  Window*WinLen) %>%
# group_by(IntronName) %>%
#  mutate(NormCounts = scale(Readcount)) %>%