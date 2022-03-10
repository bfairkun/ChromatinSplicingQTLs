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

SampleName <- strsplit(intronFile, ".IntronWindows")[[1]][1]
SampleName <- strsplit(SampleName, "/")[[1]][3]

#dir.create(file.path('IntronSlopes', 'slopes', showWarnings = FALSE)

print(SampleName)

minIntronCounts = as.integer(args[2])
minCoverageCounts = as.integer(args[3])
minCoverage = as.numeric(args[4])
WinLen = as.integer(args[5])

print('loaded params')

tryCatch(
        {
    print('maybe it crashes when loading data')
            print(intronFile)
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

print('or when processing')
            
IntronCountSums <- Counts.EqualSizeBins %>%
  group_by(IntronName) %>%
  summarise(IntronSum=sum(Readcount), IntronObsSum=sum(Readcount>0), IntronCoverage=mean(Readcount>=minCoverageCounts))

print('or when filtering')
# Hard-coded IntronObsSum minimum of 2, despite min Intron Coverage, to ensure rlm can run
IntronsToGetSlope <- IntronCountSums %>% filter((IntronSum>=minIntronCounts), (IntronObsSum>=2), (IntronCoverage>=minCoverage)) %>% pull(IntronName)

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
            
CoverageFits.df <- Counts.EqualSizeBins %>%
  filter(IntronName %in% IntronsToGetSlope) %>% #filter(Readcount >= 1) %>%
   ## See if we can make length independent
  #mutate(RelativeIntronPosInBp = Window*WinLen) %>%
  mutate(RelativeIntronPosInBp = (Window*WinLen)/(stop - start)) %>%
  # sample_n_of(9, IntronName) %>%
  group_by(IntronName) %>%
  do(CoverageFit = rlm(Readcount ~ RelativeIntronPosInBp, data = ., maxit=40)) %>%
  summarise(IntronName, tidy(CoverageFit))
            
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
            
            
print(paste('IntronSlopes/slopes/', SampleName, '.tab.gz', sep=''))

CoverageFits.df.tidy %>%
  mutate(IsSlopeNegative = Slope < 0 ) %>%
  write_delim(paste('IntronSlopes/slopes/', SampleName, '.tab.gz', sep=''), delim = '\t')
            
},
    error = function(e) {
        
        data.frame(IntronName=character(),
                  gene=character(),
                  IntChr=character(),
                  IntStart=character(),
                  IntStop=character(),
                  IntStrand=character(),
                  Intercept=character(),
                  slope=character(),
                  IntronLength=character(),
                  std.error = character(),
                  counts = character(),
                  coverageMean = character(),
                  IsSlopeNegative = character(),
                  stringsAsFactors=FALSE) %>% write_delim(paste('slopes/', SampleName, '.tab.gz', sep=''), delim = '\t')
       
        })



#### GLM NB slopes


tryCatch(
        {
set.seed(0)

CoverageFits.nb.glm.df <- Counts.EqualSizeBins %>%
  filter(IntronName %in% IntronsToGetSlope) %>% #filter(Readcount >= 1) %>%
  mutate(RelativeIntronPosInBp = (Window*WinLen)/(stop - start)) %>%
  #mutate(RelativeIntronPosInBp = Window*WinLen) %>%
  group_by(IntronName) %>%
  do(CoverageFit = possibly(glm.nb, otherwise = NA)(Readcount ~ RelativeIntronPosInBp, data = ., link="identity")) %>%
  summarise(IntronName, tidy(CoverageFit)) #tidy(CoverageFit)

CoverageFits.df.tidy <- CoverageFits.nb.glm.df %>%
  ungroup() %>%
  dplyr::select(-c("std.error", "statistic", "p.value", "x")) %>%
  spread(key="term", value="estimate") %>%
  separate(IntronName, into = c("gene", "IntChr", "IntStart", "IntStop", "IntStrand"),convert=T, sep = "_", remove = F) %>%
  mutate(IntronLength = IntStop - IntStart) %>% 
     dplyr::select(c("IntronName", "gene", "IntChr", "IntStart", "IntStop", "IntStrand", "(Intercept)",
                    "RelativeIntronPosInBp", "IntronLength"))

colnames(CoverageFits.df.tidy) <- c("IntronName", "Gene", "IntChr", "IntStart", "IntStop", "IntStrand", "Intercept",
                    "Slope", "IntronLength")

CoverageFits.df.tidy['Slope.p.value'] = CoverageFits.nb.glm.df[CoverageFits.nb.glm.df$term == "RelativeIntronPosInBp",]$p.value
CoverageFits.df.tidy['Intercept.p.value'] = CoverageFits.nb.glm.df[CoverageFits.nb.glm.df$term == "(Intercept)",]$p.value

IntronCountSums_filtered <- IntronCountSums[IntronCountSums$IntronName %in% CoverageFits.df.tidy$IntronName,]
IntronCountSums_filtered <- IntronCountSums_filtered %>% column_to_rownames(., var = "IntronName")
            
CoverageFits.df.tidy['counts'] <- IntronCountSums_filtered[CoverageFits.df.tidy$IntronName,]$IntronSum
CoverageFits.df.tidy['coverageMean'] <- IntronCountSums_filtered[CoverageFits.df.tidy$IntronName,]$IntronCoverage
            
print(paste('IntronSlopes/slopes/', SampleName, '_glm.nb.tab.gz', sep=''))

CoverageFits.df.tidy %>%
  mutate(IsSlopeNegative = Slope < 0 ) %>%
  write_delim(paste('IntronSlopes/slopes/', SampleName, '_glm.nb.tab.gz', sep=''), delim = '\t')
            
},
    error = function(e) {
        
        data.frame(IntronName=character(),
                  gene=character(),
                  IntChr=character(),
                  IntStart=character(),
                  IntStop=character(),
                  IntStrand=character(),
                  Intercept=character(),
                  slope=character(),
                  IntronLength=character(),
                  std.error = character(),
                  counts = character(),
                  coverageMean = character(),
                  IsSlopeNegative = character(),
                  stringsAsFactors=FALSE) %>% write_delim(paste('IntronSlopes/slopes/', SampleName, '_glm.nb.tab.gz', sep=''), delim = '\t')
       
        })


