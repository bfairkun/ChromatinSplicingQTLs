library(tidyverse)
library(STAN)

args <- commandArgs(trailingOnly=TRUE)

chrom <- args[1]
strand <- args[2]



df <- read.table(paste0('NonCodingRNA/tables/', chrom, '.', strand, '.counts.tab.gz'), header=1)[,1:86]
df <- mutate_all(df, function(x) as.numeric(as.character(x)))

X = list('reg'=as.matrix(df))[1]
hmm_nb = initHMM(X, nStates=2, "PoissonLogNormal")

hmm_fitted_nb = fitHMM(X, hmm_nb, maxIters=200)
viterbi_nb = getViterbi(hmm_fitted_nb, X)

df['state'] <- viterbi_nb[[1]]
   
write.table(df['state'], paste0('NonCodingRNA/tables/', chrom, '.', strand, '.predicted.tab'), sep='\t', quote=FALSE, row.names=TRUE)