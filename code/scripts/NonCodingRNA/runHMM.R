print('loading tidyverse')
library(tidyverse)
print('loading STAN')
library(STAN)

print('Running...')

args <- commandArgs(trailingOnly=TRUE)

chrom <- args[1]
strand <- args[2]
#states <- as.integer(args[3])
tables_dir <- args[3]



df <- read.table(paste0(tables_dir, chrom, '.', strand, '.counts.tab.gz'), header=1)#[,1:86]
df <- mutate_all(df, function(x) as.numeric(as.character(x)))

X = list('reg'=as.matrix(df))[1]
#hmm_nb = initHMM(X, nStates=states, "PoissonLogNormal")
hmm_nb = initHMM(X, nStates=5, "PoissonLogNormal")                 

hmm_fitted_nb = fitHMM(X, hmm_nb, maxIters=200)
viterbi_nb = getViterbi(hmm_fitted_nb, X)

df['state'] <- viterbi_nb[[1]]
   
write.table(df['state'], paste0(tables_dir, chrom, '.', strand, '.predicted_states.tab'), 
            sep='\t', quote=FALSE, row.names=TRUE)
