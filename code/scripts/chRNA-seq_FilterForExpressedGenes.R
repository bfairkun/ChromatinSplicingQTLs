library(tidyverse)


# Rscript to filter introns bed for introns in expressed host genes


args = commandArgs(trailingOnly=TRUE)
Gtf = args[1]
Intron.bed = args[2]
ExpressedGenesList = args[3]
Output = args[4]



Gtf <- read_delim(Gtf, comment='#', delim='\t', col_names=c("Chr", "Source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

GeneAnnotation <- Gtf %>%
    filter(feature=="transcript") %>%
    mutate(gene = str_replace(attribute, 'gene_id "(.+?)".+', "\\1")) %>%
    mutate(transcript = str_replace(attribute, '.+?transcript_id "(.+?)".+$', "\\1")) %>%
    dplyr::select(gene, transcript)

Intron.bed  <- read_delim(Intron.bed, delim='\t', col_names=c("Chrom", "Start", "Stop", "Name", "Score", "Strand")) %>%
    mutate(transcript = str_replace(Name, "^(.+?)_.+$", "\\1"))

ExpressedGenesList <- read_delim(ExpressedGenesList, delim='\t',
                                col_names = c(c("Chrom", "Start", "Stop", "Name", "Score", "Strand")))

Intron.bed.merged <- merge(Intron.bed,GeneAnnotation,by="transcript") %>%
   filter(gene %in% ExpressedGenesList$Name) %>%
    dplyr::select(-transcript) %>%
    write_delim(Output, col_names = F , delim = '\t')


# Vestigial, logRPKM calculation in chRNA-seq expression

# library(edgeR)

#CountTable  <- read_delim("featureCounts/chRNA.Expression/Counts.txt", delim='\t', comment='#')

# LogRPKM <- CountTable %>%
#   dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
#   column_to_rownames("Geneid") %>%
#   DGEList() %>%
#   rpkm(gene.length=CountTable$Length, prior.count=0.01, log=T)


# LogRPKM %>%
#   as.data.frame() %>%
#   rownames_to_column("Gene") %>%
#   filter(Gene %in% ProteinCodingGenes$gene) %>%
#   mutate(Mean=apply(dplyr::select(., -Gene),1,mean)) %>%
#   # filter(Mean>-5) %>% dim()
#   # pull(Mean) %>% hist()
#   arrange(desc(Mean)) %>% 
#   # head(10000) %>%
#   write_delim("featureCounts/chRNA.Expression/CountTable.MeanLogRPKM.txt.gz", delim='\t')

