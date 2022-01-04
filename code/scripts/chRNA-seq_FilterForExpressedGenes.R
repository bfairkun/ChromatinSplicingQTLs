library(tidyverse)
library(edgeR)
#library(gplots)

# """
# Rscript to filter introns bed for introns in expressed host genes
# """

CountTable  <- read_delim("featureCounts/chRNA.Expression.Splicing/Counts.txt", delim='\t', comment='#')# %>%
#    rename_if(str_detect(names(.), "/"), funs(str_replace(., ".+/(.+)$", "\\1")))

Gtf <- read_delim("ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf", comment='#', delim='\t', col_names=c("Chr", "Source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

# ProteinCodingGenes <- Gtf %>%
#     filter(str_detect(attribute, 'gene_type "protein_coding"')) %>%
#     mutate(attribute = str_replace(attribute, 'gene_id "(.+?)".+', "\\1")) %>% distinct(attribute) %>% pull(attribute)

ProteinCodingGenes <-
    Gtf %>%
    filter(feature=="transcript") %>%
#     filter(str_detect(attribute, 'gene_type "protein_coding"')) %>%
    mutate(gene = str_replace(attribute, 'gene_id "(.+?)".+', "\\1")) %>%
    mutate(transcript = str_replace(attribute, '.+?transcript_id "(.+?)".+$', "\\1")) %>%
    dplyr::select(gene, transcript)

LogRPKM <- CountTable %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  rpkm(gene.length=CountTable$Length, prior.count=0.01, log=T)


#LogRPKM %>% as.data.frame() %>% cor() %>% heatmap.2(trace="none")

#LogRPKM %>% as.data.frame() %>% apply(1, mean) %>% hist()

# Get top 10000 protein coding genes by expression.
# Write out average LogRPKM of protein coding genes
LogRPKM %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% ProteinCodingGenes$gene) %>%
  mutate(Mean=apply(dplyr::select(., -Gene),1,mean)) %>%
  # filter(Mean>-5) %>% dim()
  # pull(Mean) %>% hist()
  arrange(desc(Mean)) %>% 
  # head(10000) %>%
  write_delim("featureCounts/chRNA.Expression.Splicing/CountTable.MeanLogRPKM.txt.gz", delim='\t')

#Read in intron bed
Intron.bed  <- read_delim("ReferenceGenome/Annotations/Introns.GencodeV34.hg38.UCSC.bed.gz", delim='\t', col_names=c("Chrom", "Start", "Stop", "Name", "Score", "Strand")) %>%
    mutate(transcript = str_replace(Name, "^(.+?)_.+$", "\\1"))

# Filter introns bed for top 10K expressed host genes
ExpressedGenes <- LogRPKM %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% ProteinCodingGenes$gene) %>%
  mutate(Mean=apply(dplyr::select(., -Gene),1,mean)) %>%
  arrange(desc(Mean)) %>% 
  head(10000) %>%
  dplyr::select(Gene, Mean) %>%
  inner_join(ProteinCodingGenes, by =c("Gene"="gene")) %>%
 pull(transcript)


Intron.bed.merged <- merge(Intron.bed,ProteinCodingGenes,by="transcript")

Intron.bed.merged %>%
    filter(transcript %in% ExpressedGenes) %>%
    dplyr::select(-transcript) %>%
    write_delim("Misc/GencodeHg38_all_introns.expressedHostGenes.bed.gz", col_names = F , delim = '\t')
