library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
in_fn <- args[1]
out_fn <- args[2]

in_fn
out_fn

# in_fn <- "ReferenceGenome/Annotations/GTFTools_BasicAnnotations/gencode.v34.chromasomal.exons.bed"
# out_fn <- "scratch/internalexons.bed"

Exons <- read_tsv(in_fn, col_names=c("chr", "start", "stop" ,"strand", "transcript", "gene")) %>%
    distinct(chr, start, stop, .keep_all=T)

ExonsPlus <- Exons %>%
  filter(strand=="+") %>%
  arrange(gene, start) %>%
  group_by(gene) %>%
  mutate(ExonNumber=rank(start, ties.method = "min"))

ExonsMinus <- Exons %>%
  filter(strand=="-") %>%
  arrange(gene, desc(stop)) %>%
  group_by(gene) %>%
  mutate(ExonNumber=rank(desc(stop), ties.method = "min"))

dplyr::bind_rows(ExonsPlus, ExonsMinus) %>% head()

Out.df <- dplyr::bind_rows(ExonsPlus, ExonsMinus) %>%
  add_count(gene, name="ExonCountStart") %>%
  filter(ExonCountStart>=3) %>%
  group_by(gene) %>%
  dplyr::slice(2:n()) %>%
  dplyr::slice(1:(n()-1)) %>%
  ungroup() %>%
    # ggplot(aes(x=ExonNumber)) +
    # geom_histogram() +
    # xlim(c(0,10))
  unite(name, gene, ExonNumber) %>%
  mutate(chr = paste0("chr", chr)) %>%
  mutate(score = "0") %>%
  arrange(chr, start, stop) %>%
  distinct(.keep_all = T) %>%
  dplyr::select(chr, start, stop, name, score, strand)

head(Out.df)

write_delim(Out.df, file=out_fn, delim='\t', col_names = FALSE)

