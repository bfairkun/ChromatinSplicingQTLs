library(tidyverse)

ThousandGenomesProject.Samples <- read_tsv("config/ExternalFastqDataAccessions/igsr_samples_Phase3ReleaseSampleList.tsv")

Geuvadis <- read_delim("config/ExternalFastqDataAccessions/Geuvadis_RNA-seq.tsv", delim='\t') %>%
    select(IndID=`Source Name`, fastq_ftp=`Comment[FASTQ_URI]`) %>%
    mutate(Read=str_replace(fastq_ftp, ".+_([12]).fastq.gz", "R\\1_ftp")) %>%
    spread(key="Read", value="fastq_ftp") %>%
    mutate(Assay="RNA-seq", Phenotype="Expression.Splicing", study_accession="E-GEUV-1")

#Samples that can be downloaded from ENA project table report online
GrubertEtAl_ChIPSeq <- read_tsv("config/ExternalFastqDataAccessions/GrubertSampleList.tsv")  %>%
    select(study_accession, fastq_ftp, sample_title) %>%
    separate(sample_title, into=c("IndID", "Phenotype"), sep='_', remove=F) %>%
    filter(Phenotype %in% c("H3K27AC", "H3K4ME1", "H3K4ME3")) %>%
    mutate(IndID = str_replace(IndID, "-\\d+", "")) %>%
    mutate(Assay="ChIP-seq") %>%
    separate(fastq_ftp, into=c("R1_ftp", "R2_ftp"), sep=';') %>%
    rename(notes=sample_title) %>%
    mutate(Include=)

DingEtAl_CTCF <- read_tsv("config/ExternalFastqDataAccessions/DingSampleList.tsv") %>%
    select(study_accession, fastq_ftp, sample_alias) %>%
    mutate(IndID = str_replace(sample_alias, "(GM\\w+?)[-_].+$", "\\1")) %>%
    mutate(Assay="ChIP-seq", Phenotype="CTCF") %>%
    separate(fastq_ftp, into=c("R1_ftp", "R2_ftp"), sep=';') %>%
    rename(notes=sample_alias)

Table.out <-
    bind_rows(Geuvadis, GrubertEtAl_ChIPSeq, DingEtAl_CTCF) %>%
    mutate(IndID=case_when(
                           IndID %in% ThousandGenomesProject.Samples$`Sample name` ~ IndID,
                           TRUE ~ str_replace(IndID, "GM(\\d+)\\w*?$", "NA\\1")
                           )) %>%
    mutate(IsIn1KG.Phase3=IndID %in% ThousandGenomesProject.Samples$`Sample name`) %>%
    mutate(R1_local=NA, R2_local=NA) %>%
    mutate(Include=T) %>%
    group_by(IndID, Phenotype, Assay, study_accession) %>%
    mutate(RepNumber = row_number()) %>%
    ungroup() %>%
    mutate(Include = case_when(
                               RepNumber == 1 ~ T,
                               TRUE ~ F
                               )) %>%
    select(IndID, Assay, Phenotype, study_accession, R1_ftp, R2_ftp, R1_local, R2_local, notes, IsIn1KG.Phase3, Include, RepNumber)

write_tsv(Table.out, "config/samples.tsv")



