library(tidyverse)

ThousandGenomesProject.Samples <- read_tsv("config/ExternalFastqDataAccessions/igsr_samples_Phase3ReleaseSampleList.tsv")

Geuvadis <- read_delim("config/ExternalFastqDataAccessions/Geuvadis_RNA-seq.tsv", delim='\t') %>%
    select(IndID=`Source Name`, fastq_ftp=`Comment[FASTQ_URI]`) %>%
    mutate(Read=str_replace(fastq_ftp, ".+_([12]).fastq.gz", "R\\1_ftp")) %>%
    spread(key="Read", value="fastq_ftp") %>%
    mutate(Assay="RNA-seq", Phenotype="Expression.Splicing", study_accession="E-GEUV-1")

#Samples that can be downloaded from ENA project table report online
GrubertEtAl_ChIPSeq <- read_tsv("config/ExternalFastqDataAccessions/GrubertSampleList.tsv")  %>%
    select(study_accession, fastq_ftp, fastq_aspera, sample_title) %>%
    separate(sample_title, into=c("IndID", "Phenotype"), sep='_', remove=F) %>%
    filter(Phenotype %in% c("H3K27AC", "H3K4ME1", "H3K4ME3")) %>%
    mutate(IndID = str_replace(IndID, "-\\d+", "")) %>%
    mutate(Assay="ChIP-seq") %>%
    separate(fastq_ftp, into=c("R1_ftp", "R2_ftp"), sep=';') %>%
    separate(fastq_aspera, into=c("R1_aspera", "R2_aspera"), sep=';') %>%
    rename(notes=sample_title) %>%
    mutate(Include=T)

DingEtAl_CTCF <- read_tsv("config/ExternalFastqDataAccessions/DingSampleList.tsv") %>%
    select(study_accession, fastq_ftp, fastq_aspera, sample_alias) %>%
    mutate(IndID = str_replace(sample_alias, "(GM\\w+?)[-_].+$", "\\1")) %>%
    mutate(Assay="ChIP-seq", Phenotype="CTCF") %>%
    separate(fastq_ftp, into=c("R1_ftp", "R2_ftp"), sep=';') %>%
    separate(fastq_aspera, into=c("R1_aspera", "R2_aspera"), sep=';') %>%
    rename(notes=sample_alias)

#Only Rep1 samples will be used in QTL calling
Table.out <-
    bind_rows(Geuvadis, GrubertEtAl_ChIPSeq, DingEtAl_CTCF) %>%
    mutate(IndID=case_when(
                           IndID %in% ThousandGenomesProject.Samples$`Sample name` ~ IndID,
                           TRUE ~ str_replace(IndID, "GM(\\d+)\\w*?$", "NA\\1")
                           )) %>%
    mutate(IsIn1KG.Phase3=IndID %in% ThousandGenomesProject.Samples$`Sample name`) %>%
    mutate(R1_local=NA, R2_local=NA) %>%
    mutate(Include=T) %>%
    filter((!is.na(R1_ftp) & !is.na(R2_ftp)) | (!is.na(R1_local) & !is.na(R2_local)) | (!is.na(R1_aspera) & !is.na(R2_aspera))) %>%
    group_by(IndID, Phenotype, Assay, study_accession) %>%
    mutate(RepNumber = row_number()) %>%
    ungroup() %>%
    mutate(Include = case_when(
                               RepNumber == 1 ~ T,
                               TRUE ~ F
                               )) %>%
    select(IndID, Assay, Phenotype, study_accession, R1_ftp, R2_ftp, R1_aspera, R2_aspera, R1_local, R2_local, notes, IsIn1KG.Phase3, Include, RepNumber)

# write_tsv(Table.out, "config/samples.tsv")

####
#Append rows for homemade samples

Sequenced.202104 <- data.frame(path=Sys.glob("/cds/yangili1/bjf79/Fastq/20210427_NovaSeq_chRNASeq/*/FastQ/*.fastq.gz")) %>%
    filter(str_detect(path, "-\\w+-chRNA-")) %>%
    mutate(Lane_Line_Read = str_replace(path, ".+?/(21042.+?)/.+?-(\\w+)-chRNA-.+?_(R[12])_.+$", "\\1.\\2 \\3")) %>%
    separate(Lane_Line_Read, into=c("Lane.IndID", "Read"), sep=" ") %>%
    pivot_wider(names_from=Read, values_from=path) %>%
    separate(Lane.IndID, into=c("Lane", "IndID"), sep="\\.") %>%
    mutate(RepNumber=case_when(
                               str_detect(IndID, "A") ~ 1,
                               str_detect(IndID, "B") ~ 2,
                               TRUE ~ 1
                               )) %>%
    mutate(IndID=str_replace(IndID, "[AB]", "")) %>%
    mutate(IndID=paste0("NA", IndID)) %>%
    mutate(IsIn1KG.Phase3=IndID %in% ThousandGenomesProject.Samples$`Sample name`) %>%
    mutate(Include = case_when(
                               RepNumber == 1 ~ T,
                               TRUE ~ F
                               )) %>%
    rename(R1_local=R1, R2_local=R2) %>%
    select(-Lane) %>%
    mutate(Assay="RNA-seq", Phenotype="chRNA.Expression.Splicing")


Table.out.w.InHouseSamples <- bind_rows(Table.out, Sequenced.202104)
write_tsv(Table.out.w.InHouseSamples, "config/samples.tsv")

### Add in house cut and tag tests
data.frame(path=Sys.glob("/cds/yangili1/bjf79/Fastq/20210427_NovaSeq_chRNASeq/*/FastQ/*BF-[HL]*.fastq.gz")) %>%
    mutate(Lane_Antibody_Amount_Read = str_replace(path, ".+?/(21042.+?)/.+?-BF-[HL]-(.+?)-(.+?)_S\\d+_(R[12])_.+$", "\\1.NA19210 \\2 \\3 \\4")) %>%
    separate(Lane_Antibody_Amount_Read, into=c("Lane.IndID", "Phenotype","notes", "Read"), sep=" ") %>%
    pivot_wider(names_from=Read, values_from=path) %>%
    separate(Lane.IndID, into=c("Lane", "IndID"), sep="\\.") %>%
    mutate(RepNumber = 1) %>%
    rename(R1_local=R1, R2_local=R2) %>%
    select(-Lane) %>%
    mutate(Phenotype=toupper(Phenotype), Assay="CutAndTag", Include=F) %>%
    bind_rows(Table.out.w.InHouseSamples, .) %>%
    write_tsv("config/samples.tsv")
