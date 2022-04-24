#!/usr/bin/env sh

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : Download
# @created     : Monday Apr 19, 2021 21:26:18 CDT
#
# @description : Download files w/ accession numbers
######################################################################

curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA268086&result=read_run&fields=study_accession,sample_accession,run_accession,library_name,nominal_length,library_layout,library_selection,fastq_ftp,fastq_aspera,sample_alias,sample_title&format=tsv&download=true" > GrubertSampleList.tsv

curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB1350&result=read_run&fields=study_accession,sample_accession,run_accession,library_name,nominal_length,library_layout,library_selection,fastq_ftp,fastq_aspera,sample_alias,sample_title&format=tsv&download=true" > DingSampleList.tsv

curl "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" > Geuvadis_RNA-seq.tsv

