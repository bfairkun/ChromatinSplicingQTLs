import pandas as pd
import os

###### Config file and sample sheets #####
configfile: "config.yaml"

autosomes = [str(i) for i in range(1,23)]

## GEUVADIS RNA-seq info
GEUVADIS_samples = pd.read_csv("../data/E-GEUV-1.sdrf.txt", sep='\t')
GEUVADIS_sample_links_dict = dict(zip(GEUVADIS_samples["Scan Name"], GEUVADIS_samples["Comment[FASTQ_URI]"]))
GEUVADIS_line_fastq_dict = dict([(k.split(".")[0], k[:-10]) for k in GEUVADIS_sample_links_dict.keys()])

## Grubert et al ChIP-seq info
Grubert_samples = pd.read_csv("../data/PRJNA268086_SraRunTable.GrubertEtAl.csv")
# df containing only antibody, SRR accession, and cell line
Grubert_ChIP_seq = Grubert_samples[Grubert_samples['Assay Type']=="ChIP-Seq"][['Run', 'Antibody', 'Cell_Line']]
Grubert_ChIP_seq_dict = dict(zip(zip(Grubert_ChIP_seq['Antibody'], Grubert_ChIP_seq['Cell_Line']), Grubert_ChIP_seq['Run']))

