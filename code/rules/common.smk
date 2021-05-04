import pandas as pd
import os

autosomes = [str(i) for i in range(1,23)]

# ## GEUVADIS RNA-seq info
# GEUVADIS_samples = pd.read_csv("../data/E-GEUV-1.sdrf.txt", sep='\t')
# GEUVADIS_sample_links_dict = dict(zip(GEUVADIS_samples["Scan Name"], GEUVADIS_samples["Comment[FASTQ_URI]"]))
# GEUVADIS_line_fastq_dict = dict([(k.split(".")[0], k[:-10]) for k in GEUVADIS_sample_links_dict.keys()])

# ## Grubert et al ChIP-seq info
# Grubert_samples = pd.read_csv("../data/PRJNA268086_SraRunTable.GrubertEtAl.csv")
# # df containing only antibody, SRR accession, and cell line
# Grubert_ChIP_seq = Grubert_samples[Grubert_samples['Assay Type']=="ChIP-Seq"][['Run', 'Antibody', 'Cell_Line']]
# Grubert_ChIP_seq_dict = dict(zip(zip(Grubert_ChIP_seq['Antibody'], Grubert_ChIP_seq['Cell_Line']), Grubert_ChIP_seq['Run']))

## All Fastq samples
Fastq_samples = pd.read_csv("config/samples.tsv", sep='\t')
# Retrieve samples where
# Fastq_samples.loc[(Fastq_samples['IndID'] == "NA18489") & (Fastq_samples['Phenotype'] == "H3K27AC") & (Fastq_samples['Include']==True)]
# Get unique samples by group
# Fastq_samples.groupby(['IndID', 'Phenotype']).size()
#Get all unique IndID Phenotype groups that are included:
#Fastq_samples.loc[(Fastq_samples['Include']==True) ,['IndID', 'Phenotype'] ].drop_duplicates
# If any lines are identical for IndID, Phenotype, and RepNumber, combine fastq

def require_at_least_one(filelist):
    """
    https://stackoverflow.com/questions/65227729/how-to-make-snakemake-input-optional-but-not-empty
    """
    try:
        existing = [file for file in filelist if os.path.isfile(file)]
        return existing if len(existing) else "non_existing_file"
    except:
        return "non_existing_file"

#Input functions

def GetFastqLocal(wildcards):
    """
    Return list of local fastq filepaths defined in Fastq_samples df if they exist.
    """
    df_subset = Fastq_samples.loc[
            (Fastq_samples['IndID'] == wildcards.IndID) &
            (Fastq_samples['Phenotype'] == wildcards.Phenotype) &
            (Fastq_samples['RepNumber'] == wildcards.Rep)]
    return require_at_least_one(df_subset[['R1_local', 'R2_local']][0].values())


def GetFastqLocalFuncs(Read):
    def F(wildcards):
        df_subset = Fastq_samples.loc[
                (Fastq_samples['IndID'] == wildcards.IndID) &
                (Fastq_samples['Phenotype'] == wildcards.Phenotype) &
                (Fastq_samples['RepNumber'] == int(wildcards.Rep))]
        return require_at_least_one(df_subset[Read].tolist())
    return F

def GetDownloadLinkFuncs(ColName):
    def F(wildcards):
        df_subset = Fastq_samples.loc[
                (Fastq_samples['IndID'] == wildcards.IndID) &
                (Fastq_samples['Phenotype'] == wildcards.Phenotype) &
                (Fastq_samples['RepNumber'] == int(wildcards.Rep))]
        return df_subset[ColName].fillna('').tolist()
    return F

def GetFastpParams(wildcards):
    if wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "--umi --umi_loc read1 --umi_len 12"

def GetReadsForAlignmentFuncs(Read):
    """
    Not used
    Rules to align reads need cutadapt reads, depending on the phenotype
    """
    def F(wildcards):
        if wildcards.Phenotype in ['CTCF']:
            return "FastqCutadapt/{wildcards.Phenotype}/{wildcards.IndID}/{wildcards.Rep}.{Read}.fastq.gz".format(wildcards=wildcards, Read=Read)
        else:
            return "Fastq/{wildcards.Phenotype}/{wildcards.IndID}/{wildcards.Rep}.{Read}.fastq.gz".format(wildcards=wildcards, Read=Read)
    return F

