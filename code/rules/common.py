import pandas as pd
import numpy as np
import os
import glob

autosomes = [str(i) for i in range(1,23)]

N_PermutationChunks = 50
MyPhenotypes = ["chRNA.IR", "Expression.Splicing", "chRNA.Expression.Splicing",  "H3K27AC", "CTCF", "H3K4ME3", "chRNA.Splicing", "polyA.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "Expression.Splicing.Subset_YRI", "polyA.Splicing.Subset_YRI"]
PhenotypesToColoc = [p for p in MyPhenotypes if p not in ["chRNA.Splicing", "polyA.Splicing"]]

## All Fastq samples
Fastq_samples = pd.read_csv("config/samples.tsv", sep='\t', comment='#')

# Define some df subsets from the samples df as useful global variables
PhenotypeSet = Fastq_samples['Phenotype'].unique().tolist() + ["chRNA.IR",  "polyA.Splicing", "chRNA.Expression"]
ChromatinProfilingPhenotypes = Fastq_samples.loc[ (Fastq_samples['Assay'].isin(["ChIP-seq", "CutAndTag"]))  ]['Phenotype'].unique().tolist()
RNASeqPhenotypes = Fastq_samples.loc[ (Fastq_samples['Assay']=="RNA-seq")  ]['Phenotype'].unique().tolist()
RNASeqExpressionPhenotypes = ['polyA.Expression', 'chRNA.Expression', 'Expression.Splicing.Subset_YRI']
RNASeqPhenotypes_extended = RNASeqPhenotypes + RNASeqExpressionPhenotypes
ChromatinProfilingSamples_df = Fastq_samples.loc[ (Fastq_samples['Assay'].isin(["ChIP-seq", "CutAndTag"])) , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()
RNASeqSamples_df = Fastq_samples.loc[ (Fastq_samples['Assay']=="RNA-seq") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()

#Other useful lists
AllChromatinProfilingBams = expand("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam", zip, Phenotype= ChromatinProfilingSamples_df['Phenotype'], IndID=ChromatinProfilingSamples_df['IndID'], Rep=ChromatinProfilingSamples_df['RepNumber'])
AllRNASeqBams = expand("Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam", zip, Phenotype=RNASeqSamples_df['Phenotype'], IndID=RNASeqSamples_df['IndID'], Rep=RNASeqSamples_df['RepNumber'])
AllBams = AllChromatinProfilingBams + AllRNASeqBams
AllBais = [fn + ".bai" for fn in AllChromatinProfilingBams]
chRNASeqSamples_df = Fastq_samples.loc[ (Fastq_samples['Phenotype']=="chRNA.Expression.Splicing") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()

# Get random sample for each phenotype
Fastq_samples.groupby('Phenotype').apply(lambda x: x.sample(2, random_state=1)).reset_index(drop=True)

# Read in table of GWAS summary stat links (obtained from GWAS catalog's online
# interface)
gwas_df = pd.read_csv("../data/list_gwas_summary_statistics_PMID27863252.csv", index_col='Study accession')


# print(ChromatinProfilingPhenotypes)

# Retrieve samples where
# Fastq_samples.loc[(Fastq_samples['IndID'] == "NA18489") & (Fastq_samples['Phenotype'] == "H3K27AC") & (Fastq_samples['Include']==True)]
# Get unique samples by group
# Fastq_samples.groupby(['IndID', 'Phenotype']).size()
#Get all unique IndID Phenotype groups that are included:
#Fastq_samples.loc[(Fastq_samples['Include']==True) ,['IndID', 'Phenotype'] ].drop_duplicates
# If any lines are identical for IndID, Phenotype, and RepNumber, combine fastq

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 40000

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

def GetDownloadLinkFuncs(LinkType):
    def F(wildcards):
        if wildcards.Read == "SE":
            wildcards.Read = "R1"
        df_subset = Fastq_samples.loc[
                (Fastq_samples['IndID'] == wildcards.IndID) &
                (Fastq_samples['Phenotype'] == wildcards.Phenotype) &
                (Fastq_samples['RepNumber'] == int(wildcards.Rep))]
        return df_subset[wildcards.Read + '_' + LinkType].fillna('').tolist()
    return F

def GetFastpParamsUmi(wildcards):
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

def GetBamForBigwig(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
        return "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam"


def GetBaiForBigwig(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
        return "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam.bai"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam.bai"

def GetUnfilteredBamForBigwig(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
        return "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam"


def GetUnfilteredBaiForBigwig(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
        return "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam.bai"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai"


def GetBigwigParams(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
        return "-split"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return ""
        # return "-pc"
    else:
        return ""

def GetLibStrandForRegtools(wildcards):
    if wildcards.Phenotype == "Expression.Splicing":
        return 0
    elif wildcards.Phenotype == "chRNA.Expression.Splicing":
        return 2

def GetQualimapLibtype(wildcards):
    if wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "-p strand-specific-forward"

def GetAnnotationsForPhenotype(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes_extended:
        return "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf"
    elif wildcards.Phenotype in ["H3K27AC", "H3K4ME1", "CTCF"]:
        return "PeakCalling/{Phenotype}_peaks.narrowPeak.saf"
    elif wildcards.Phenotype in ["H3K4ME3", "H3K36ME3", "POL2S2", "POL2S5", "H3K9ME3", "H3K79ME2"]:
        return "PeakCalling/{Phenotype}_peaks.broadPeak.saf"

def GetFeatureCountsParams(wildcards):
    if wildcards.Phenotype == "chRNA.Expression":
        return "-s 1"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return "-F SAF"
    else:
        return ""


def GetBamForPhenotype(wildcards):
    df_subset = Fastq_samples.loc[
                (Fastq_samples['Phenotype'] == wildcards.Phenotype)]
    if wildcards.Phenotype in RNASeqPhenotypes:
        return expand("Alignments/STAR_Align/{{Phenotype}}/{IndID}/{Rep}/Filtered.bam", zip, IndID=df_subset['IndID'], Rep=df_subset['RepNumber'])
    elif wildcards.Phenotype == "polyA.Expression":
        df_subset = Fastq_samples.loc[
                    (Fastq_samples['Phenotype'] == "Expression.Splicing")]
        return expand("Alignments/STAR_Align/Expression.Splicing/{IndID}/{Rep}/Filtered.bam", zip, IndID=df_subset['IndID'], Rep=df_subset['RepNumber'])
    elif wildcards.Phenotype == "chRNA.Expression":
        df_subset = Fastq_samples.loc[
                    (Fastq_samples['Phenotype'] == "chRNA.Expression.Splicing")]
        return expand("Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/{Rep}/Filtered.bam", zip, IndID=df_subset['IndID'], Rep=df_subset['RepNumber'])
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return expand("Alignments/Hisat2_Align/{{Phenotype}}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam", zip, IndID=df_subset['IndID'], Rep=df_subset['RepNumber'])

def GetAnnotationsForRegion(wildcards):
    if wildcards.Region == "AtTSS":
        return "ReferenceGenome/Annotations/TSSRegions.saf"
