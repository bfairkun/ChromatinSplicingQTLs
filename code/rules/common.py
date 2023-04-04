import pandas as pd
import numpy as np
import os
import glob

autosomes = [str(i) for i in range(1,23)]

N_PermutationChunks = 100
NumPvalsForPi1Chunks = 40

ChunkNumbers = range(0, 1+N_PermutationChunks) 

MyPhenotypes = ["Expression.Splicing", "Expression.Splicing.Subset_YRI", "chRNA.Expression.Splicing", 
                "MetabolicLabelled.30min", "MetabolicLabelled.60min",
                "CTCF", "H3K27AC", "H3K4ME3", "H3K4ME1", "H3K36ME3", "H3K36ME3_ncRNA", "ProCap",
                "polyA.Splicing", "polyA.Splicing.Subset_YRI", "chRNA.Splicing", 
                'MetabolicLabelled.30min.Splicing', 'MetabolicLabelled.60min.Splicing', 
                #"polyA.Expression_ncRNA", "polyA.Expression_ncRNA.Subset_YRI", 
                "chRNA.Expression_ncRNA", 
                #'MetabolicLabelled.30min_ncRNA', 'MetabolicLabelled.60min_ncRNA',
                'APA_Nuclear', 'APA_Total',
                'polyA.IER', 'polyA.IER.Subset_YRI', 'chRNA.IER',
                'MetabolicLabelled.30min.IER', 'MetabolicLabelled.60min.IER',
                "chRNA.Slopes", 'chRNA.Splicing.Order', "DNaseISensitivity"]

#, 'chRNA.RNA.Editing']
# "polyA.Expression.AllRNA.Subset_YRI", "MetabolicLabelled.30min.AllRNA.Subset_YRI", 
#                 "MetabolicLabelled.60min.AllRNA.Subset_YRI", "chRNA.Expression.AllRNA.Subset_YRI"]

######## We can use this later to add ncRNAs to QTL analysis
#MyPhenotypes += ncRNA_Phenotypes

#PhenotypesToColoc = [p for p in MyPhenotypes if ((p not in ["Expression.Splicing", "polyA.Splicing", "polyA.IR"]) and (p not in ncRNA_Phenotypes))]

# PhenotypesToColoc = MyPhenotypes #[p for p in MyPhenotypes if (p not in ["Expression.Splicing", "polyA.Splicing", "polyA.IR"])]
PhenotypesToColoc = ["Expression.Splicing", "Expression.Splicing.Subset_YRI", "chRNA.Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "H3K27AC", "H3K4ME3", "H3K4ME1", "H3K36ME3",'APA_Nuclear', 'APA_Total', "polyA.Splicing", "polyA.Splicing.Subset_YRI", "chRNA.Splicing", "ProCap"]

## All Fastq samples
Fastq_samples = pd.read_csv("config/samples.tsv", sep='\t', comment='#')

## Small Molecule treatment RNA seq samples
SM_samples = pd.read_csv("config/SmallMoleculeRNASeq.Samples.tsv", sep='\t')

# Define some df subsets from the samples df as useful global variables
PhenotypeSet = Fastq_samples['Phenotype'].unique().tolist() + ["chRNA.IR",  "polyA.Splicing", "chRNA.Expression"]
ChromatinProfilingPhenotypes = Fastq_samples.loc[ (Fastq_samples['Assay'].isin(["ChIP-seq", "CutAndTag", "DNase-seq"]))  ]['Phenotype'].unique().tolist()
RNASeqPhenotypes = Fastq_samples.loc[ (Fastq_samples['Assay']=="RNA-seq")  ]['Phenotype'].unique().tolist()
RNASeqExpressionPhenotypes = ['polyA.Expression', 'chRNA.Expression', 'Expression.Splicing.Subset_YRI']
RNASeqPhenotypes_extended = RNASeqPhenotypes + RNASeqExpressionPhenotypes
ChromatinProfilingSamples_df = Fastq_samples.loc[ (Fastq_samples['Assay'].isin(["ChIP-seq", "CutAndTag", "DNase-seq"])) , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()
RNASeqSamples_df = Fastq_samples.loc[ (Fastq_samples['Assay']=="RNA-seq") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()
RNASeqSamplesNoProcap_df = Fastq_samples.loc[ (Fastq_samples['Assay']=="RNA-seq") & (Fastq_samples["Phenotype"] != "ProCap") , ['Phenotype', 'IndID', 'RepNumber']  ].drop_duplicates()
ProCapSamples_df = Fastq_samples.loc[ (Fastq_samples['Phenotype']=="ProCap") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()

DNaseSamples_df = Fastq_samples.loc[ (Fastq_samples['Phenotype']=="DNaseISensitivity") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()

#Other useful lists
AllChromatinProfilingBams = expand("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam", zip, Phenotype= ChromatinProfilingSamples_df['Phenotype'], IndID=ChromatinProfilingSamples_df['IndID'], Rep=ChromatinProfilingSamples_df['RepNumber'])
AllRNASeqBams = expand("Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam", zip, Phenotype=RNASeqSamples_df['Phenotype'], IndID=RNASeqSamples_df['IndID'], Rep=RNASeqSamples_df['RepNumber'])
AllBams = AllChromatinProfilingBams + AllRNASeqBams
AllBais = [fn + ".bai" for fn in AllChromatinProfilingBams]
chRNASeqSamples_df = Fastq_samples.loc[ (Fastq_samples['Phenotype']=="chRNA.Expression.Splicing") , ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()


igsr_samples = pd.read_csv('../data/igsr_samples.tsv.gz', sep='\t', index_col=0)
chRNASeqSamples = list(pd.Index(chRNASeqSamples_df['IndID'].unique()).intersection(igsr_samples.index))

# Get random sample for each phenotype
Fastq_samples.groupby('Phenotype').apply(lambda x: x.sample(2, random_state=1)).reset_index(drop=True)

gwas_df = pd.read_csv("config/gwas_table.tsv", index_col='gwas', sep='\t', usecols=['gwas', 'trait', 'FTPPath', "SummaryStatsLocalFilepath", 'ProcessingMethod', "Continuous"])

colocs_df = pd.read_csv("config/ColocRunWildcards.tsv", index_col=0, comment='#', sep='\t', keep_default_na=False)
colocs_gwas = colocs_df.loc[colocs_df['FeatureCoordinatesRedefinedFor']=='ForGWASColoc']
colocs_genewise = colocs_df.loc[colocs_df['FeatureCoordinatesRedefinedFor']=='ForColoc']


chRNAsamples = list(Fastq_samples.loc[ (Fastq_samples['Phenotype']=='chRNA.Expression.Splicing') ].drop_duplicates().IndID)
polyAsamples = list(Fastq_samples.loc[ (Fastq_samples['Phenotype']=='Expression.Splicing') ].drop_duplicates().IndID)
metabolic30samples = list(Fastq_samples.loc[ (Fastq_samples['Phenotype']=='MetabolicLabelled.30min') ].drop_duplicates().IndID)
metabolic60samples = list(Fastq_samples.loc[ (Fastq_samples['Phenotype']=='MetabolicLabelled.60min') ].drop_duplicates().IndID)



proseq_samples = ['GM18505', 'GM19239', 'GM19238', 'GM19222', 'GM19193', 'GM19131', 'GM19099', 'GM18522', 'GM18520', 'GM18517']
chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
        'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
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
        return 52000


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
    if wildcards.Phenotype in [i for i in RNASeqPhenotypes if i != "ProCap"]:
        return "-split"
#     elif wildcards.Phenotype == "ProCap":
#         return "-5"
    elif wildcards.Phenotype in ChromatinProfilingPhenotypes:
        return ""
    else:
        return ""

def GetLibStrandForRegtools(wildcards):
    if wildcards.Phenotype == "Expression.Splicing":
        return 0
    elif wildcards.Phenotype == "chRNA.Expression.Splicing":
        return 2

def GetQualimapLibtype(wildcards):
    if wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "-p strand-specific-reverse"

def GetAnnotationsForPhenotype(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes_extended:
        return "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf"
    elif wildcards.Phenotype in ["H3K27AC", "H3K4ME1", "CTCF", "DNaseISensitivity"]:
        return "PeakCalling/{Phenotype}_peaks.narrowPeak.saf"
    elif wildcards.Phenotype in ["H3K4ME3", "H3K36ME3", "POL2S2", "POL2S5", "H3K9ME3", "H3K79ME2"]:
        return "PeakCalling/{Phenotype}_peaks.broadPeak.saf"

def PairedEndParams(wildcards):
    if wildcards.Phenotype in ["MetabolicLabelled.30min", "MetabolicLabelled.60min", 
                               "ProCap", "DNaseISensitivity"]:
        return ""
    else:
        return "-p"

def GetFeatureCountsParams(wildcards):
    try:
        if (wildcards.Phenotype=="ProCap"): 
            return  "-s 1 -F SAF --read2pos 5"
        if wildcards.Phenotype == "chRNA.Expression":
            if  wildcards.Region == "AtInternalExons":
                return "-s 2 -F SAF"
            else:
                return "-s 2"
        elif (wildcards.Phenotype in ChromatinProfilingPhenotypes) or (wildcards.Region == 'AtTSS'):
            return "-F SAF"
        elif wildcards.Phenotype == "polyA.Expression" and wildcards.Region == "AtInternalExons":
            return "-F SAF"
        else:
            return ""
    except AttributeError:
        return ""

def GetFeatureCountsParams2(wildcards):
    if wildcards.Phenotype == "chRNA.Expression":
        return "-s 2"
    else:
        return ""

def GetSTARJunctionScoreParams(wildcards):
    if wildcards.Phenotype == "ProCap":
        return "--scoreGap -1000000"
    else:
        return ""

def GetBamForPhenotype(wildcards):
    df_subset = Fastq_samples.loc[
                (Fastq_samples['Phenotype'] == wildcards.Phenotype)]
    if "MetabolicLabelled" in wildcards.Phenotype:
        return expand("Alignments/STAR_Align/{{Phenotype}}/{IndID}/1/Filtered.bam", IndID=df_subset['IndID'].unique())
    elif wildcards.Phenotype in RNASeqPhenotypes:
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
        if wildcards.Phenotype == 'DNaseISensitivity':
            return expand("Alignments/Hisat2_Align/DNaseISensitivity/{IndID}.merged.wasp_filterd.markdup.sorted.bam", zip, IndID=df_subset['IndID'].unique())
        else:
            return expand("Alignments/Hisat2_Align/{{Phenotype}}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam", zip, IndID=df_subset['IndID'], Rep=df_subset['RepNumber'])

def GetAnnotationsForRegion(wildcards):
    if wildcards.Region == "AtTSS":
        return "ReferenceGenome/Annotations/TSSRegions.saf"
    if wildcards.Region == "AtInternalExons":
        return "SplicingAnalysis/Annotations/InternalExons.noOverlapping.saf"
    if wildcards.Region == 'uaRNA_TSS':
        return "NonCodingRNA/annotation/uaRNA.TSS.saf"

def GetMolPhenotypesToColoc(wildcards):
    ProvidedMolPhenotypeList = colocs_df.loc[wildcards.ColocName]['MolPhenotypesToColoc']
    if ProvidedMolPhenotypeList == '':
        return ' '.join(PhenotypesToColoc)
    else:
        return ProvidedMolPhenotypeList


def GetColocTsvFormattedString(string):
    """
    return a snakemake input function that returns formatted string with wildcard values that match colocs_df based on ColocName index wildcard
    """
    def F(wildcards):
        return string.format(**colocs_df.loc[wildcards.ColocName].to_dict())
    return F

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 52000

def GetBigwigForDeeptToolscheRNA(wildcards):
    template = 'bigwigs/chRNA.Expression.Splicing_stranded/{IndID}.1.{strand}.bw'
    # needed for as long as strands are swapped
    if wildcards.strand == 'plus':
        strand = 'minus'
    else:
        strand = 'plus'
    bw = template.format(IndID = wildcards.IndID, strand = strand)
    return bw


long_read_samples = ['GM' + str(i) for i in range(1, 11)]
#long_read_samples += ['SRR1163655', 'SRR1163657', 'SRR1163658']

long_read_samples_df = pd.read_csv("config/samples.longreads.tsv", sep='\t', comment='#')

long_read_samples_to_download = list(long_read_samples_df.loc[ (long_read_samples_df['Assay']=="ONT") , ['Phenotype', 'IndID'] ].drop_duplicates().index)

long_read_samples_ONT = [row.Phenotype + '.' + row.IndID for idx, row in long_read_samples_df.iterrows()]

all_long_reads = long_read_samples + long_read_samples_ONT

#expand("{Phenotype}.{IndID}", zip, Phenotype=long_read_samples_df['Phenotype'], IndID=long_read_samples_df['IndID'])

# long_read_samples += long_read_samples_from_df

