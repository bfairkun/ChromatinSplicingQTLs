rule GatherNormFactors:
    input:
        expand("featureCounts/{Phenotype}/NormFactors.tsv", Phenotype=["polyA.Expression", "chRNA.Expression"])

def GetGeneList(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes_extended:
        return "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    else:
        return []

def GetGeneListCol(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes_extended:
        return "4"
    else:
        return ""

rule CalculateNormFactorsForBigwig:
    input:
        Counts = "featureCounts/{Phenotype}/Counts.txt",
        OptionalGeneList = GetGeneList
    output:
        "featureCounts/{Phenotype}/NormFactors.tsv"
    log:
        "logs/CalculateNormFactorsForBigwig/{Phenotype}.log"
    params:
        OptionalGeneListCol = GetGeneListCol
    conda:
        "../envs/r_2.yaml"
    priority: 1
    shell:
        """
        Rscript scripts/CalculateNormFactorsForBigwig.R {input.Counts} {output} {input.OptionalGeneList} {params.OptionalGeneListCol} &> {log}
        """


def GetFeatureCountsNormFactorsFile(wildcards):
    if wildcards.Phenotype == "Expression.Splicing":
        return "featureCounts/polyA.Expression/NormFactors.tsv"
    elif wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "featureCounts/chRNA.Expression.Splicing/NormFactors.tsv"
    else:
        return "featureCounts/{Phenotype}/NormFactors.tsv"


rule MakeBigwigs_NormalizedToEdgeRFeatureCounts:
    """
    Scale bigwig to base coverage per billion chromosomal reads
    """
    input:
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = GetBamForBigwig,
        bai = GetBaiForBigwig,
        NormFactorsFile = GetFeatureCountsNormFactorsFile
    params:
        GenomeCovArgs=GetBigwigParams,
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
        BamToBigwigScript = "scripts/GenometracksByGenotype/BamToBigwig.sh"
        # Region = ""
    # wildcard_constraints:
    #     Phenotype = "|".join(RNASeqPhenotypes)
    shadow: "shallow"
    output:
        bw = "bigwigs/{Phenotype}/{IndID}.{Rep}.bw",
        bw_minus = []
    log:
        "logs/MakeBigwigs/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = much_more_mem_after_first_attempt
    shell:
        """
        ScaleFactor=$(bc <<< "scale=3;1000000000/$(grep '{input.bam}' {input.NormFactorsFile} | awk 'NR==1 {{print $2}}')")
        {params.BamToBigwigScript} {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs} -scale ${{ScaleFactor}}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """

use rule MakeBigwigs_NormalizedToEdgeRFeatureCounts as MakeBigwigs_NormalizedToEdgeRFeatureCounts_stranded with:
    output:
        bw = "bigwigs/{Phenotype}_stranded/{IndID}.{Rep}.plus.bw",
        bw_minus = "bigwigs/{Phenotype}_stranded/{IndID}.{Rep}.minus.bw"
    log:
        "logs/MakeBigwigs_stranded/{Phenotype}/{IndID}.{Rep}.log"

use rule MakeBigwigs_NormalizedToEdgeRFeatureCounts as MakeBigwigs_NormalizedToEdgeRFeatureCounts_stranded_nosplit with:
    params:
        GenomeCovArgs="",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
        BamToBigwigScript = "scripts/GenometracksByGenotype/BamToBigwig.sh"
    output:
        bw = "bigwigs_nosplit/{Phenotype}_stranded/{IndID}.{Rep}.plus.bw",
        bw_minus = "bigwigs_nosplit/{Phenotype}_stranded/{IndID}.{Rep}.minus.bw"
    log:
        "logs/MakeBigwigs_stranded_nosplit/{Phenotype}/{IndID}.{Rep}.log"

use rule MakeBigwigs_NormalizedToEdgeRFeatureCounts as MakeBigwigs_NormalizedToEdgeRFeatureCounts_stranded_nosplit_filtered with:
    params:
        GenomeCovArgs="",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
        BamToBigwigScript = "scripts/BamToBigwig_PrefilterBySize.sh"
    output:
        bw = "bigwigs_nosplit_filtered/{Phenotype}_stranded/{IndID}.{Rep}.plus.bw",
        bw_minus = "bigwigs_nosplit_filtered/{Phenotype}_stranded/{IndID}.{Rep}.minus.bw"
    log:
        "logs/MakeBigwigs_stranded_nosplit_filtered/{Phenotype}/{IndID}.{Rep}.log"

rule GatherUnsplitBigwigs:
    input:
        "bigwigs_nosplit/chRNA.Expression_stranded/NA18486.1.minus.bw",
        "bigwigs_nosplit_filtered/chRNA.Expression_stranded/NA18486.1.minus.bw"

rule GatherAllBigwigs:
    input:
        expand(
            "bigwigs/{Phenotype}/{IndID}.{Rep}.bw",
            zip,
            Phenotype=Fastq_samples["Phenotype"],
            IndID=Fastq_samples["IndID"],
            Rep=Fastq_samples["RepNumber"],
        ),
        expand(
            "bigwigs/{Phenotype}_stranded/{IndID}.{Rep}.plus.bw",
            zip,
            Phenotype=chRNASeqSamples_df["Phenotype"],
            IndID=chRNASeqSamples_df["IndID"],
            Rep=chRNASeqSamples_df["RepNumber"],
        ),
        expand(
            "bigwigs/{Phenotype}_stranded/{IndID}.{Rep}.minus.bw",
            zip,
            Phenotype=chRNASeqSamples_df["Phenotype"],
            IndID=chRNASeqSamples_df["IndID"],
            Rep=chRNASeqSamples_df["RepNumber"],
        ),
