rule GatherNormFactors:
    input:
        expand("featureCounts/{Phenotype}/NormFactors.tsv", Phenotype=["polyA.Expression", "chRNA.Expression"])

def GetGeneList(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes_extended:
        return "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    else:
        return []

def GetGeneListCol(wildcards):
    if wildcards.Phenotype in RNASeqPhenotypes:
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
        "../envs/r_essentials.yml"
    priority: 1
    shell:
        """
        Rscript scripts/CalculateNormFactorsForBigwig.R {input.Counts} {output} {input.OptionalGeneList} {params.OptionalGeneListCol} &> {log}
        """

# rule MakeBigwigs:
#     """
#     Scale bigwig to base coverage per billion chromosomal reads
#     """
#     input:
#         fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
#         bam = GetBamForBigwig,
#         bai = GetBaiForBigwig
#     params:
#         GetBigwigParams
#     wildcard_constraints:
#         Phenotype = "|".join(ChromatinProfilingPhenotypes)
#     output:
#         bw = "bigwigs/{Phenotype}/{IndID}.{Rep}.bw"
#     log:
#         "logs/MakeBigwigs/{Phenotype}/{IndID}.{Rep}.log"
#     resources:
#         mem = 48000
#     shell:
#         """
#         scripts/GenometracksByGenotype/BamToBigwig.sh {input.fai} {input.bam} {output.bw} {params} -scale $(bc <<< "scale=3;1000000000/$(samtools idxstats {input.bam} | awk '$1 ~ "^chr" {{sum+=$2}} END{{printf sum}}')") &> {log}
#         """

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
        GetBigwigParams
    # wildcard_constraints:
    #     Phenotype = "|".join(RNASeqPhenotypes)
    output:
        bw = "bigwigs/{Phenotype}/{IndID}.{Rep}.bw"
    log:
        "logs/MakeBigwigs/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 52000
    shell:
        """
        scripts/GenometracksByGenotype/BamToBigwig.sh {input.fai} {input.bam} {output.bw} {params} -scale $(bc <<< "scale=3;1000000000/$(grep '{input.bam}' {input.NormFactorsFile} | awk '{{print $2}}')") &> {log}
        """
