def GetMaxNumFeaturesForQTLCalling(wildcards):
    if wildcards.Phenotype == "CTCF":
        return 50000
    else:
        return 1000000

rule PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks:
    input:
        "featureCounts/{Phenotype}/Counts.txt"
    output:
        AllSamples = "QTLs/QTLTools/{Phenotype}/AllReps.qqnorm.bed.gz",
        FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        Phenotype = "|".join(ChromatinProfilingPhenotypes)
    params:
        max_features = GetMaxNumFeaturesForQTLCalling
    log:
        "logs/PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks/{Phenotype}.log"
    shell:
        """
        Rscript scripts/PreparePhenotypeTablesFromFeatureCounts_ChromatinProfilingPeaks.R {input} {params.max_features} {output.AllSamples} {output.FirstReps} &> {log}
        """


use rule PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks as MakeChromatinSpliceSitePhenotypes with:
    wildcard_constraints:
        Phenotype = 'H3K4ME1.5PrimeSS|H3K4ME1.3PrimeSS|H3K4ME3.5PrimeSS|H3K4ME3.3PrimeSS|H3K27AC.5PrimeSS|H3K27AC.3PrimeSS'

# rule GetGenotypePCs:
#     input:
#         FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
#         vcf = "Genotypes/1KG_GRCh38/7.vcf.gz"
