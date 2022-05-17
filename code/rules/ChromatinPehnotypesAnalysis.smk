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
        Phenotype = "|".join([p for p in ChromatinProfilingPhenotypes if p != 'H3K36ME3'])
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
        Phenotype = 'H3K4ME1.5PrimeSS|H3K4ME1.3PrimeSS|H3K4ME3.5PrimeSS|H3K4ME3.3PrimeSS|H3K27AC.5PrimeSS|H3K27AC.3PrimeSS|H3K36ME3.5PrimeSS|H3K36ME3.3PrimeSS'

# rule GetGenotypePCs:
#     input:
#         FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
#         vcf = "Genotypes/1KG_GRCh38/7.vcf.gz"

H3K36ME3_IndIDs = Fastq_samples.loc[(Fastq_samples['Phenotype']=='H3K36ME3') & (Fastq_samples['RepNumber']==1)]['IndID'].tolist()
rule PrepareH3K36ME3PhenotypeTable:
    input:
        bams = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam", IndID = H3K36ME3_IndIDs),
        bais = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam.bai", IndID = H3K36ME3_IndIDs),
        bed = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        "MiscCountTables/H3K36ME3.bed"
    log:
        "logs/PrepareH3K36ME3PhenotypeTable.log"
    params:
        header = '\t'.join(["#Chr", "start", "end", "pid", "gid", "strand"] + H3K36ME3_IndIDs) + '\n'
    shell:
        """
        printf '{params.header}' > {output}
        (bedtools multicov -bed {input.bed} -bams {input.bams} >> {output}) 2> {log}
        """

rule PreparePhenotypeTableForQTLToolsFromWellFormatedCountTable:
    input:
        "MiscCountTables/H3K36ME3.bed"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        Phenotype = "H3K36ME3"
    log:
        "logs/PreparePhenotypeTableForQTLToolsFromWellFormatedCountTable/{Phenotype}.log"
    shell:
        """
        Rscript scripts/PreparePhenotypeTableFromWellFormattedTable.R {input} {output} &> {log}
        """

