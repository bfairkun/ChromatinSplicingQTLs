def GetPSIBed(wildcards):
    PhenotypeToFileWildcard = dict(zip(["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"],["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"]))
    return f"SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{PhenotypeToFileWildcard[wildcards.Phenotype]}.bed.gz"

rule LeafcutterPSI_table:
    """
    PSI Junction excision ratio. beta will measure delta PSI
    """
    input:
        PSI = GetPSIBed,
        standardized = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"]),
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        #Script to impute NA with median
        Rscript scripts/PrepareUnstandardizedPSIPhenotypeTables.R {input.PSI} {input.standardized} {output}
        """

rule Subset_YRI_phenotype_table_unstandardized:
    input:
        input_file = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz",
        igsr = '../data/igsr_samples.tsv.gz'
    output:
        "QTLs/QTLTools/{Phenotype}.Subset_YRI/OnlyFirstRepsUnstandardized.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "polyA.IR", "polyA.IER", "polyA.Splicing.5PrimeSS", "polyA.Splicing.3PrimeSS", "Expression.Splicing"])
    log:
        "logs/Subsample_YRI_unstandardized/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/Subset_YRI.R {input.input_file} {output}  &> {log}
        """


rule Unstandardized_CPM_Tables:
    """
    LogCPM. beta will measure fold change (additive change in log scale)
    """
    input:
        Counts = "featureCounts/{Phenotype}/Counts.txt",
        standardized = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join([i for i in ChromatinProfilingPhenotypes if i != "H3K36ME3"])
    log:
        "logs/Unstandardized_RPKM_Tables/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/PrepareLogCPM_PhenotypeTables.R {input.Counts} {input.standardized} {output} &> {log}
        """

rule Unstandardized_H3K36ME3_CPM_Table:
    input:
        Counts = "MiscCountTables/H3K36ME3.bed",
        standardized = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "H3K36ME3"
    log:
        "logs/Unstandardized_H3K36ME3_Tables/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/PrepareLogRPKM_H3K36ME3_PhenotypeTables.R {input.Counts} {input.standardized} {output} &> {log}
        """

def GetCountTable(wildcards):
    if wildcards.Phenotype == "polyA.Expression.Splicing":
        return "featureCounts/polyA.Expression/Counts.txt"
    elif wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "featureCounts/chRNA.Expression/Counts.txt"
    else:
        return "featureCounts/{Phenotype}/Counts.txt"

rule Unstandardized_RPKM_Tables:
    """
    LogRPKM. beta will measure fold change (additive change in log scale)
    """
    input:
        Counts = GetCountTable,
        standardized = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join([i for i in RNASeqPhenotypes if i != "ProCap"])
    log:
        "logs/Unstandardized_RPKM_Tables/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/PrepareLogRPKM_PhenotypeTables.R {input.Counts} {input.standardized} {output} &> {log}
        """

rule CollectUnstandardizedTables:
    input:
        expand("QTLs/QTLTools/{Phenotype}/OnlyFirstRepsUnstandardized.sorted.qqnorm.bed.gz", Phenotype=[i for i in RNASeqPhenotypes if i != "ProCap"] + ChromatinProfilingPhenotypes + ["Expression.Splicing.Subset_YRI", "polyA.Splicing.Subset_YRI"] + ["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"] ),
        # expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}_{StandardizedOrUnstandardized}.txt.gz",
        #         Phenotype=[i for i in RNASeqPhenotypes if i != "ProCap"] + ChromatinProfilingPhenotypes + ["Expression.Splicing.Subset_YRI", "polyA.Splicing.Subset_YRI"] + ["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"],
        #         Pass = "NominalPass",
        #         QTLsGenotypeSet="",
        #         FeatureCoordinatesRedefinedFor="",
        #         StandardizedOrUnstandardized="Unstandardized"
        #         ),
        expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}_{StandardizedOrUnstandardized}.txt.gz",
                Phenotype=ChromatinProfilingPhenotypes + ["polyA.Splicing.Subset_YRI"] + ["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"],
                Pass = "NominalPass",
                QTLsGenotypeSet="",
                FeatureCoordinatesRedefinedFor="",
                StandardizedOrUnstandardized="Unstandardized"
                )

            # ChromatinProfilingPhenotypes +
            # [i for i in RNASeqPhenotypes if i != "ProCap"] +
            # ["polyA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.Splicing"] +
            # ["polyA.Splicing.Subset_YRI", "Expression.Splicing.Subset_YRI"]
            # )


