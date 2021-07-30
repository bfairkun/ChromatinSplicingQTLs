rule ParseQTLtoolsOutputAndGetSE:
    """
    hyprcoloc needs beta and standard error for each snp/phenotype pair. QTLtools outputs beta, and a nominal P-value, from which we can calculate t-statistic and get SE. The covariates used in QTL mapping are needed to count the degrees of freedom to get appropriate t-distribution.
    """
    input:
        QTLtools_nominal_output = "QTLs/QTLTools/{Phenotype}/NominalPass_ForColoc.txt.gz",
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        "hyprcoloc/{Phenotype}/summarystats.txt.gz",
    log:
        "logs/ParseQTLtoolsOutputAndGetSE/{Phenotype}.log"
    shell:
        """
        python scripts/AddSEToQTLtoolsOutput.py {input.QTLtools_nominal_output} {input.cov} {input.vcf} {output} &> {log}
        """

#TODO: make a script to split/combine all summary statistics into a separate for each gene (each hyprcoloc attampt)
rule SplitAndCombineSummaryStatsPerGene:
    input:
        QTLtools_nominal_output = expand("QTLs/QTLTools/{Phenotype}/NominalPass_ForColoc.txt.gz", Phenotype=PhenotypesToColoc)
    output:
        directory("hyprcoloc/GenewiseSummaryStatsInput")
    log:
        "logs/SplitAndCombineSummaryStatsPerGene.log"
    shell:
        """
        python scripts/CombineAndSplitSummaryStatsForColoc.py {output}/ {input} &> {log}
        """
