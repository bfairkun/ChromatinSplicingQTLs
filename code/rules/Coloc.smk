rule ParseQTLtoolsOutputAndGetSE:
    """
    hyprcoloc needs beta and standard error for each snp/phenotype pair. QTLtools outputs beta, and a nominal P-value, from which we can calculate t-statistic and get SE. The covariates used in QTL mapping are needed to count the degrees of freedom to get appropriate t-distribution.
    """
    input:
        QTLtools_nominal_output = "QTLs/QTLTools/{Phenotype}/NominalPass{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz",
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        "hyprcoloc/summarystats/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/{Phenotype}.txt.gz",
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor="|".join(["ForColoc", "ForGWASColoc"]),
    log:
        "logs/ParseQTLtoolsOutputAndGetSE/{FeatureCoordinatesRedefinedFor}/{Phenotype}/{QTLsGenotypeSet}.log"
    shell:
        """
        python scripts/AddSEToQTLtoolsOutput.py {input.QTLtools_nominal_output} {input.cov} {input.vcf} {output} &> {log}
        """

rule SplitAndCombineSummaryStatsPerGene:
    input:
        QTLtools_nominal_output = expand("hyprcoloc/summarystats/{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}/{Phenotype}.txt.gz", Phenotype=PhenotypesToColoc)
    output:
        directory("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}")
    log:
        "logs/SplitAndCombineSummaryStatsPerGene.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    shell:
        """
        python scripts/CombineAndSplitSummaryStatsForColoc.py {output}/ {input} &> {log}
        """

rule InstallHyprcoloc:
    """
    hyprcoloc r package is not on conda. This rule installs it on the conda environment. Here is a command to recreate the conda environment with dependencies (without hyprcoloc)
    mamba create --name r_hyprcoloc -c r r-rmpfr r-iterpc r-tidyverse r-devtools r-pheatmap r-rcppeigen
    """
    output:
        touch("hyprcoloc/hyprcoloc_installed.touchfile")
    log:
        "logs/InstallHyprcoloc.log"
    conda:
        "../envs/r_hyprcoloc.yml"
    shell:
        """
        Rscript -e 'devtools::install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE, dependencies=F)' &> {log}
        """
