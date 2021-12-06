rule SplitAndCombineSummaryStatsPerLocus:
    """
    since hyprcoloc will require reading in summary stat data in R, it will be
    more convenient if the summary stats were organized into smaller files so
    we don't have to read in huge files into R. Also, it is sort of pointless
    to attempt colocalization on molecular traits that don't have an QTL
    signal. Therefore, this script reads the QTLtools output files for multiple
    phenotypes, and writes out the necessary summary stats for all phenotypes
    to new smaller files (one file for each gwas trait or gene) if they pass
    some minimum Pvalue treshold from the QTLtools permutation test.  TODO:
    This rule takes a long time to run. It might be convenient to parallelize
    computation somehow
    """
    input:
        zip( expand(
            "QTLs/QTLTools/{Phenotype}/{Pass}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}.txt.gz",
            Pass="PermutationPass",
            Phenotype=PhenotypesToColoc,
        ),
        expand(
            "QTLs/QTLTools/{Phenotype}/{Pass}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}.txt.gz",
            Pass="NominalPass",
            Phenotype=PhenotypesToColoc,
        ))
    output:
        directory("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}")
    log:
        "logs/SplitAndCombineSummaryStatsPerGene.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    params:
        MinNominalP = 0.01
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc|ForColoc"
    shell:
        """
        python scripts/CombineAndSplitSummaryStatsForGWASColoc.py {output}/ {params.MinNominalP} {input} &> {log}
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
