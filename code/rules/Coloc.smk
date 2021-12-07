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
    This rule takes a long time to run. If we are going to be changing the
    pipeline a lot and rerunning this rule a lot, it might be convenient to
    parallelize this rule somehow computation somehow
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
    mamba create --name r_hyprcoloc -c r r-rmpfr r-iterpc r-tidyverse r-devtools r-pheatmap r-rcppeigen r-essentials
    For reasons I don't understand, conda won't export this environment to yaml
    with `conda export`. So I manually created an environment with the command
    above, then ran `conda list -e` and manually indented lines to conform to
    yaml to create the conda-compatible yaml file specified.
    """
    input:
        "envs/r_hyprcoloc.yml"
    output:
        touch("hyprcoloc/hyprcoloc_installed.touchfile")
    log:
        "logs/InstallHyprcoloc.log"
    conda:
        "../envs/r_hyprcoloc.yml"
    shell:
        """
        Rscript -e 'devtools::install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE, dependencies=F); install.packages("R.utils", repos = "http://cran.us.r-project.org")' &> {log}
        """

rule create_gwascoloc_bash_scripts:
    output:
        expand("hyprcoloc/Results/{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}/Chunks/{n}.sh", n=range(0,config["gwas_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc"
    run:
        from itertools import cycle
        gwas_bashscript_pairs = zip(gwas_df.index, cycle(output))
        for accession, out_f in gwas_bashscript_pairs:
            with open(out_f, 'a') as f:
                _ = f.write(f'Rscript scripts/hyprcoloc_gwas.R hyprcoloc/LociWiseSummaryStatsInput/ForGWASColoc/{accession}.txt.gz gwas_summary_stats/leadSnpWindowStats/{accession}.tsv.gz {out_f.rstrip(".sh")}.txt.gz\n')
        # If there are more output files than accession numbers, the extra
        # output files won't get made in the previous loop and snakemake will
        # complain of missing output files. as a fail safe, let's append to
        # each file in output, in effect making an empty file if a file wasn't
        # made in the for loop above
        for f in output:
            open(f, 'a').close()



rule gwas_coloc_chunk:
    input:
        bashscript = "hyprcoloc/Results/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/Chunks/{n}.sh",
        MolQTLSummaryStats = "hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}",
        Gwas_summary_stats =  expand("gwas_summary_stats/leadSnpWindowStats/{accession}.tsv.gz", accession=gwas_df.index),
        ModifiedCondaEnvConfirmation = "hyprcoloc/hyprcoloc_installed.touchfile"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc"
    output:
        "hyprcoloc/Results/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/Chunks/{n}.txt.gz"
    resources:
        mem_mb = 16000
    log:
        "logs/gwas_coloc_chunk/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/{n}.log"
    conda:
        "../envs/r_hyprcoloc.yml"
    shell:
        """
        bash {input.bashscript} &> {log}
        """

rule Gather_gwas_coloc_chunks:
    input:
        expand("hyprcoloc/Results/{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}/Chunks/{n}.txt.gz", n=range(0, config["gwas_coloc_chunks"]))
    output:
        "../output/hyprcoloc_results/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/hyprcoloc.results.txt.gz"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc"
    shell:
        """
        cat <(echo "GWASLeadSnpChrom_Pos_RefAllele_AltAllele_rsID_trait\tHyprcolocIteration\tColocalizedTraits\tPosteriorColocalizationPr\tRegionalAssociationPr\tTopCandidateSNP\tProportionPosteriorPrExplainedByTopSNP\tDroppedTrait") <(zcat {input}) | gzip - > {output}
        """
