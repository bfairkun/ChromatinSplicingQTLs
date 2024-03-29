rule SplitAndCombineSummaryStatsPerLocus:
    """
    since hyprcoloc will require reading in summary stat data in R, it will be
    more convenient if the summary stats were organized into smaller files so
    we don't have to read in huge files into R. Also, it is sort of pointless
    to attempt colocalization on molecular traits that don't have an QTL
    signal. Therefore, this script reads the QTLtools output files for multiple
    phenotypes, and writes out the necessary summary stats for all phenotypes
    to new smaller files (one file for each gwas trait or gene) if they pass
    some minimum Pvalue treshold from the QTLtools permutation test.
    """
    input:
        zip( expand(
            "QTLs/QTLTools/{Phenotype}/{Pass}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}.txt.gz",
            Pass="PermutationPass",
            Phenotype=PhenotypesToColoc,
        ),
        expand(
            "QTLs/QTLTools/{Phenotype}/{Pass}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}Chunks/{{n}}.txt",
            Pass="NominalPass",
            Phenotype=PhenotypesToColoc,
        )) 
    output:
        directory("hyprcoloc/LociWiseSummaryStatsInput/Chunks/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/{n}")
    log:
        "logs/SplitAndCombineSummaryStatsPerGene.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.{n}.log"
    params:
        MinNominalP = 0.01
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc|ForColoc"
    shell:
        """
        python scripts/CombineAndSplitSummaryStatsForGWASColoc.py {output}/ {params.MinNominalP} {input} &> {log}
        """

rule GatherSummaryStatsFileChunks:
    input:
        expand("hyprcoloc/LociWiseSummaryStatsInput/Chunks/{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}/{n}", n=range(1, 1+N_PermutationChunks) )
    output:
        directory("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}")
    log:
        "logs/GatherSummaryStatsFileChunks/{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc|ForColoc"
    shell:
        """
        python scripts/MergeSummaryStatChunks.py {output} {input} &> {log}
        """

CalculateLDChunks = 100
rule CalculateLociwiseLDMat:
    input:
        MolQTLSummaryStats = "hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}",
    output:
        MolQTLSummaryStats = touch("hyprcoloc/LociWiseSummaryStatsLD/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/.touch_{n}"),
    log:
        "logs/CalculateLociwiseLDMat/{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}/{n}"
    resources:
        mem_mb = 16000
    shell:
        """
        python scripts/CalculateLD_PerGeneWindow_Chunks.py {input.MolQTLSummaryStats}/ hyprcoloc/LociWiseSummaryStatsLD/{wildcards.QTLsGenotypeSet}{wildcards.FeatureCoordinatesRedefinedFor}/ {wildcards.n} {CalculateLDChunks} &> {log}
        """

rule GatherLDMatrices:
    input:
        expand("hyprcoloc/LociWiseSummaryStatsLD/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}/.touch_{n}", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc" ,n=range(0,CalculateLDChunks))

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


rule create_genewise_summarystats_listchunks:
    input:
        MolQTLSummaryStats = GetColocTsvFormattedString("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}"),
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    output:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.list.txt", n=range(0,config["genewise_coloc_chunks"]))
    run:
        from itertools import cycle
        import glob
        genefile_outfile_pairs = zip(glob.glob(input.MolQTLSummaryStats + "/*.txt.gz"), cycle(output))
        for genefile, out_f in genefile_outfile_pairs:
            with open(out_f, 'a') as f:
                _ = f.write(f'{genefile}\n')
        # If there are more output files than accession numbers, the extra
        # output files won't get made in the previous loop and snakemake will
        # complain of missing output files. as a fail safe, let's append to
        # each file in output, in effect making an empty file if a file wasn't
        # made in the for loop above
        for f in output:
            open(f, 'a').close()

rule create_gwascoloc_bash_scripts:
    output:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.sh", n=range(0,config["gwas_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc",
        ColocName = '|'.join(colocs_gwas.index)
    params:
        PhenotypesToColoc = GetMolPhenotypesToColoc
    run:
        from itertools import cycle
        gwas_bashscript_pairs = zip(gwas_df.index, cycle(output))
        for accession, out_f in gwas_bashscript_pairs:
            with open(out_f, 'a') as f:
                _ = f.write(f'Rscript scripts/hyprcoloc_gwas2.R hyprcoloc/LociWiseSummaryStatsInput/ForGWASColoc/{accession}.txt.gz gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz {out_f.rstrip(".sh")}.txt.gz "{params.PhenotypesToColoc}"\n')
        # If there are more output files than accession numbers, the extra
        # output files won't get made in the previous loop and snakemake will
        # complain of missing output files. as a fail safe, let's append to
        # each file in output, in effect making an empty file if a file wasn't
        # made in the for loop above
        for f in output:
            open(f, 'a').close()


rule gwas_coloc_chunk:
    input:
        bashscript = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.sh",
        MolQTLSummaryStats = GetColocTsvFormattedString("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}"),
        Gwas_summary_stats =  expand("gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz", accession=gwas_df.index),
        ModifiedCondaEnvConfirmation = "hyprcoloc/hyprcoloc_installed.touchfile"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc",
        ColocName = '|'.join(colocs_gwas.index)
    output:
        "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.txt.gz"
    resources:
        mem_mb = lambda wildcards, attempt: 32000 if int(attempt) == 1 else 48000
    log:
        "logs/gwas_coloc_chunk/{FeatureCoordinatesRedefinedFor}/{ColocName}/{n}.log"
    conda:
        "../envs/r_hyprcoloc.yml"
    shell:
        """
        bash {input.bashscript} &> {log}
        """

rule genewise_coloc_chunk:
    input:
        MolQTLSummaryStats = GetColocTsvFormattedString("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}"),
        summarystatslist = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.list.txt",
        ModifiedCondaEnvConfirmation = "hyprcoloc/hyprcoloc_installed.touchfile"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    output:
        clusters = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.txt.gz",
        snpscores = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.snpscores.txt.gz",
        pairwisecor = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/Chunks/{n}.pairwisecor.txt.gz"
    resources:
        mem_mb = lambda wildcards, attempt: 32000 if int(attempt) == 1 else 42000
    params:
        PhenotypesToColoc = GetMolPhenotypesToColoc,
        Threshold = lambda wildcards: colocs_df.loc[wildcards.ColocName]['Threshold']
    log:
        "logs/gwas_coloc_chunk/{FeatureCoordinatesRedefinedFor}/{ColocName}/{n}.log"
    conda:
        "../envs/r_hyprcoloc.yml"
    shell:
        """
        Rscript scripts/hyprcoloc_genewise2.R {input.summarystatslist} {output.clusters} {output.snpscores} {output.pairwisecor} {params.Threshold} '{params.PhenotypesToColoc}' &> {log}
        """

rule Gather_gwas_coloc_chunks:
    input:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.txt.gz", n=range(0, config["gwas_coloc_chunks"]))
    output:
        "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/results.txt.gz"
    params:
        header = "GWASLeadSnpChrom_Pos_RefAllele_AltAllele_rsID_trait\tHyprcolocIteration\tColocalizedTraits\tPosteriorColocalizationPr\tRegionalAssociationPr\tTopCandidateSNP\tProportionPosteriorPrExplainedByTopSNP\tDroppedTrait"
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc",
        ColocName = '|'.join(colocs_gwas.index)
    shell:
        """
        cat <(echo "{params.header}") <(zcat {input}) | gzip - > {output}
        """

use rule Gather_gwas_coloc_chunks as Gather_genewise_coloc_chunks with:
    input:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.txt.gz", n=range(0, config["genewise_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    params:
        header = "GeneLocus\tHyprcolocIteration\tPosteriorColocalizationPr\tRegionalAssociationPr\tTopCandidateSNP\tProportionPosteriorPrExplainedByTopSNP\tTrait"

use rule Gather_gwas_coloc_chunks as Gather_genewise_coloc_chunks_snpscores with:
    input:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.snpscores.txt.gz", n=range(0, config["genewise_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    output:
        "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/snpscores.txt.gz"
    params:
        header = "snp\tColocalizedCluster\tFinemapPr\tLocus"

use rule Gather_gwas_coloc_chunks as Gather_genewise_coloc_chunks_pairwisecor with:
    input:
        expand("hyprcoloc/Results/{{FeatureCoordinatesRedefinedFor}}/{{ColocName}}/Chunks/{n}.pairwisecor.txt.gz", n=range(0, config["genewise_coloc_chunks"]))
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    output:
        "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/pairwisecor.txt.gz"
    params:
        header = "Trait1\tTrait2\tcor.z.pearson\tcor.z.spearman\tcor.logp.pearson\tGeneLocus"

rule TidyColocalizedTraitsAndSummaryStats:
    input:
        ColocResults = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/results.txt.gz",
        MolQTLSummaryStats = GetColocTsvFormattedString("hyprcoloc/LociWiseSummaryStatsInput/{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}"),
    log:
        "logs/TidyColocalizedTraitsAndSummaryStats/{FeatureCoordinatesRedefinedFor}/{ColocName}.log"
    output:
        "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/tidy_results_OnlyColocalized.txt.gz"
    resources:
        mem_mb = 16000
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/TidyGenewiseColocs.R {input.ColocResults} {input.MolQTLSummaryStats}/ {output} &> {log}
        """

rule CopyGenewiseHyprcolocResultsToTrackedOutputDir:
    input:
        results = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/results.txt.gz",
        snpscores = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/snpscores.txt.gz",
        tidy_results = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/tidy_results_OnlyColocalized.txt.gz",
    output:
        results = "../output/hyprcoloc_results/{FeatureCoordinatesRedefinedFor}/{ColocName}/hyprcoloc.results.txt.gz",
        snpscores = "../output/hyprcoloc_results/{FeatureCoordinatesRedefinedFor}/{ColocName}/snpscores.txt.gz",
        tidy_results = "../output/hyprcoloc_results/{FeatureCoordinatesRedefinedFor}/{ColocName}/hyprcoloc.results.OnlyColocalized.Stats.txt.gz",
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForColoc",
        ColocName = '|'.join(colocs_genewise.index)
    shell:
        """
        cp {input.results} {output.results}
        cp {input.snpscores} {output.snpscores}
        cp {input.tidy_results} {output.tidy_results}
        """

rule CopyGwasHyprcolocResultsToTrackedOutputDir:
    input:
        results = "hyprcoloc/Results/{FeatureCoordinatesRedefinedFor}/{ColocName}/results.txt.gz",
    output:
        results = "../output/hyprcoloc_results/{FeatureCoordinatesRedefinedFor}/{ColocName}/hyprcoloc.results.txt.gz",
    wildcard_constraints:
        FeatureCoordinatesRedefinedFor = "ForGWASColoc",
        ColocName = '|'.join(colocs_gwas.index)
    shell:
        """
        cp {input.results} {output.results}
        """

rule GatherAllHyprcolocResultsToTrackedOutputDir:
    input:
        expand("../output/hyprcoloc_results/{FeatureCoordinatesRedefinedFor}/{ColocName}/hyprcoloc.results.txt.gz", FeatureCoordinatesRedefinedFor = "ForColoc", ColocName = colocs_genewise.index)
