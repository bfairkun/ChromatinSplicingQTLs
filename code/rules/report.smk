## QTL checks
# snakemake GatherReportFiles --report report.html

# rule GatherExistingReportFiles:
#     input:
#         glob.glob("QTLs/QTLTools/*/FeatureSizes.png")

rule GatherReportFiles:
    input:
        expand("QTLs/QTLTools/{Phenotype}/FeatureSizes.pdf", Phenotype=MyPhenotypes),
        "QC/NumQTLsPerPhenotype.pdf",
        expand("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.Pvals_{Plot}.pdf", Phenotype=MyPhenotypes, Plot=["hist", "qq"])

rule PlotFeatureSize:
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
    output:
        report("QTLs/QTLTools/{Phenotype}/FeatureSizes.pdf", category="QTLMapping_QC")
    log:
        "logs/PlotFeatureSize/{Phenotype}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PlotFeatureSizes_FromPhenotypeTable.R {input} {output} &> {log}
        """


rule PlotNumQTLs:
    input:
        # expand("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz", Phenotype=MyPhenotypes)
        expand("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz", Phenotype=MyPhenotypes)
    output:
        pdf = report("QC/NumQTLsPerPhenotype.pdf", category="QTLMapping_QC"),
        bed = "QC/AllQTLPhenotypes.PermutationTest.bed.gz"
    log:
        "logs/PlotNumQTLs.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Plot_QTLPermutationTestNumSig.R {output.pdf}
        {output.bed} {input} &> {log}
        """

rule PlotQTLsPermutationPvalsHist:
    input:
       table = "QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz",
    output:
        hist = report("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.Pvals_hist.pdf", category="QTLMapping_QC"),
        qq = report("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.Pvals_qq.pdf", category="QTLMapping_QC")
    log:
        "logs/PlotQTLsPermutationPvalsHist/{Phenotype}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Plot_QTLPermutationTestPvals.R {input} {output.qq} {output.hist} &> {log}
        """

Rmds, = glob_wildcards("../analysis/{fn}.Rmd")
rule BuildRmd:
    """
    Build site for easy sharing of analysis done with Rmarkdown. If possible,
    make sure the Rmd only references relatively small files that are tracked
    with git (eg in `../output/` or `../data/`) and use relative filepaths in
    the Rmd. This way  one could succesfully run the Rmd files or use this rule
    this rule to build Rmd files after pulling/cloning the repo. I think saving
    could be useful to do exploratory data analysis between Carlos and I
    without necessarily running the computationally intensive parts of the
    Snakemake pipeline on each of our clones.  without running the rest of the
    snakemake. I have been writing my Rmd files assuming `../analysis` is the
    working directory.  Note that input files that might be referenced within an
    Rmd are not specified in the snakemake. This rule might have to be manually
    rerun if important input files referenced in the Rmd get updated.
    """
    input:
        "../analysis/{fn}.Rmd",
    output:
        "../docs/{fn}.html"
    log:
        "logs/BuildRmd/{fn}.log"
    wildcard_constraints:
        fn = "|".join(Rmds)
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript -e 'setwd("../analysis"); workflowr::wflow_build("{wildcards.fn}.Rmd")' &> {log}
        """

rule CollectBuiltRmds:
    input:
        expand("../docs/{fn}.html", fn=Rmds)
