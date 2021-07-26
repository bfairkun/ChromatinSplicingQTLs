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
        expand("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz", Phenotype=MyPhenotypes)
    output:
        report("QC/NumQTLsPerPhenotype.pdf", category="QTLMapping_QC")
    log:
        "logs/PlotNumQTLs.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Plot_QTLPermutationTestNumSig.R {output} {input} &> {log}
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
