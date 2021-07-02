## QTL checks
# snakemake GatherReportFiles --report report.html

# rule GatherExistingReportFiles:
#     input:
#         glob.glob("QTLs/QTLTools/*/FeatureSizes.png")

rule GatherReportFiles:
    input:
        expand("QTLs/QTLTools/{Phenotype}/FeatureSizes.pdf", Phenotype=["chRNA.Expression.Splicing", "polyA.IR", "H3K4ME3"])

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

