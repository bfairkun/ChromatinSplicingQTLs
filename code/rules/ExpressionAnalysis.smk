rule Prepare_polyA_RNA_seq_PhenotypeTable:
    input:
        featureCounts = "featureCounts/polyA.Expression/Counts.txt",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        genes = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed",
    output:
        GeneList = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        PhenotypesBed = "QTLs/QTLTools/Expression.Splicing/AllReps.qqnorm.bed.gz",
        FirstReps = "QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/Prepare_polyA_RNA_seq_PhenotypeTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTablesFromFeatureCounts.R {input.featureCounts} rename_STAR_alignment_samples {input.gtf} {input.genes} {output.GeneList} {output.PhenotypesBed} {output.FirstReps}
        """

rule Prepare_chrRNA_RNA_seq_ExpressionPhenotypeTable:
    input:
        featureCounts = "featureCounts/chRNA.Expression/Counts.txt",
        GeneList = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
    output:
        FirstReps = "QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/Prepare_chrRNA_RNA_seq_ExpressionPhenotypeTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTableFromFeatureCounts_SubsetGeneList.R {input.featureCounts} {input.GeneList} {output.FirstReps}
        """

