rule Prepare_polyA_RNA_seq_PhenotypeTable_AndMakeGeneList:
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
        Rscript scripts/PreparePhenotypeTablesFromFeatureCounts.R {input.featureCounts} rename_STAR_alignment_samples {input.gtf} {input.genes} {output.GeneList} {output.PhenotypesBed} {output.FirstReps} &> {log}
        """

def GetFeatureCountsForRNASeqExpressionPhenotype(wildcards):
    if wildcards.Phenotype == "chRNA.Expression.Splicing":
        return "featureCounts/chRNA.Expression/Counts.txt"
    elif wildcards.Phenotype == "Expression.Splicing.Subset_YRI":
        return "featureCounts/Expression.Splicing/Counts.txt"
    else:
        return "featureCounts/{Phenotype}/Counts.txt"

def Get1KGMetadataFile(wildcards):
    if "YRI" in wildcards.Phenotype:
        return "../data/igsr_samples.tsv.gz"
    else:
        return []

rule Prepare_RNA_seq_ExpressionPhenotypeTable_ForGenesInList:
    input:
        featureCounts = GetFeatureCountsForRNASeqExpressionPhenotype,
        GeneList = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        YRI_List = Get1KGMetadataFile
    output:
        FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join([x for x in RNASeqPhenotypes_extended if x != 'chRNA.Expression.Splicing'])
    log:
        "logs/Prepare_chrRNA_RNA_seq_ExpressionPhenotypeTable/{Phenotype}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTableFromFeatureCounts_SubsetGeneList.R {input.featureCounts} {input.GeneList} {output.FirstReps} {input.YRI_List} &> {log}
        """




rule Prepare_chRNA_ExpressionPhenotypes:
    input:
        "featureCounts/chRNA.Expression/Counts.txt",
        "featureCounts/chRNA.Expression_ncRNA/Counts.txt",
        "featureCounts/chRNA.Expression_annotated_ncRNA/Counts.txt",
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
        "NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "RPKM_tables/chRNA.RPKM.bed.gz",
    log:
        "logs/Prepare_chRNA_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem = 48000,
    shell:
        """
        mkdir -p RPKM_tables/;
        Rscript scripts/Prepare_chRNA_Phenotypes.R &> {log}
        """
        
rule tabix_genes_bed:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
    output:
        bed = "ExpressionAnalysis/polyA/ExpressedGeneList.bed.gz",
        tbi = "ExpressionAnalysis/polyA/ExpressedGeneList.bed.gz.tbi"
    shell:
        """
        cat {input} | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

