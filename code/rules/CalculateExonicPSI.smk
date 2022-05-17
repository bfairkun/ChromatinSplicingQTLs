rule GetInternalExons:
    input:
        exons_in = "ReferenceGenome/Annotations/GTFTools_BasicAnnotations/gencode.v34.chromasomal.exons.bed"
    output:
        "SplicingAnalysis/Annotations/InternalExons.bed"
    log:
        "logs/GetInternalExons.log"
    conda:
        "../envs/r_slopes.yml"
    shell:
        """
        Rscript scripts/GetInternalExons.R {input} {output} &> {log}
        """

rule RemoveOverlappingInternalExons:
    input:
        "SplicingAnalysis/Annotations/InternalExons.bed"
    output:
        "SplicingAnalysis/Annotations/InternalExons.noOverlapping.saf"
    shell:
        """
        bedtools intersect -a {input} -b {input} -s -c -sorted | awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }} $7==1 {{ print $4, $1, $2, $3, $6 }}'  > {output}
        """

rule CountSkippedReadsForInternalExons:
    """
    slop expand the exon bed by 20bp on each side and sum the number of splice
    junction intron reads that completely overlap the exon (skipped reads)
    """
    input:
        saf = "SplicingAnalysis/Annotations/InternalExons.noOverlapping.saf",
        JuncBed = "SplicingAnalysis/leafcutter/regtools_annotate/comprehensive/{Phenotype}_{IndID}_{Rep}.bed.gz"
    output:
        "SplicingAnalysis/Annotations/InternalExons.noOverlapping.SkippedReads/{Phenotype}_{IndID}_{Rep}.bed.gz"
    shell:
        """
        bedtools map -null 0 -o sum -f 1 -s -a <(awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $2, $3-20, $4+20, $1, ".", $5}}' {input.saf})  -b <(zcat {input.JuncBed} | awk 'NR>1' | bedtools sort -i -) | gzip - > {output}
        """

rule GatherRegtoolsAnnotatechRNA:
    input:
        expand("SplicingAnalysis/Annotations/InternalExons.noOverlapping.SkippedReads/{Phenotype}_{IndID}_{Rep}.bed.gz", zip,Phenotype=chRNASeqSamples_df["Phenotype"],IndID=chRNASeqSamples_df["IndID"],Rep=chRNASeqSamples_df["RepNumber"]),
        "featureCounts/AtInternalExons/chRNA.Expression/Counts.txt",
        "featureCounts/AtInternalExons/polyA.Expression/Counts.txt"


def GetRepresentativeBwForAssay(wildcards):
    if wildcards.Assay == "H3K36ME3":
        return "bigwigs/H3K36ME3/NA18486.1.bw"
    elif wildcards.Assay == "HEK36ME3_ENCODE":
        return "/project2/yangili1/bjf79/20201123_CheckH3K36me3_CutAndTag/code/ENCODE/H3K36me3.rep1.bw"
    elif wildcards.Assay == "H3K4ME3":
        return "bigwigs/H3K4ME3/NA18489.1.bw"

rule ExonMetaPlot:
    input:
        beds = expand("SplicingAnalysis/Annotations/InternalExons.PSI.{RangeList}.bed", RangeList=["0-20", "20-40", "40-60", "60-80", "80-100"] ),
        bw = GetRepresentativeBwForAssay
    output:
        MetaProfile = "deepTools/Metplots/CassetteExonic/{Assay}.png",
        Metamatrix = "deepTools/Metplots/CassetteExonic/{Assay}.gz"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/deepTools/ExonMetaPlot/{Assay}.log"
    shell:
        """
        computeMatrix reference-point -R {input.beds} -S {input.bw} -o {output.Metamatrix} -a 1000 -b 1000 --missingDataAsZero --referencePoint center
        plotHeatmap -m {output.Metamatrix} -o {output.MetaProfile} -z 0-20 20-40 40-60 60-80 80-100 --heatmapHeight 10 --refPointLabel CenterOfExon
        """

rule MetagenePlot:
    input:
        beds = expand("ExpressionAnalysis/chRNAExpression_Quartile{ntile}.bed", ntile=[1,2,3,4] ),
        bw = GetRepresentativeBwForAssay
    output:
        MetaProfile = "deepTools/Metplots/Metagene/{Assay}.png",
        Metamatrix = "deepTools/Metplots/Metagene/{Assay}.gz"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/deepTools/Metagene/{Assay}.log"
    shell:
        """
        computeMatrix scale-regions -R {input.beds} -S {input.bw} -o {output.Metamatrix} --regionBodyLength 5000 -a 4000 -b 2000 --missingDataAsZero
        plotHeatmap -m {output.Metamatrix} -o {output.MetaProfile} -z Q1 Q2 Q3 Q4 --heatmapHeight 10
        """

rule CollectMetaExonPlots:
    input:
        expand("deepTools/Metplots/CassetteExonic/{Assay}.gz", Assay=["H3K36ME3", "HEK36ME3_ENCODE", "H3K4ME3"]),
        expand("deepTools/Metplots/Metagene/{Assay}.gz", Assay=["H3K36ME3", "HEK36ME3_ENCODE", "H3K4ME3"])
