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
        "featureCounts/AtInternalExons/chRNA.Expression/Counts.txt"


