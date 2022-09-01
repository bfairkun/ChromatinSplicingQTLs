rule featureCounts:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForPhenotype
    output:
        "featureCounts/{Phenotype}/Counts.txt"
    params:
        extraParams = GetFeatureCountsParams,
        paired = PairedEndParams
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(PhenotypeSet)
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.log"
    shell:
        """
        featureCounts {params.paired} {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

use rule featureCounts as featureCountsAtRegion with:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForRegion
    output:
        "featureCounts/{Region}/{Phenotype}/Counts.txt"
    log:
        "logs/featureCountsAtRegion/{Region}/{Phenotype}.log"
    wildcard_constraints:
        Phenotype = "chRNA.Expression|polyA.Expression|H3K27AC|H3K4ME3|H3K4ME1|ProCap"

rule Get2kbTSSRegions:
    """
    The GTFTools TSS bed file output from basic annotations GTF still has multiple TSS for some genes, presumably because sometimes more than one transript in the GTF is annotated as basic, and sometimes they have different TSS. To get just one TSS per gene (for simplicity of downstream analysis/interpretation) I will choose the TSS that is used in the most annotated transcripts, picking at random in case of ties
    """
    input:
        bed = "ReferenceGenome/Annotations/GTFTools_BasicAnnotations/gencode.v34.chromasomal.tss.bed",
        faidx = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        bed = "ReferenceGenome/Annotations/TSSRegions.bed",
        saf = "ReferenceGenome/Annotations/TSSRegions.saf"
    shell:
        """
        awk '{{print $1, $2,$3,$4, $6}}' {input.bed} | sort | uniq -c | shuf | sort -k6,6 -k1,1nr | sort -u -k6,6 | awk -v OFS='\\t' '{{print "chr"$2, $3, $4, $6, $5}}' | bedtools slop -b 999 -i - -g {input.faidx} | bedtools sort -i - > {output.bed}
        awk -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }} {{ print $4, $1,$2,$3,$6 }}' {output.bed} > {output.saf}
        """
