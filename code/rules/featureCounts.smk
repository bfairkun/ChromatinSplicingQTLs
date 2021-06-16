rule featureCounts:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForPhenotype
    output:
        "featureCounts/{Phenotype}/Counts.txt"
    params:
        extraParams = GetFeatureCountsParams,
    threads:
        8
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.log"
    shell:
        """
        featureCounts -p {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """
