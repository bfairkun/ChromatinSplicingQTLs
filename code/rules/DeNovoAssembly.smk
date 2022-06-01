rule Stringtie:
    input:
        bam = GetBamForPhenotype,
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    output:
        "ReferenceGenome/Annotations/Assembly/{Phenotype}.stringtie.assembly.gtf"
    log:
        "logs/stringtie.{Phenotype}.log"
    resources:
        mem_mb = 58000
    threads:
        8
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing"
    shell:
        """
        {config[stringtie]} --rf -c 10 -f 0.99 -t -m 2000 -v -o {output} -p {threads} -G {input.gtf} Alignments/STAR_Align/chRNA.Expression.Splicing/NA19153/1/Filtered.bam &> {log}
        """