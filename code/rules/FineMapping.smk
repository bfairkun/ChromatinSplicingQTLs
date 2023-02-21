rule VCF2Txt:
    input:
        "Genotypes/1KG_GRCh38/{chrom}.vcf.gz"
    output:
        "FineMapping/Genotypes/1KG_GRCh38/{chrom}.txt.gz"
    log:
        "logs/FineMapping/vcf2txt.{chrom}.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        chrom = '|'.join([str(x) for x in range(1, 23)])
    shell:
        """
        python scripts/vcf2txt.py {input} {output} &> {log}
        """