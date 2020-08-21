rule MergeGEUVADIS_ForSTAR_WASP:
    input:
        expand("Genotypes/GEUVADIS_1000G/{Chrom}.vcf.gz", Chrom=["chr" + i for i in autosomes]),
    output:
        vcf = "Genotypes/GEUVADIS_1000G/All.vcf.gz",
        tbi = "Genotypes/GEUVADIS_1000G/All.vcf.gz.tbi",
    log:
        "logs/Genotypes/MergeGEUVADIS_ForSTAR_WASP.log"
    shell:
        """
        bcftools concat -a -O z {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """
