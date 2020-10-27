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

rule MakeChrTranslationFile:
    output:
        chrNames = "ReferenceGenome/Other/AutosomesTranslation_UCSC_to_Ensemble.txt"
    run:
        with open(output[0], "w") as f_out:
            for i in list(range(1,23)) + ["X", "Y"]:
                _ = f_out.write("{}\tchr{}\n".format(i,i))

rule MakeVcfForWASP:
    """
    as per this discussion (https://github.com/alexdobin/STAR/issues/772), STAR
    requires a vcf with the 10th field marked as HET (0|1) for every variant
    """
    input:
        vcf = "Genotypes/GEUVADIS_1000G/All.vcf.gz",
        tbi = "Genotypes/GEUVADIS_1000G/All.vcf.gz.tbi",
        chrNames = "ReferenceGenome/Other/AutosomesTranslation_UCSC_to_Ensemble.txt"
    output:
        DummyVcf = "Genotypes/GEUVADIS_1000G/All.Dummy.vcf"
    shell:
        """
        cat <(bcftools view -h -G {input.vcf} | bcftools annotate --rename-chrs {input.chrNames} | awk -F'\\t' -v OFS='\\t' '$0!~/^#CHROM/ {{print}} $0~/^#CHROM/ {{print $0, "FORMAT"}}' ) <(bcftools annotate --rename-chrs {input.chrNames} {input.vcf} | bcftools view -H -G - | sed 's/$/\\t1\|0/' ) > {output.DummyVcf}
        """


