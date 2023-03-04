rule MakeUniversalVcfForEUR:
    """
    Make EUR vcf.
    """
    input:
        vcf = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz.tbi",
        metadata = "../data/igsr_samples.tsv.gz"
    output:
        vcf = "Genotypes/1KG_GRCh38_SubsetEUR/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38_SubsetEUR/{chrom}.vcf.gz.tbi"
    resources:
        mem_mb = 8000
    log:
        "logs/MakeUniversalVcfForTesting/{chrom}.log"
    params:
        filters = "-q 0.01:minor"
    shell:
        """
        (bcftools view --force-samples -S <(zcat {input.metadata} | awk -F'\\t' '$6=="EUR" {{ print $1 }}') -O z {params} {input.vcf} > {output.vcf} ) &> {log}
        tabix -p vcf {output.vcf}
        """
        
        
rule GatherWholeGenomeVcfEUR:
    input:
        vcf = expand("Genotypes/1KG_GRCh38_SubsetEUR/{chrom}.vcf.gz", chrom=autosomes),
        tbi = expand("Genotypes/1KG_GRCh38_SubsetEUR/{chrom}.vcf.gz.tbi", chrom=autosomes)
    output:
        vcf = "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz.tbi"
    log:
        "logs/GatherWholeGenomeVcfByPhenotype/EUR.log"
    resources:
        mem_mb = 16000
    shell:
        """
        bcftools concat -O z {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf}
        """

rule EUR_VCF2Txt:
    input:
        "QTLs/QTLTools/Expression.Splicing.Subset_EUR/Genotypes/WholeGenome.vcf.gz"
    output:
        gz = temp("FineMapping/Genotypes/1KG_GRCh38/EUR.txt.gz"),
        bgz = "FineMapping/Genotypes/1KG_GRCh38/EUR.txt.bgz"
    log:
        "logs/FineMapping/vcf2txt.EUR.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/vcf2txt.py {input} {output.gz} &> {log};
        (gunzip -c {output.gz} | bgzip > {output.bgz}) &>> {log}
        """

rule TabixSNPTables:
    input:
        "FineMapping/Genotypes/1KG_GRCh38/EUR.txt.bgz"
    output:
        "FineMapping/Genotypes/1KG_GRCh38/EUR.txt.bgz.tbi"
    log:
        "logs/FineMapping/vcftabix.EUR.log"
    resources:
        mem_mb = 12000
    shell:
        """
        tabix -s 1 -b 2 -e 2 {input} &> {log}
        """

rule run_susie:
    input:
        nominal = "QTLs/QTLTools/Expression.Splicing.Subset_EUR/NominalPass.txt.gz",
        permutation = "QTLs/QTLTools/Expression.Splicing.Subset_EUR/PermutationPass.FDR_Added.txt.gz", 
        genotype = "FineMapping/Genotypes/1KG_GRCh38/EUR.txt.bgz"
    output:
        "FineMapping/susie_runs/susie_output.tab"
    log:
        "logs/FineMapping/run_susie.log"
    resources:
        mem_mb = 32000
    shell:
        """
        Rscript scripts/run_susie.R {input.nominal} {input.permutation} {input.genotype} {output} &> {log}
        """
        
        
        
        
        
        
        
        
        
