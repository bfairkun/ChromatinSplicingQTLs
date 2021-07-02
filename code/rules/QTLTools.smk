rule SortQTLtoolsPhenotypeTable:
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz.tbi"
    log:
        "logs/SortQTLtoolsPhenotypeTable/{Phenotype}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule PhenotypePCs:
    """
    QTLtools format expression PCs as covariates
    including the number of PCs that explain more
    variance then when the phenotype table is
    permuted
    """
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca",
    log:
        "logs/PhenotypePCs/{Phenotype}.log"
    resources:
        mem_mb = 24000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """

rule PlotPhenotypePCs:
    """
    Plot PCs manually to check for outliers
    """
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca.pdf"
    conda:
        "../envs/r_essentials.yml"
    log:
        "logs/PlotPhenotypePCs/{Phenotype}.log"
    shell:
        """
        Rscript {input} {output} &> {log}
        """


rule GetSamplesVcfByChrom:
    input:
        vcf = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz.tbi",
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
    output:
        vcf = "QTLs/QTLTools/{Phenotype}/Genotypes/{chrom}.vcf.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/Genotypes/{chrom}.vcf.gz.tbi"
    resources:
        mem_mb = 8000
    log:
        "logs/GetSamplesVcfByChrom/{Phenotype}/{chrom}.log"
    params:
        filters = "-q 0.01:minor -c 3:minor"
    shell:
        """
        (bcftools view --force-samples -S <(zcat {input.bed} | head -1 | scripts/transpose | awk 'NR>6') -O z {params} {input.vcf} > {output.vcf} ) &> {log}
        tabix -p vcf {output.vcf}
        """

rule GatherWholeGenomeVcfByPhenotype:
    input:
        vcf = expand("QTLs/QTLTools/{{Phenotype}}/Genotypes/{chrom}.vcf.gz", chrom=autosomes),
        tbi = expand("QTLs/QTLTools/{{Phenotype}}/Genotypes/{chrom}.vcf.gz.tbi", chrom=autosomes)
    output:
        vcf = "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz.tbi"
    log:
        "logs/GatherWholeGenomeVcfByPhenotype/{Phenotype}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        bcftools concat -O z {input.vcf} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf}
        """

# bcftools view --force-samples -S <(zcat QTLs/QTLTools/chRNA.IR/OnlyFirstReps.sorted.qqnorm.bed.gz | head -1 | scripts/transpose | awk 'NR>6') -O z -q 0.01:minor -c 3:minor Genotypes/1KG_GRCh38/21.vcf.gz | vcftools --gzvcf - --hwe 1E-7 --stdout

def GetQTLtoolsVcf(wildcards):
    #Because can reuse vcf for some phenotypes
    if wildcards.Phenotype == "polyA.Splicing":
        return "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz"
    else:
        return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz"

def GetQTLtoolsVcfTbi(wildcards):
    #Because can reuse vcf for some phenotypes
    if wildcards.Phenotype == "polyA.Splicing":
        return "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz.tbi"
    else:
        return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz.tbi"

# def GetQTLtoolsWindowSize(wildcards):
#     if wildcards.Phenotype in ChromatinProfilingPhenotypes:
#         return

def GetQTLtoolsFlags(wildcards):
    if wildcards.Phenotype in ["polyA.Splicing", "chRNA.Splicing"]:
        return "--grp-best"
    else:
        return ""

rule QTLtools_cis_permutation_pass:
    input:
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        bed_tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz.tbi",
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        "QTLs/QTLTools/{Phenotype}/PermutationPassChunks/{n}.txt"
    log:
        "logs/QTLtools_cis_permutation_pass/{Phenotype}/{n}.log"
    params:
        Flags = GetQTLtoolsFlags
    shell:
        """
        QTLtools_1.2_CentOS7.8_x86_64 cis --chunk {wildcards.n} {N_PermutationChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.Flags} --permute 1000 --window 100000 &> {log}
        """

rule Gather_QTLtools_cis_permutation_pass:
    input:
        expand( "QTLs/QTLTools/{{Phenotype}}/PermutationPassChunks/{n}.txt", n=range(0, 1+N_PermutationChunks) )
    output:
        "QTLs/QTLTools/{Phenotype}/PermutationPass.txt.gz"
    log:
        "logs/Gather_QTLtools_cis_permutation_pass/{Phenotype}.log"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """

rule AddQValueToPermutationPass:
    input:
        "QTLs/QTLTools/{Phenotype}/PermutationPass.txt.gz"
    output:
        table = "QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz",
        PvalPlot = "QTLs/QTLTools/{Phenotype}/PermutationPass.Pvals.pdf"
    log:
        "logs/AddQValueToPermutationPass/{Phenotype}.log"
    conda:
        "../envs/r_essentials.yml"
    priority:
        10
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} {output.PvalPlot} &> {log}
        """
