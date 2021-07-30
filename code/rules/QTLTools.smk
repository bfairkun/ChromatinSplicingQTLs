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
        return "--grp-best --window 10000"
    elif wildcards.Phenotype in ["chRNA.IR", "polyA.IR", "chRNA.Expression.Splicing"]:
        return "--window 10000"
    else:
        return "--window 100000"

rule QTLtools_cis_permutation_pass:
    input:
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        bed_tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz.tbi",
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        temp("QTLs/QTLTools/{Phenotype}/PermutationPassChunks/{n}.txt")
    log:
        "logs/QTLtools_cis_permutation_pass/{Phenotype}/{n}.log"
    params:
        Flags = GetQTLtoolsFlags
    shell:
        """
        QTLtools_1.2_CentOS7.8_x86_64 cis --chunk {wildcards.n} {N_PermutationChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.Flags} --permute 1000  &> {log}
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
    log:
        "logs/AddQValueToPermutationPass/{Phenotype}.log"
    conda:
        "../envs/r_essentials.yml"
    priority:
        10
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} &> {log}
        """


rule MakePhenotypeTableToColocPeaksWithGenes:
    input:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        genes = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz.tbi"
    params:
        #max distance between gene and peaks to attempt coloc
        cis_window = 100000,
        bedtools_intersect_params = "",
        #coloc window from gene
        coloc_window = 100000
    log:
        "logs/MakePhenotypeTableToColocPeaksWithGenes/{Phenotype}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        Phenotype = "H3K4ME3|H3K27AC|H3K4ME1|CTCF"
    shell:
        """
        cat <(zcat {input.bed} | head -1) <(bedtools slop -i {input.genes} -b {params.cis_window} -g {input.fai} | bedtools intersect -b {input.bed} -a - -sorted -wo {params.bedtools_intersect_params} | awk -F'\\t' -v OFS='\\t' '{{$10="{wildcards.Phenotype}:"$4":"$10; $11=$4; print $0}}' | rev | cut -d$'\\t' -f2- | rev | cut -d$'\\t' -f 4,7- | sort | join -t$'\\t' <(awk -F'\\t' -v OFS='\\t' '{{print $4, $0}}' {input.genes} | sort) - | cut -d$'\\t' -f 2-4,11- | bedtools slop -i - -b {params.coloc_window} -g {input.fai} ) | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocIntronsWithGenes with:
    wildcard_constraints:
        Phenotype = "chRNA.IR|chRNA.Splicing|polyA.Splicing|polyA.IR"
    params:
        cis_window = 0,
        bedtools_intersect_params = "-s",
        coloc_window = 100000,

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocGenes with:
    wildcard_constraints:
        Phenotype = "Expression.Splicing|chRNA.Expression.Splicing"
    params:
        cis_window = 0,
        bedtools_intersect_params = "-s -f 1 -r",
        coloc_window = 100000,

rule QTLtools_cis_nominal_pass_for_coloc:
    input:
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz",
        bed_tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz.tbi",
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        "QTLs/QTLTools/{Phenotype}/NominalPass_ForColoc.txt.gz",
    log:
        "logs/QTLtools_cis_nominal_pass_for_coloc/{Phenotype}.log"
    params:
        extra = ""
        # extra = "--chunk 1 500"
    shell:
        """
        QTLtools_1.2_CentOS7.8_x86_64 cis  --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out QTLs/QTLTools/{wildcards.Phenotype}/NominalPass_ForColoc.txt  --nominal 1 --window 0 {params.extra} &> {log}
        gzip QTLs/QTLTools/{wildcards.Phenotype}/NominalPass_ForColoc.txt
        """

use rule QTLtools_cis_permutation_pass as QTLtools_cis_permutation_pass_for_coloc with:
    input:
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz",
        bed_tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz.tbi",
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        temp("QTLs/QTLTools/{Phenotype}/PermutationPassForColocChunks/{n}.txt")
    log:
        "logs/QTLtools_cis_permutation_passForColoc/{Phenotype}/{n}.log"
    params:
        Flags = "--window 0"

use rule Gather_QTLtools_cis_permutation_pass as Gather_QTLtools_cis_permutation_pass_ForColoc with:
    input:
        expand( "QTLs/QTLTools/{{Phenotype}}/PermutationPassForColocChunks/{n}.txt", n=range(0, 1+N_PermutationChunks) )
    output:
        "QTLs/QTLTools/{Phenotype}/PermutationPassForColoc.txt.gz"
    log:
        "logs/Gather_QTLtools_cis_permutation_pass_ForColoc/{Phenotype}.log"

# rule ShortenQTLToolsOutput:
#     input:
#         "QTLs/QTLTools/{Phenotype}/NominalPass_ForColoc.txt.gz",
#     output:

# use rule QTLtools_cis_permutation_pass as QTLtools_cis_nominal_pass_for_coloc with:
