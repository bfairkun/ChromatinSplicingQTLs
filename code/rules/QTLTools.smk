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

rule MakeUniversalVcfForTesting:
    """
    In my QTL calling, I consider SNPs for each molecular phenotype based a MAF cutoff based on the samples for the molecular assay. Alternatively, I consider a universal set of SNPs for all molecular assays based on MAF among all YRI samples in 1KG
    """
    input:
        vcf = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz.tbi",
        metadata = "../data/igsr_samples.tsv.gz"
    output:
        vcf = "Genotypes/1KG_GRCh38_SubsetYRI/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38_SubsetYRI/{chrom}.vcf.gz.tbi"
    resources:
        mem_mb = 8000
    log:
        "logs/MakeUniversalVcfForTesting/{chrom}.log"
    params:
        filters = "-q 0.01:minor"
    shell:
        """
        (bcftools view --force-samples -S <(zcat {input.metadata} | awk -F'\\t' '$4=="YRI" {{ print $1 }}') -O z {params} {input.vcf} > {output.vcf} ) &> {log}
        tabix -p vcf {output.vcf}
        """


use rule GatherWholeGenomeVcfByPhenotype as GatherWholeGenomeVcfUniversal with:
    input:
        vcf = expand("Genotypes/1KG_GRCh38_SubsetYRI/{chrom}.vcf.gz", chrom=autosomes),
        tbi = expand("Genotypes/1KG_GRCh38_SubsetYRI/{chrom}.vcf.gz.tbi", chrom=autosomes)
    output:
        vcf = "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz.tbi"
    log:
        "logs/GatherWholeGenomeVcfUniversal.log"


def GetQTLtoolsVcf(wildcards):
    #Because can reuse vcf for some phenotypes
    if wildcards.QTLsGenotypeSet == "":
        if wildcards.Phenotype == "polyA.Splicing":
            return "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz"
        elif wildcards.Phenotype == "polyA.Splicing.Subset_YRI":
            return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz"
        else:
            return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz"
    else:
        return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz"

def GetQTLtoolsVcfTbi(wildcards):
    #Because can reuse vcf for some phenotypes
    if wildcards.QTLsGenotypeSet == "":
        if wildcards.Phenotype == "polyA.Splicing":
            return "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz.tbi"
        elif wildcards.Phenotype == "polyA.Splicing.Subset_YRI":
            return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz.tbi"
        else:
            return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz.tbi"
    else:
        return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz.tbi"

def GetQTLtoolsBed(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == "":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForGWASColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForGWASColoc.sorted.qqnorm.bed.gz"

def GetQTLtoolsBedTbi(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == "":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz.tbi"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForColoc.sorted.qqnorm.bed.gz.tbi"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForGWASColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForGWASColoc.sorted.qqnorm.bed.gz.tbi"

def GetQTLtoolsFlags(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor in ["ForColoc", "ForGWASColoc" ]:
        # using --window 0 sometimes results in errors (exit code 139). No idea why
        return "--window 1"
    else:
        if wildcards.Phenotype in ["polyA.Splicing", "chRNA.Splicing", "polyA.Splicing.Subset_YRI"]:
            return "--grp-best --window 10000"
        elif wildcards.Phenotype in ["chRNA.IR", "polyA.IR", "chRNA.Expression.Splicing"]:
            return "--window 10000"
        else:
            return "--window 100000"


def GetQTLtoolsPassFlags(wildcards):
    if wildcards.Pass == "PermutationPass":
        return "--permute 1000"
    elif wildcards.Pass == "NominalPass":
        return "--nominal 1"
        
def GetExcludeFile(wildcards):
    if wildcards.Phenotype.split('.')[0] == 'chRNA':
        return '--exclude-samples config/chRNA.exc'
    else:
        return ''

rule QTLtools_generalized:
    input:
        vcf = GetQTLtoolsVcf,
        tbi = GetQTLtoolsVcfTbi,
        bed = GetQTLtoolsBed,
        bed_tbi = GetQTLtoolsBedTbi,
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}Chunks/{n}.txt"
    log:
        "logs/QTLtools_cis_permutation_pass/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}/{n}.log"
    resources:
        mem_mb = 58000
    envmodules:
        "gsl/2.5"
    params:
        Flags = GetQTLtoolsFlags,
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = GetExcludeFile,
        N_PermutationChunks = Get_N_PermutationChunks
    shell:
        """
        {config[QTLtools]} cis --std-err --chunk {wildcards.n} {params.N_PermutationChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.Flags} {params.PassFlags} {params.ExcFlag} &> {log}
        """

rule Gather_QTLtools_cis_pass:
    input:
        expand( "QTLs/QTLTools/{{Phenotype}}/{{Pass}}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}Chunks/{n}.txt", n=range(0, 101) )
    output:
        "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz"
    wildcard_constraints:
        Phenotype = "|".join([x for x in MyPhenotypes if ((x != 'chRNA.Expression_eRNA') and (x != 'chRNA.Expression_cheRNA'))])
    log:
        "logs/Gather_QTLtools_cis_pass/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        
use rule Gather_QTLtools_cis_pass as Gather_QTLtools_cis_pass_ncRNA with:
    input:
        expand( "QTLs/QTLTools/{{Phenotype}}/{{Pass}}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}Chunks/{n}.txt", n=range(0, 11) )
    wildcard_constraints:
        Phenotype = "chRNA.Expression_eRNA|chRNA.Expression_cheRNA"

rule AddQValueToPermutationPass:
    input:
        "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz"
    output:
        table = "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.FDR_Added.txt.gz",
    log:
        "logs/AddQValueToPermutationPass/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    conda:
        "../envs/r_essentials.yml"
    priority:
        10
    wildcard_constraints:
        Pass = "PermutationPass"
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} {wildcards.Pass} &> {log}
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
        cat <(zcat {input.bed} | head -1) <(bedtools slop -i {input.genes} -b {params.cis_window} -g {input.fai} | bedtools intersect -b {input.bed} -a - -sorted -wo {params.bedtools_intersect_params} | awk -F'\\t' -v OFS='\\t' '{{$10=$10":"$4; $11=$4; print $0}}' | rev | cut -d$'\\t' -f2- | rev | cut -d$'\\t' -f 4,7- | sort | join -t$'\\t' <(awk -F'\\t' -v OFS='\\t' '{{print $4, $0}}' {input.genes} | sort) - | cut -d$'\\t' -f 2-4,11- | bedtools slop -i - -b {params.coloc_window} -g {input.fai} ) | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocIntronsWithGenes with:
    wildcard_constraints:
        Phenotype = "chRNA.IR|chRNA.Splicing|polyA.Splicing|polyA.IR|polyA.Splicing.Subset_YRI|polyA.IR.Subset_YRI|chRNA.Slopes"
    params:
        cis_window = 0,
        bedtools_intersect_params = "-s",
        coloc_window = 100000,

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocGenes with:
    wildcard_constraints:
        Phenotype = "Expression.Splicing.Subset_YRI|chRNA.Expression.Splicing|Expression.Splicing|MetabolicLabelled.30min|MetabolicLabelled.60min|chRNA.Expression_eRNA|chRNA.Expression_cheRNA|chRNA.Expression_lncRNA|chRNA.Expression_snoRNA"
    params:
        cis_window = 0,
        bedtools_intersect_params = "-s -f 1 -r",
        coloc_window = 100000,

#"|".join(['Expression.Splicing.Subset_YRI', 'chRNA.Expression.Splicing', 'Expression.Splicing', 'MetabolicLabelled.30min', 
#          'MetabolicLabelled.60min', 'polyA.Expression.AllRNA.Subset_YRI', 'MetabolicLabelled.30min.AllRNA.Subset_YRI',
#          'MetabolicLabelled.60min.AllRNA.Subset_YRI', 'chRNA.Expression.AllRNA.Subset_YRI'])

rule MakePhenotypeTableToColocFeaturesWithGWASLoci:
    """
    To colocalize molQTLs with gwas we need summary stats (betas and se) for
    every snp in the same window as a gwas locus. Here, I will take a phenotype
    table, find the features that intersect the gwas locus (1MB window centered
    on lead snp, defined in the input file), and output a phenotype table with
    every intersection with coordinates redefined to be the 1MB gwas locus
    window. Running QTLtools on this will help me get the necessary summary
    stats in those windows
    """
    input:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        loci = "gwas_summary_stats/LeadSnpWindows.bed"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForGWASColoc.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstRepsForGWASColoc.sorted.qqnorm.bed.gz.tbi"
    resources:
        mem_mb = 48000
    log:
        "logs/MakePhenotypeTableToColocFeaturesWithGWASLoci/{Phenotype}.log"
    shell:
        """
        (cat <(zcat {input.bed} | head -1) <(  bedtools intersect  -wo -a {input.bed} -b {input.loci} -sorted | awk -F'\\t' -v OFS='\\t' '{{$4=$4":"$(NF-1); $5=$(NF-1); $2=$(NF-3); $3=$(NF-2); print $0}}' | rev | cut -f 6- | rev ) | tr ' ' '\\t' | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

