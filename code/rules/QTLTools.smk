rule SortQTLtoolsPhenotypeTable:
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz.tbi"
    log:
        "logs/SortQTLtoolsPhenotypeTable/{Phenotype}/{StandardizedOrUnstandardized}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        (bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}) &> {log}
        (tabix -p bed {output.bed}) &>> {log}
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
        mem_mb = 16000
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
        elif wildcards.Phenotype == "polyA.Splicing.Subset_EUR":
            return "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz"
        else:
            return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz"
    elif wildcards.QTLsGenotypeSet == "EUR":
        return "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz"
    else:
        return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz"

def GetQTLtoolsVcfTbi(wildcards):
    #Because can reuse vcf for some phenotypes
    if wildcards.QTLsGenotypeSet == "":
        if wildcards.Phenotype == "polyA.Splicing":
            return "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz.tbi"
        elif wildcards.Phenotype == "polyA.Splicing.Subset_YRI":
            return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz.tbi"
        elif wildcards.Phenotype == "polyA.Splicing.Subset_EUR":
            return "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz.tbi"
        else:
            return "QTLs/QTLTools/{Phenotype}/Genotypes/WholeGenome.vcf.gz.tbi"
    elif wildcards.QTLsGenotypeSet == "EUR":
        return "Genotypes/1KG_GRCh38_SubsetEUR/WholeGenome.vcf.gz.tbi"
    else:
        return "Genotypes/1KG_GRCh38_SubsetYRI/WholeGenome.vcf.gz.tbi"

def GetQTLtoolsBed(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == "":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForColoc.sorted.qqnorm.bed.gz"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForGWASColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForGWASColoc.sorted.qqnorm.bed.gz"

def GetQTLtoolsBedTbi(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == "":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz.tbi"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForColoc.sorted.qqnorm.bed.gz.tbi"
    elif wildcards.FeatureCoordinatesRedefinedFor == "ForGWASColoc":
        return "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForGWASColoc.sorted.qqnorm.bed.gz.tbi"


def GetQTLtoolsWindowFlag(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor in ["ForColoc", "ForGWASColoc" ]:
        # using --window 0 sometimes results in errors (exit code 139). No idea why
        return "--window 1"
    elif wildcards.Pass in ["PermutationPass500kb",  "NominalPass500kb"]:
        return "--window 500000"
    elif wildcards.Pass == "PermutationPass250kb":
        return "--window 250000"
    else:
        if wildcards.Phenotype in ["polyA.Splicing", "chRNA.Splicing", "polyA.Splicing.Subset_YRI", "polyA.Splicing.Subset_EUR", 
        "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.RNA.Editing", "chRNA.Splicing.Order"]:
            return "--window 10000"
        elif wildcards.Phenotype.split('.')[-1] in ['5PrimeSS', '3PrimeSS']:
            return "--window 0"
        elif ('5PrimeSS' in wildcards.Phenotype) or ('3PrimeSS' in wildcards.Phenotype):
            return "--window 0"
        elif wildcards.Phenotype in ["chRNA.IR", "polyA.IR", "chRNA.IER", "polyA.IER", "polyA.IER.Subset_YRI",
        "MetabolicLabelled.30min.IER", "MetabolicLabelled.60min.IER"]:
            return "--window 10000"
        else:
            return "--window 100000"

def GetQTLtoolsOtherFlags(wildcards):
    if wildcards.Phenotype in ["polyA.Splicing", "chRNA.Splicing", "polyA.Splicing.Subset_YRI", "polyA.Splicing.Subset_EUR", 
    "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.RNA.Editing", "chRNA.Splicing.Order", "APA_Nuclear", "APA_Total"] and wildcards.Pass == "GroupedPermutationPass":
        return "--grp-best"
    else:
        return ""


def GetQTLtoolsPassFlags(wildcards):
    if wildcards.Pass in ["PermutationPass", "GroupedPermutationPass", "PermutationPass500kb", "PermutationPass250kb"]:
        return "--permute 1000"
    elif wildcards.Pass in ["NominalPass", "NominalPass500kb"]:
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
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz",
        bed_tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}{FeatureCoordinatesRedefinedFor}.sorted.qqnorm.bed.gz.tbi",
        cov = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca"
    output:
        temp("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}Chunks/{QTLTools_chunk_n}{StandardizedOrUnstandardized}.txt")
    log:
        "logs/QTLtools_cis_permutation_pass/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}/{StandardizedOrUnstandardized}{QTLTools_chunk_n}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    envmodules:
        "gsl/2.5"
    params:
        WindowFlag = GetQTLtoolsWindowFlag,
        OtherFlags = GetQTLtoolsOtherFlags,
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = GetExcludeFile,
    wildcard_constraints:
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "PermutationPass|GroupedPermutationPass|PermutationPass500kb|PermutationPass250kb|NominalPass|NominalPass500kb"
    shell:
        """
        {config[QTLtools]} cis --std-err --chunk {wildcards.QTLTools_chunk_n} {N_PermutationChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.OtherFlags} {params.WindowFlag} {params.PassFlags} {params.ExcFlag} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        """


rule Gather_QTLtools_cis_pass:
    input:
        expand( "QTLs/QTLTools/{{Phenotype}}/{{Pass}}{{QTLsGenotypeSet}}{{FeatureCoordinatesRedefinedFor}}Chunks/{QTLTools_chunk_n}{{StandardizedOrUnstandardized}}.txt", QTLTools_chunk_n=ChunkNumbers )
    output:
        "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}{StandardizedOrUnstandardized}.txt.gz"
    log:
        "logs/Gather_QTLtools_cis_pass/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.{StandardizedOrUnstandardized}.log"
    wildcard_constraints:
        Pass = "PermutationPass|GroupedPermutationPass|PermutationPass500kb|PermutationPass250kb|NominalPass|NominalPass500kb"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        

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
        Pass = "PermutationPass|GroupedPermutationPass|PermutationPass500kb|PermutationPass250kb"
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.R {input} {output.table} {wildcards.Pass} &> {log}
        """


rule MakePhenotypeTableToColocPeaksWithGenes:
    input:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        genes = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForColoc.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForColoc.sorted.qqnorm.bed.gz.tbi"
    params:
        #max distance between gene and peaks to attempt coloc
        cis_window = 100000,
        bedtools_intersect_params = "",
        #coloc window from gene
        coloc_window = 100000
    log:
        "logs/MakePhenotypeTableToColocPeaksWithGenes/{Phenotype}.{StandardizedOrUnstandardized}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        Phenotype = "H3K4ME3|H3K27AC|H3K4ME1|CTCF|ProCap|chRNA.Expression_ncRNA"
    shell:
        """
        cat <(zcat {input.bed} | head -1) <(bedtools slop -i {input.genes} -b {params.cis_window} -g {input.fai} | bedtools intersect -b {input.bed} -a - -sorted -wo {params.bedtools_intersect_params} | awk -F'\\t' -v OFS='\\t' '{{$10=$10":"$4; $11=$4; print $0}}' | rev | cut -d$'\\t' -f2- | rev | cut -d$'\\t' -f 4,7- | sort | join -t$'\\t' <(awk -F'\\t' -v OFS='\\t' '{{print $4, $0}}' {input.genes} | sort) - | cut -d$'\\t' -f 2-4,11- | bedtools slop -i - -b {params.coloc_window} -g {input.fai} ) | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocIntronsWithGenes with:
    wildcard_constraints:
        Phenotype = "MetabolicLabelled.30min.IRjunctions|MetabolicLabelled.60min.IRjunctions|polyA.IRjunctions|chRNA.IRjunctions|MetabolicLabelled.30min.IER|MetabolicLabelled.60min.IER|polyA.IER|chRNA.IER|MetabolicLabelled.30min.IR|MetabolicLabelled.30min.Splicing|MetabolicLabelled.60min.IR|MetabolicLabelled.60min.Splicing|chRNA.IR|chRNA.Splicing|polyA.Splicing|polyA.IR|polyA.Splicing.Subset_YRI|polyA.IER.Subset_YRI|polyA.IR.Subset_YRI|chRNA.Slopes|chRNA.Splicing.Order|chRNA.RNA.Editing|APA_Nuclear|APA_Total"
    params:
        cis_window = 0,
        bedtools_intersect_params = "-s",
        coloc_window = 100000,

use rule MakePhenotypeTableToColocPeaksWithGenes as MakePhenotypeTableToColocGenes with:
    wildcard_constraints:
        Phenotype = "Expression.Splicing.Subset_YRI|chRNA.Expression.Splicing|Expression.Splicing|MetabolicLabelled.30min|MetabolicLabelled.60min|H3K36ME3|H3K36ME3_ncRNA"
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
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}.sorted.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        loci = "gwas_summary_stats/LeadSnpWindows.bed"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForGWASColoc.sorted.qqnorm.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps{StandardizedOrUnstandardized}ForGWASColoc.sorted.qqnorm.bed.gz.tbi"
    resources:
        mem_mb = 48000
    log:
        "logs/MakePhenotypeTableToColocFeaturesWithGWASLoci/{Phenotype}.{StandardizedOrUnstandardized}.log"
    shell:
        """
        (cat <(zcat {input.bed} | head -1) <(  bedtools intersect  -wo -a {input.bed} -b {input.loci} -sorted | awk -F'\\t' -v OFS='\\t' '{{$4=$4":"$(NF-1); $5=$(NF-1); $2=$(NF-3); $3=$(NF-2); print $0}}' | rev | cut -f 6- | rev ) | tr ' ' '\\t' | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

rule tabixNominalPassQTLResults:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}{StandardizedOrUnstandardized}.txt.gz"
    wildcard_constraints:
        Pass = "NominalPass"
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
        # sort_temp = ""
    output:
        txt = "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}{StandardizedOrUnstandardized}.txt.tabix.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}{StandardizedOrUnstandardized}.txt.tabix.gz.tbi"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "logs/tabixNominalPassQTLResults/{Phenotype}.{Pass}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.{StandardizedOrUnstandardized}.log"
    shadow: "shallow"
    shell:
        """
        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
        """



#rule tabixNominalPass:
#    """
#    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
#    for easy access with tabix
#    """
#    input:
#        "QTLs/QTLTools/{Phenotype}/NominalPass.txt.gz"
#    wildcard_constraints:
#        Phenotype = "Expression.Splicing"
#    params:
#        sort_temp = '-T ' + config['scratch'][:-1]
#        # sort_temp = ""
#    output:
#        txt = "QTLs/QTLTools/{Phenotype}/NominalPass.txt.tabix.gz",
#        tbi = "QTLs/QTLTools/{Phenotype}/NominalPass.txt.tabix.gz.tbi"
#    resources:
#        mem_mb = much_more_mem_after_first_attempt
#    log:
#        "logs/tabixNominalPass/{Phenotype}.log"
#    shadow: "shallow"
#    shell:
#        """
#        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
#        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
#        """


rule SampleRandomPvals:
    input:
        txt = "QTLs/QTLTools/{Phenotype}/NominalPass{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.tabix.gz",
    output:
        "QTLs/QTLTools/{Phenotype}/NominalPass{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.RandomSamplePvals.txt.gz"
    log:
        "logs/SampleRandomPvals/{Phenotype}.{QTLsGenotypeSet}.{FeatureCoordinatesRedefinedFor}.log"
    params:
        NumVals = 1000000
    shell:
        """
        (zcat {input.txt} | awk -F'\\t' 'NR>1 {{print $12}}' | shuf -n {params.NumVals} | gzip - > {output} ) &> {log}
        """

rule GatherSampleRandomPvals:
    input:
        expand("QTLs/QTLTools/{Phenotype}/NominalPass{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.RandomSamplePvals.txt.gz", Phenotype=PhenotypesToColoc, QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc")

rule CountPCs_and_SampleSize:
    input:
        beds = expand("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz", Phenotype=MyPhenotypes),
        PCs = expand("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.pca", Phenotype=MyPhenotypes)
    output:
        N = "QC/SampleSize.PerPhenotype.tsv",
        PCs = "QC/NumPCs.PerPhenotype.tsv"
    shell:
        """
        set +o pipefail;
        for fn in {input.beds};
        do
            printf "$fn\\tn\\t" >> {output.N}
            zcat $fn | head -1 | tr '\\t' '\\n' | wc -l | awk '{{print $1-6}}' >> {output.N}
        done
        for fn in {input.PCs};
        do
            printf "$fn\\tn_PCs\\t" >> {output.PCs}
            wc -l $fn | awk '{{print $1-1}}' >> {output.PCs}
        done
        """

rule Gather_NominalPassTabixForUnstandardizedScatters:
    input:
        expand("QTLs/QTLTools/{Phenotype}/NominalPassForColocUnstandardized.txt.tabix.gz.tbi", Phenotype=["chRNA.Expression.Splicing", "Expression.Splicing", "polyA.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min"]),
        expand("QTLs/QTLTools/{Phenotype}/NominalPassForColocUnstandardized.txt.tabix.gz", Phenotype=["chRNA.Expression.Splicing", "Expression.Splicing", "polyA.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
