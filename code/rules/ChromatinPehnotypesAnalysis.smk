def GetMaxNumFeaturesForQTLCalling(wildcards):
    if wildcards.Phenotype == "CTCF":
        return 50000
    else:
        return 1000000

rule PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks:
    input:
        "featureCounts/{Phenotype}/Counts.txt"
    output:
        AllSamples = "QTLs/QTLTools/{Phenotype}/AllReps.qqnorm.bed.gz",
        FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        Phenotype = "|".join([p for p in ChromatinProfilingPhenotypes if p != 'H3K36ME3'])
    params:
        max_features = GetMaxNumFeaturesForQTLCalling
    log:
        "logs/PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks/{Phenotype}.log"
    shell:
        """
        Rscript scripts/PreparePhenotypeTablesFromFeatureCounts_ChromatinProfilingPeaks.R {input} {params.max_features} {output.AllSamples} {output.FirstReps} {wildcards.Phenotype} &> {log}
        """


use rule PrepareQTLToolsPhenotypeTable_FromFeatureCountsPeaks as MakeChromatinSpliceSitePhenotypes with:
    wildcard_constraints:
        Phenotype = 'H3K4ME1.5PrimeSS|H3K4ME1.3PrimeSS|H3K4ME3.5PrimeSS|H3K4ME3.3PrimeSS|H3K27AC.5PrimeSS|H3K27AC.3PrimeSS|H3K36ME3.5PrimeSS|H3K36ME3.3PrimeSS'

# rule GetGenotypePCs:
#     input:
#         FirstReps = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
#         vcf = "Genotypes/1KG_GRCh38/7.vcf.gz"

H3K36ME3_IndIDs = Fastq_samples.loc[(Fastq_samples['Phenotype']=='H3K36ME3') & (Fastq_samples['RepNumber']==1)]['IndID'].tolist()
rule PrepareH3K36ME3PhenotypeTable:
    input:
        bams = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam", IndID = H3K36ME3_IndIDs),
        bais = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam.bai", IndID = H3K36ME3_IndIDs),
        bed = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        "MiscCountTables/H3K36ME3.bed"
    log:
        "logs/PrepareH3K36ME3PhenotypeTable.log"
    params:
        header = '\t'.join(["#Chr", "start", "end", "pid", "gid", "strand"] + H3K36ME3_IndIDs) + '\n'
    shell:
        """
        printf '{params.header}' > {output}
        (bedtools multicov -bed {input.bed} -bams {input.bams} >> {output}) 2> {log}
        """
        
use rule PrepareH3K36ME3PhenotypeTable as PrepareH3K36ME3NonCodingPhenotypeTable with:
    input:
        bams = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam", IndID = H3K36ME3_IndIDs),
        bais = expand("Alignments/Hisat2_Align/H3K36ME3/{IndID}.1.wasp_filterd.markdup.sorted.bam.bai", IndID = H3K36ME3_IndIDs),
        bed = "NonCodingRNA/annotation/NonCodingRNA.bed.gz"
    output:
        "MiscCountTables/H3K36ME3_ncRNA.bed"
    log:
        "logs/PrepareH3K36ME3_ncRNAPhenotypeTable.log"

rule PreparePhenotypeTableForQTLToolsFromWellFormatedCountTable:
    input:
        "MiscCountTables/{Phenotype}.bed"
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        Phenotype = "H3K36ME3|H3K36ME3_ncRNA"
    log:
        "logs/PreparePhenotypeTableForQTLToolsFromWellFormatedCountTable/{Phenotype}.log"
    shell:
        """
        Rscript scripts/PreparePhenotypeTableFromWellFormattedTable.R {input} {output} &> {log}
        """

rule CalculatePeaksClosestToTSS:
    """
    take annotated TSSs (there may be multiple per gene) and find closest
    existing peak call for which we called QTLs. Also do this after using
    bedtools shuffle as a control.
    """
    input:
        chromatin_peaks = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz",
        TSS = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        real = "Misc/PeaksClosestToTSS/{Phenotype}_real.bed.gz",
        permuted = "Misc/PeaksClosestToTSS/{Phenotype}_permuted.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["H3K4ME3", "H3K27AC", "H3K4ME1"])
    shell:
        """
        awk -v OFS='\\t' '{{print "chr"$1,$2,$3, $6, ".", $4}}' {input.TSS} | sort | uniq  |  bedtools sort -i - | bedtools closest -a - -b <(zcat {input.chromatin_peaks} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,$4}}' | bedtools shuffle -i - -g {input.fai} | bedtools sort -i -  ) -D a | gzip - > {output.permuted}
        awk -v OFS='\\t' '{{print "chr"$1,$2,$3, $6, ".", $4}}' {input.TSS} | sort | uniq  |  bedtools sort -i - | bedtools closest -a - -b <(zcat {input.chromatin_peaks} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,$4}}' ) -D a | gzip - >  {output.real}
        """

rule AssignPeaksToTSS:
    """
    Based on manually inspecting the distance between peaks and TSS
    """
    input:
        real = "Misc/PeaksClosestToTSS/{Phenotype}_real.bed.gz",
    output:
        mappings = "Misc/PeaksClosestToTSS/{Phenotype}_assigned.tsv.gz"
    shell:
        """
        zcat {input.real} | awk -F '\\t' -v OFS='\\t' 'BEGIN {{print "chrom", "TSS_start", "gene", "strand", "peak", "distance"}} $NF<=500 && $NF>=-500 {{print $1, $2, $4, $6, $10, $NF}}' | gzip - > {output}
        """

rule GatherPeaksClosestToTSS:
    input:
        expand("Misc/PeaksClosestToTSS/{Phenotype}_assigned.tsv.gz", Phenotype=["H3K4ME3", "H3K27AC", "H3K4ME1"])

