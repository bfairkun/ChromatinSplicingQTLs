
rule ExtractJuncs:
    input:
        bam = GetBamForBigwig,
        bai = GetBaiForBigwig
    output:
        junc = temp(expand("SplicingAnalysis/leafcutter/juncfiles/chr{chrom}/{{Phenotype}}_{{IndID}}_{{Rep}}.junc", chrom=autosomes)),
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",
    params:
        # strand = GetLibStrandForRegtools
        strand = 0
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{Phenotype}/{IndID}.{Rep}.log"
    shell:
        """
        for chrom in {autosomes}
        do
            (regtools junctions extract -m 20 -s {params.strand} -r chr${{chrom}} {input.bam} > SplicingAnalysis/leafcutter/juncfiles/chr${{chrom}}/{wildcards.Phenotype}_{wildcards.IndID}_{wildcards.Rep}.junc ) &> {log}
        done
        cat {output.junc} > {output.junc_autosomes}
        """

rule make_leafcutter_juncfile:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
    output:
        "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    params:
        SamplesToRemove = ""
    run:
        import os
        if params.SamplesToRemove:
            SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        else:
            SamplesToRemove=[]
        with open(output[0], "w") as out:
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter_cluster:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_cluster/autosomes.log"
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input.juncfile_list} -r SplicingAnalysis/leafcutter/clustering/autosomes/ &> {log}
        """

rule annotate_juncfiles:
    input:
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        basic_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",
    output:
        basic = "SplicingAnalysis/leafcutter/regtools_annotate/basic/{Phenotype}_{IndID}_{Rep}.bed.gz",
        comprehensive = "SplicingAnalysis/leafcutter/regtools_annotate/comprehensive/{Phenotype}_{IndID}_{Rep}.bed.gz"
    log:
        "logs/annotate_juncfiles/{Phenotype}_{IndID}_{Rep}.log"
    shell:
        """
        (regtools junctions annotate {input.junc_autosomes} {input.fa} {input.basic_gtf} | gzip - > {output.basic} ) &> {log}
        (regtools junctions annotate {input.junc_autosomes} {input.fa} {input.Comprehensive_gtf} | gzip - > {output.comprehensive} ) &>> log
        """

rule GetIntronFeatures:
    """
    independent introns that don't overlap any annotated exons, and are in
    expressed genes to quantify intron retention phenotypes
    """
    input:
        AllIntron = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.introns.bed",
        IndependentIntrons = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.independentintron.bed",
        eQTL_genes_qqnorm = "QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",

    log:
        "logs/GetIntronFeatures.log"
    output:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.saf"
    shell:
        """
        (bedtools sort -i {input.AllIntron} | bedtools intersect -a {input.IndependentIntrons} -b - -f 1 -r  -sorted -wo | awk -F'\\t' '$4==$10 {{print $1, $2, $3, $4, $5, $11}}' | sort | uniq | grep -F -f <(zcat {input.eQTL_genes_qqnorm} | awk -F'\\t' '{{print $4}}') - | awk -v OFS='\\t' '{{$1="chr"$1; print $0}}' |  bedtools slop -b -3 -i - -g {input.fai} | awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }} {{print $4"_IntID."NR, $1, $2, $3, $6}}' > {output} ) &> {log}
        """

rule GetIntronRetentionSpliceSitesAndIntronsBed:
    input:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.saf"
    output:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.bed"
    shell:
        """
        awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $2, $3,$4,$1, ".", $5}}' {input} > {output}
        """

rule CountSpliceSitesOverIR_Intron_Features:
    input:
        IRFeatures = "SplicingAnalysis/Annotations/IntronRetentionTargets.bed",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        vcf = "QTLs/QTLTools/Expression.Splicing.Subset_YRI/Genotypes/WholeGenome.vcf.gz"
    output:
        ss_bed = "SplicingAnalysis/Annotations/IntronRetentionTargets_ss.bed",
        vcf = "SplicingAnalysis/Annotations/YRI.Snps.OverInts.vcf.gz"
    shell:
        """
        bedtools slop -b 1 -i {input.IRFeatures} -g {input.fai} | bedtools flank -b 2 -g {input.fai} -i - | bedtools sort -i - > {output.ss_bed}
        bcftools view -O z -R {output.ss_bed} {input.vcf} > {output.vcf}
        """

rule featureCounts_IR:
    input:
        bam = GetBamForPhenotype,
        annotations = "SplicingAnalysis/Annotations/IntronRetentionTargets.saf",
    output:
        "SplicingAnalysis/IR/{Phenotype}/Counts.txt"
    threads:
        8
    params:
        extraParams = GetFeatureCountsParams
    resources:
        mem = 12000,
        cpus_per_node = 9
    log:
        "logs/featureCounts_IR/{Phenotype}.log"
    shell:
        """
        featureCounts -p {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

def Get_intron_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "SplicingAnalysis/IR/chRNA.Expression.Splicing/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "SplicingAnalysis/IR/Expression.Splicing/Counts.txt"

def Get_gene_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "featureCounts/chRNA.Expression/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "featureCounts/polyA.Expression/Counts.txt"

def SampleMinimumIR_Counts(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "100000"
    elif wildcards.Phenotype == "polyA.IR":
        return "1"

rule featureCounts_IR_to_bedgz:
    input:
        IR_counts = Get_intron_feature_Counts_ForIR,
        Gene_counts = Get_gene_feature_Counts_ForIR
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
        plot = "SplicingAnalysis/IR/{Phenotype}/IR.Ratio.plot.pdf"
    wildcard_constraints:
        Phenotype = "|".join(["chRNA.IR", "polyA.IR"])
    log:
        "logs/featureCounts_IR_to_bedgz/{Phenotype}.log"
    params:
        SampleMinimumIR_Counts = SampleMinimumIR_Counts
    resources:
        mem_mb = 16000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/ProcessIRFeatureCounts.R {input.IR_counts} {input.Gene_counts} {output.bed} {output.plot} {params.SampleMinimumIR_Counts} &> {log}
        """

rule SplitLeafcutter_countsTable:
    input:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz"
    output:
        expand("SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.{Phenotype}.gz", Phenotype=RNASeqPhenotypes)
    log:
        "logs/SplitLeafcutter_countsTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/SplitLeafcutterPerindCounts.R {input} &> {log}
        """

def GetSplitLeafcutterCountTablesForPhenotype(wildcards):
    if wildcards.Phenotype == "polyA.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.Expression.Splicing.gz"
    elif wildcards.Phenotype == "chRNA.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.chRNA.Expression.Splicing.gz"


rule leafcutter_PreparePhenotypes:
    """
    scripts in main leafcutter repo is buggy. Use fork from https://github.com/mdshw5/leafcutter.git
    """
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing"])
    conda:
        "../envs/py27.yaml"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_PreparePhenotypes/{Phenotype}.log"
    shadow: "shallow"
    shell:
        """
        python2.7 scripts/prepare_phenotype_table.py {input} -p 20 &> {log}
        (awk -F'\\t' -v OFS='\\t' 'NR==1 {{$4="pid\\tgid\\tstrand"; print $0}} FNR!=1 {{$1="chr"$1; split($4, a, ":"); split(a[4], b, "_"); $4=$4"\\t"$1"_"a[4]"\\t"b[3]; print $0}}' {input}.qqnorm_chr* | gzip - > {output} ) &>> {log}
        """

rule Subset_YRI_leafcutter_phenotype_table:
    input:
        "QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/polyA.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"
    shell:
        """
        python scripts/subsample_polyA.Splicing_YRI.py
        """

# rule ScoreSpliceSiteSNPs:
#     input:
#         vcf = ,
#         introns = ,
#     output:
#     conda:
#     shell:
