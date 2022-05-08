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
    conda:
        "../envs/regtools.yml"
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
        extraParams = GetFeatureCountsParams,
        pairedEndParams = PairedEndParams
    resources:
        mem = 12000,
        cpus_per_node = 9
    log:
        "logs/featureCounts_IR/{Phenotype}.log"
    shell:
        """
        featureCounts {params.pairedEndParams} {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

def Get_intron_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "SplicingAnalysis/IR/chRNA.Expression.Splicing/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "SplicingAnalysis/IR/Expression.Splicing/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.30min.IR":
        return "SplicingAnalysis/IR/MetabolicLabelled.30min/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.IR":
        return "SplicingAnalysis/IR/MetabolicLabelled.60min/Counts.txt"

def Get_gene_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "featureCounts/chRNA.Expression/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "featureCounts/polyA.Expression/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.30min.IR":
        return "featureCounts/MetabolicLabelled.30min/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.IR":
        return "featureCounts/MetabolicLabelled.60min/Counts.txt"

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
        plot = "SplicingAnalysis/IR/{Phenotype}/IR.Ratio.plot.pdf",
        ir = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.IR.bed.gz",
    wildcard_constraints:
        Phenotype = "|".join(["chRNA.IR", "polyA.IR", "MetabolicLabelled.30min.IR", "MetabolicLabelled.60min.IR"])
    log:
        "logs/featureCounts_IR_to_bedgz/{Phenotype}.log"
    params:
        SampleMinimumIR_Counts = "1" #SampleMinimumIR_Counts
    resources:
        mem_mb = 32000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/ProcessIRFeatureCounts.R {input.IR_counts} {input.Gene_counts} {output.bed} {output.plot} {output.ir} {params.SampleMinimumIR_Counts} &> {log}
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
    elif wildcards.Phenotype == "MetabolicLabelled.30min.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.MetabolicLabelled.30min.gz"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.MetabolicLabelled.60min.gz"
        

rule Leafcutter_countsTable_toPSI:
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        temp("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed")
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
    log:
        "logs/Leafcutter_countsTable_toPSI/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/LeafcutterToPSIBed.R {input} {output} &> {log}
        """

rule sortAndIndexPSITable:
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed.gz.tbi"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
    shell:
        """
        bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule leafcutter_PreparePhenotypes:
    """
    scripts in main leafcutter repo is buggy. Use fork from https://github.com/mdshw5/leafcutter.git
    """
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
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

rule Subset_YRI_phenotype_table:
    input:
        input_file = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
        igsr = '../data/igsr_samples.tsv.gz'
    output:
        "QTLs/QTLTools/{Phenotype}.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "polyA.IR"])
    log:
        "logs/Subsample_YRI.{Phenotype}.log"
    shell:
        """
        python scripts/subsample.Splicing_YRI.py --input {input.input_file} --output {output} &> {log}
        """
        
        
def GetSkippedSamples(wildcards):
    if wildcards.Phenotype == 'chRNA.Splicing':
        return '--skip_samples NA18855'
    else:
        return '--skip_samples NOSAMPLE'
        
        
rule GetSpliceSitePhenotypes:
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        "QTLs/QTLTools/{Phenotype}.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/{Phenotype}.5PrimeSS/OnlyFirstReps.PSI.bed.gz",
        "QTLs/QTLTools/{Phenotype}.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/{Phenotype}.3PrimeSS/OnlyFirstReps.PSI.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing"])
    log:
        "logs/leafcutter_PreparePhenotypes/{Phenotype}.SpliceSites.log"
    params:
        flags = GetSkippedSamples
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/GetSpliceSitesFromLeafcutter.py --counts {input} --output QTLs/QTLTools/{wildcards.Phenotype} {params.flags} &> {log}
        """
    
    
# Ugly rule, but merging quick processes to avoid clogging slurm
rule MakeSpliceSiteAnnotation:
    input:
        "QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz"
    output:
        "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed.gz",
        "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed.gz"
    params:
        up = "94",
        dn = "97"
    shell:
        """
        zcat QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' > ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ; 
        zcat QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ;
        zcat QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ; 
        zcat QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ;
        sort -u ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed > ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed;
        gzip ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed;
        rm ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed;
        zcat QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' > ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ; 
        zcat QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ;
        zcat QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ; 
        zcat QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ;
        sort -u ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed > ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed;
        gzip ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed;
        rm ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed;
        """
        
rule MakeChromatinSpliceSitePeakPhenotypes:
    input:
        SpliceSites = "ReferenceGenome/Annotations/SpliceSites.{Prime}.bed.gz",
        Chromatin = "QTLs/QTLTools/{chromatinPhenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{chromatinPhenotype}.{Prime}.Peaks/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'H3K4ME1|H3K4ME3|H3K27AC'
    params:
        start = "85",
        end = "86",
        name = "82",
        ncols = "78"
    shell:
        """
        zcat {input.Chromatin} | head -n 1 > QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed;
        bedtools intersect -wo -F 0.1 -a {input.Chromatin} -b {input.SpliceSites} | awk -F"\\t" '{{printf($1"\\t"${params.start}"\\t"${params.end}"\\t"${params.name}"\\t"); for (x=5; x<={params.ncols}; x++) printf("%s\\t", $x);printf("\\n"); }}' >> QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed;
        gzip QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed
        """
    
use rule MakeChromatinSpliceSitePeakPhenotypes as MakeChromatinSpliceSitePeakH3K36ME3 with:
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'H3K36ME3'
    params:
        start = "107",
        end = "108",
        name = "104",
        ncols = "100"

    
use rule MakeChromatinSpliceSitePeakPhenotypes as MakeChromatinSpliceSiteCTCF with:
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'CTCF'
    params:
        start = "64",
        end = "65",
        name = "61",
        ncols = "57"


rule MakeSpliceSiteSAF:
    input:
        in5 = "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed.gz",
        in3 = "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed.gz"
    output:
        out5 = "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.saf",
        out3 = "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.saf",
    shell:
        """
        zcat {input.in5} | awk -F'\\t' '{{print $4"\\t"$1"\\t"$2"\\t"$3"\\t."}}' > {output.out5};
        zcat {input.in3} | awk -F'\\t' '{{print $4"\\t"$1"\\t"$2"\\t"$3"\\t."}}' > {output.out3};
        """

def  GetAnnotationsForSpliceSite(wildcards):
    if wildcards.SpliceSite == '5PrimeSS':
        return "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.saf"
    elif wildcards.SpliceSite == '3PrimeSS':
        return "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.saf"

rule featureCountsForSpliceSite:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForSpliceSite
    output:
        "featureCounts/{Phenotype}.{SpliceSite}/Counts.txt"
    params:
        extraParams = GetFeatureCountsParams,
        paired = PairedEndParams
    threads:
        8
    wildcard_constraints:
        Phenotype = 'H3K4ME1|H3K4ME3|H3K27AC|H3K36ME3',
        SpliceSite = '5PrimeSS|3PrimeSS'
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.{SpliceSite}.log"
    shell:
        """
        featureCounts {params.paired} {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

rule MakeNormalizedPsiTables:
    input:
        numers = "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    output:
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"])
    log:
        "logs/MakeNormalizedPsiTables.log"
    conda:
        "envs/r_2.yaml"
    shell:
        """
        Rscript scripts/MakeNormalizedPSI.Tables.R SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI. {input.numers} &> {log}
        """

rule BgzipAndTabixPsiTables:
    input:
        "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed"
    output:
        bed = "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz",
        tbi = "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz.tbi"
    log:
        "logs/BgzipAndTabixPsiTables/{Phenotype}.log"
    shell:
        """
        bgzip {input}
        tabix -p bed {output.bed}
        """



    
def GetSPLICEqModeParams(wildcards):
    if wildcards.mode == "IER":
        return "--IERatio"
    else:
        return "-f 3"

rule RunSPLICEq_chRNA:
    input:
        bam = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/chRNA.{mode}/{IndID}.spliceq.tab"
    params:
        GetSPLICEqModeParams
    conda:
        "../envs/spliceq.yml"
    wildcard_constraints:
        IndID = '|'.join(chRNAsamples),
        mode = 'IRjunctions|IER'
    log:
        "logs/spliceq/chRNA.{mode}/{IndID}.log"
    shell:
        """
        SPLICE-q.py -b {input.bam} -g {input.Comprehensive_gtf} -o {output} {params} --MinCoverage 0 &> {log}
        """

use rule RunSPLICEq_chRNA as rule RunSPLICEq_polyA with:
    input:
        bam = 'Alignments/STAR_Align/Expression.Splicing/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/Expression.Splicing/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/polyA.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/polyA.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(polyAsamples),
        mode = 'IRjunctions|IER'
        
use rule RunSPLICEq_chRNA as rule RunSPLICEq_M30 with:
    input:
        bam = 'Alignments/STAR_Align/MetabolicLabelled.30min/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/MetabolicLabelled.30min/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/MetabolicLabelled.30min.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/MetabolicLabelled.30min.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(metabolic30samples),
        mode = 'IRjunctions|IER'
        
use rule RunSPLICEq_chRNA as rule RunSPLICEq_M60 with:
    input:
        bam = 'Alignments/STAR_Align/MetabolicLabelled.60min/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/MetabolicLabelled.60min/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/MetabolicLabelled.60min.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/MetabolicLabelled.60min.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(metabolic60samples),
        mode = 'IRjunctions|IER'
        
rule CollectSPLICEq_chRNA:
    input: 
        expand("SplicingAnalysis/SPLICEq/chRNA.{{mode}}/{IndID}.spliceq.tab", IndID=chRNAsamples)
    output:
        out_qqnorm = "QTLs/QTLTools/{Phenotype}.{mode}/OnlyFirstReps.qqnorm.bed.gz",
        out = "QTLs/QTLTools/{Phenotype}.{mode}/OnlyFirstReps.ratios.bed.gz",
    log:
        "logs/spliceq/{Phenotype}.{mode}.log"
    wildcard_constraints:
        mode = 'IRjunctions|IER',
        Phenotype = 'chRNA'
    params:
        "NA18855"
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/collect_spliceq.py --spliceq_dir SplicingAnalysis/SPLICEq/{wildcards.Phenotype}.{wildcards.mode}/ --spliceq_mode {wildcards.mode} --output {output.out} --output_qqnorm {output.out_qqnorm} --skip_samples {params} &> {log}
        """

use rule CollectSPLICEq_chRNA as CollectSPLICEq_polyA with:
    input: 
        expand("SplicingAnalysis/SPLICEq/polyA.{{mode}}/{IndID}.spliceq.tab", IndID=polyAsamples)
    wildcard_constraints:
        mode = 'IRjunctions|IER',
        Phenotype = 'polyA'
    params:
        "None"
        
use rule CollectSPLICEq_chRNA as CollectSPLICEq_M30 with:
    input: 
        expand("SplicingAnalysis/SPLICEq/MetabolicLabelled.30min.{{mode}}/{IndID}.spliceq.tab", IndID=metabolic30samples)
    wildcard_constraints:
        mode = 'IRjunctions|IER',
        Phenotype = 'MetabolicLabelled.30min'
    params:
        "None"


use rule CollectSPLICEq_chRNA as CollectSPLICEq_M60 with:
    input: 
        expand("SplicingAnalysis/SPLICEq/MetabolicLabelled.60min.{{mode}}/{IndID}.spliceq.tab", IndID=metabolic60samples)
    wildcard_constraints:
        mode = 'IRjunctions|IER',
        Phenotype = 'MetabolicLabelled.60min'
    params:
        "None"




