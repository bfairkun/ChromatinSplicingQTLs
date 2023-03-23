rule GatherSMFastq:
    input:
        # expand("SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam.bai", Sample=SM_samples['SampleName'].unique()),
        "SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        "SmallMolecule/featureCounts/Counts.txt",
        "SmallMolecule/FitModels/polyA_genes.tsv.gz"

use rule CopyFastqFromLocal as SM_CopyFastqFromLocal with:
    input:
        R1 = lambda wildcards: SM_samples.loc[SM_samples['SampleName'] == wildcards.Sample]['R1'],
        R2 = lambda wildcards: SM_samples.loc[SM_samples['SampleName'] == wildcards.Sample]['R2']
    output:
        R1 = temp("SmallMolecule/Fastq/{Sample}.R1.fastq.gz"),
        R2 = temp("SmallMolecule/Fastq/{Sample}.R2.fastq.gz")
    log:
        "logs/SM_MoveFastqFromLTS/{Sample}.log"

use rule fastp as SM_fastp with:
    input:
        R1 = "SmallMolecule/Fastq/{Sample}.R1.fastq.gz",
        R2 = "SmallMolecule/Fastq/{Sample}.R2.fastq.gz"
    output:
        R1 = "SmallMolecule/FastqFastp/{Sample}.R1.fastq.gz",
        R2 = "SmallMolecule/FastqFastp/{Sample}.R2.fastq.gz",
        html = "SmallMolecule/FastqFastp/{Sample}.fastp.html",
        json = "SmallMolecule/FastqFastp/{Sample}.fastp.json"
    params:
        I = "-I",
        O = "-O",
        # extra = "--reads_to_process 10000"
        extra = ""
    log:
        "logs/SM_fastp/{Sample}.log"

rule STAR_MultiSample2PassMapping_Pass1:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "SmallMolecule/FastqFastp/{Sample}.R1.fastq.gz",
        R2 = "SmallMolecule/FastqFastp/{Sample}.R2.fastq.gz",
        SJ = []
    output:
        bam = "SmallMolecule/AlignmentsPass1/{Sample}/Aligned.sortedByCoord.out.bam",
        align_log = "SmallMolecule/AlignmentsPass1/{Sample}/Log.final.out",
        SJ = "SmallMolecule/AlignmentsPass1/{Sample}/SJ.out.tab",
    threads: 8
    log: "logs/SM_Pass1/STAR_Align/{Sample}.log"
    params:
        ParentDirPrefix = lambda wildcards: f"SmallMolecule/AlignmentsPass1/{wildcards.Sample}/" ,
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        sjdb = ""
    resources:
        cpus_per_node = 9,
        mem = 58000,
    shell:
        """
        # module load STAR/2.7.7a
        STAR --readMapNumber {params.readMapNumber}  --outFileNamePrefix {params.ParentDirPrefix} --genomeDir ReferenceGenome/STARIndex/ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1  --limitBAMsortRAM 16000000000 {params.ENCODE_params} --outSAMstrandField intronMotif {params.sjdb} {input.SJ}  &> {log}
        """

use rule STAR_MultiSample2PassMapping_Pass1 as STAR_MultiSample2PassMapping_Pass2 with:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "SmallMolecule/FastqFastp/{Sample}.R1.fastq.gz",
        R2 = "SmallMolecule/FastqFastp/{Sample}.R2.fastq.gz",
        SJ = expand("SmallMolecule/AlignmentsPass1/{Sample}/SJ.out.tab", Sample=SM_samples['SampleName'].unique())
    output:
        bam = "SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam",
        align_log = "SmallMolecule/AlignmentsPass2/{Sample}/Log.final.out",
        SJ = "SmallMolecule/AlignmentsPass2/{Sample}/SJ.out.tab",
    params:
        ParentDirPrefix = lambda wildcards: f"SmallMolecule/AlignmentsPass2/{wildcards.Sample}/" ,
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        sjdb = "--limitSjdbInsertNsj 2494398 --sjdbFileChrStartEnd"
    log:
        "logs/SM_Pass2/STAR_Align/{Sample}.log"

rule SM_indexBam:
    input:
        bam = "SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam",
    output:
        bai = "SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """

rule SM_ExtractJuncs:
    input:
        bam = "SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam",
        bai = "SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        junc = temp(expand("SmallMolecule/leafcutter/juncfiles/chr{chrom}/{{Sample}}.junc", chrom=autosomes)),
        junc_autosomes = "SmallMolecule/leafcutter/juncfiles/autosomes/{Sample}.junc",
    params:
        # strand = GetLibStrandForRegtools
        strand = 0
    conda:
        "../envs/regtools.yml"
    log:
        "logs/SM_ExtractJuncs/{Sample}.log"
    shell:
        """
        for chrom in {autosomes}
        do
            (regtools junctions extract -m 20 -s {params.strand} -r chr${{chrom}} {input.bam} > SmallMolecule/leafcutter/juncfiles/chr${{chrom}}/{wildcards.Sample}.junc ) &> {log}
        done
        cat {output.junc} > {output.junc_autosomes}
        """

use rule make_leafcutter_juncfile as SM_make_leafcutter_juncfile with:
    input:
        expand("SmallMolecule/leafcutter/juncfiles/autosomes/{Sample}.junc", Sample=SM_samples['SampleName'].unique())
    output:
        "SmallMolecule/leafcutter/juncfilelist.autosomes.txt"

use rule leafcutter_cluster as SM_leafcutter_cluster with:
    input:
        juncs = expand("SmallMolecule/leafcutter/juncfiles/autosomes/{Sample}.junc", Sample=SM_samples['SampleName'].unique()),
        juncfile_list = "SmallMolecule/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    params:
        rundir = "-r SmallMolecule/leafcutter/clustering/autosomes/"
    log:
        "logs/SM_leafcutter_cluster/autosomes.log"

rule SM_collapse_juncfiles:
    """
    In preparation for regtools annotate using a juncfile with just one line for each unique junction among all junctions discovered across all samples
    """
    input:
        juncs = expand("SmallMolecule/leafcutter/juncfiles/autosomes/{Sample}.junc", Sample=SM_samples['SampleName'].unique()),
    output:
        "SmallMolecule/FullSpliceSiteAnnotations/ALL_SAMPLES.junc"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/Collapse_Juncfiles.R {output} {input}
        """

rule SM_annotate_juncfiles_allsamples:
    input:
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        basic_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        juncs = "SmallMolecule/FullSpliceSiteAnnotations/ALL_SAMPLES.junc"
    output:
        basic = "SmallMolecule/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz",
        comprehensive = "SmallMolecule/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.comprehensive.bed.gz"
    log:
        "logs/SM_annotate_juncfiles_allsamples.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.juncs} {input.fa} {input.basic_gtf} | gzip - > {output.basic} ) &> {log}
        (regtools junctions annotate {input.juncs} {input.fa} {input.Comprehensive_gtf} | gzip - > {output.comprehensive} ) &>> log
        """

rule Get5ssSeqs:
    """
    Filtered out entries with N in sequence
    """
    input:
        basic = "SmallMolecule/FullSpliceSiteAnnotations/JuncfilesMerged.annotated.basic.bed.gz",
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
    output:
        temp("SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab")
    shell:
        """
        zcat {input.basic} | awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $1, $2, $3, $1"_"$2"_"$3"_"$6, ".", $6}}' | sort -u | awk -v OFS='\\t' -F'\\t'  '$6=="+" {{$2=$2-4; $3=$2+11; print $0}} $6=="-" {{$3=$3+3; $2=$3-11; print $0}}' | bedtools getfasta -tab -bed - -s -name -fi {input.fa} | grep -v 'N' > {output}
        """

rule Score5ss:
    input:
        "SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab" 
    output:
        Tab = "SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
        WeblogoImg = "../docs/assets/5ssPWM.png"
    conda:
        "../envs/biopython.yml"
    shell:
        """
        python scripts/ScorePWM.py {input} {output.Tab} {output.WeblogoImg}
        """

rule SM_featurecounts:
    input:
        bam = expand("SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam", Sample=SM_samples['SampleName'].unique()),
        bai = expand("SmallMolecule/AlignmentsPass2/{Sample}/Aligned.sortedByCoord.out.bam.bai", Sample=SM_samples['SampleName'].unique()),
        basic_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
    output:
        "SmallMolecule/featureCounts/Counts.txt"
    threads:
        8
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts.log"
    params:
        extra = "-p"
    shell:
        """
        featureCounts {params.extra} -T {threads} --ignoreDup --primary -a {input.basic_gtf} -o {output} {input.bam} &> {log}
        """

rule FitDoseResponseLogLogisticModel:
    input:
        GeneCounts = "SmallMolecule/featureCounts/Counts.txt",
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        SplicingCounts = "SmallMolecule/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz",
        SpliceDonorSeqs = "SmallMolecule/leafcutter/JuncfilesMerged.annotated.basic.bed.5ss.tab.gz",
    output:
        GeneModelParams = "SmallMolecule/FitModels/polyA_genes.tsv.gz",
        IntronModelParams = "SmallMolecule/FitModels/polyA_GAGTIntrons.tsv.gz",
    log:
        "logs/FitDoseResponseLogLogisticModel.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/FitSmallMoleculeModels.R &> {log}
        """


rule Get_Fit_GAGTInts_bed:
    input:
        IntronModelParams = "SmallMolecule/FitModels/polyA_GAGTIntrons.tsv.gz",
    output:
        bed = "SmallMolecule/FitModels/polyA_GAGTIntrons.bed.gz",
        tbi = "SmallMolecule/FitModels/polyA_GAGTIntrons.bed.gz.tbi",
    shell:
        """
        zcat {input.IntronModelParams} | awk -F'\\t' 'NR>1 {{print $1}}' | sort | uniq |  awk -F'\\t' -v OFS='\\t' '{{split($1, a, ":"); print a[1], a[2], a[3], $1, ".", a[4]}}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule IndexCassetteExons:
    """
    input bed file from ../analysis/20230314_ProcessSM_ForInterpretableSplicingEffectSizes.Rmd
    """
    input:
        bed = "../output/SmallMoleculeGAGT_CassetteExonclusters.bed",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        FlankingBases = "SmallMolecule/CassetteExons/FlankingBases.tsv"
    shell:
        """
        paste -d'\\t' <(grep 'junc.skipping' {input.bed} | awk -F'\\t' -v OFS='\\t' '{{$3=$3-1; print $0}}' | bedtools flank -s -l 1 -r 0 -i - -g {input.fai}) <(grep 'junc.skipping' {input.bed} | awk -F'\\t' -v OFS='\\t' '{{$3=$3-1; print $0}}' | bedtools flank -s -l 0 -r 1 -i - -g {input.fai}) | awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "SkipJuncName", "Chrom", "strand", "UpstreamFlankBaseStart", "UpstreamFlankBaseEnd", "DownstreamFlankBaseStart", "DownstreamFlankBaseEnd"}} $4==$13 {{print $4, $1, $6, $2, $3, $11, $12}}' > {output}
        """


# rule featureCountInducedCassetteExons:
#     input:
#     output:
#         saf = "SmallMolecule/CassetteExons/CassetteExons.saf",
#         counts = "SmallMolecule/CassetteExons/Counts.txt"
