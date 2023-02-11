ReadLengthMapExperiment_samples = pd.read_csv("config/RNASeqReadLengthMapExperiment.samples.tsv", sep='\t', comment='#')

# rule sjdb_ForMultiSampleTwoPass:
#     """
#     because of temporary files, something kind of weird is happening with the DAG and unnecessarily realigning some things... to avoid this, just comment this out after the file is made with all juncs across all RNAseq experiments, but before calling rules with ReadLengthMapExperiment* in rule name.
#     """
#     input:
#         "SplicingAnalysis/regtools_annotate_combined/basic.bed.gz"
#     output:
#         "SplicingAnalysis/STAR_Multisample_sjdb.tab"
#     shell:
#         """
#         zcat {input} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1, $2+1, $3}}' | bedtools sort -i - | uniq > {output}
#         """

use rule fastp as ReadLengthMapExperiment_fastp with:
    input:
        R1 = "Fastq/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "Fastq/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz"
    output:
        R1 = temp("ReadLengthMapExperiment/FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz"),
        R2 = temp("ReadLengthMapExperiment/FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz"),
        html = "ReadLengthMapExperiment/FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.html",
        json = "ReadLengthMapExperiment/FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.json"
    wildcard_constraints:
        Phenotype = "Expression.Splicing|chRNA.Expression.Splicing"
    params:
        I = "-I",
        O = "-O",
        extra = "-b 50 -B 50"
    log:
        "logs/ReadLengthMapExperiment_fastp/{Phenotype}.{IndID}.{Rep}.log"

use rule fastp as ReadLengthMapExperiment_fastp_SE with:
    input:
        R1 = "FastqSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz",
        R2 = []
    wildcard_constraints:
        Phenotype = "MetabolicLabelled.30min|MetabolicLabelled.60min|ProCap"
    output:
        R1 = temp("ReadLengthMapExperiment/FastqFastpSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz"),
        R2 = [],
        html = "ReadLengthMapExperiment/FastqFastpSE/{Phenotype}/{IndID}/{Rep}.fastp.html",
        json = "ReadLengthMapExperiment/FastqFastpSE/{Phenotype}/{IndID}/{Rep}.fastp.json"
    log:
        "logs/ReadLengthMapExperiment_fastp_SE/{Phenotype}.{IndID}.{Rep}.log"
    params:
        I = "",
        O = "",
        extra = "-b 50"

use rule STAR_Align_WASP as ReadLengthMapExperiment_STAR_Align_WASP with:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "ReadLengthMapExperiment/FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = [],
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf",
        sjdb = "SplicingAnalysis/STAR_Multisample_sjdb.tab"
    output:
        bam = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam",
        align_log = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Log.final.out"
    log:
        "logs/ReadLengthMapExperiment_STAR_Align_WASP/{Phenotype}/{IndID}.{Rep}.log"
    wildcard_constraints:
        Phenotype = "Expression.Splicing|chRNA.Expression.Splicing"
    params:
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbFileChrStartEnd SplicingAnalysis/STAR_Multisample_sjdb.tab --limitSjdbInsertNsj 3845004",
        WASP_params = "--waspOutputMode SAMtag --outSAMattributes NH HI AS nM XS vW --varVCFfile",
        JunctionScore = GetSTARJunctionScoreParams,
        PrefixPrefix =  "ReadLengthMapExperiment/"

use rule ReadLengthMapExperiment_STAR_Align_WASP as ReadLengthMapExperiment_STAR_Align_WASP_SE with:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "ReadLengthMapExperiment/FastqFastpSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz",
        R2 = [],
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf",
        sjdb = "SplicingAnalysis/STAR_Multisample_sjdb.tab"
    wildcard_constraints:
        Phenotype = "MetabolicLabelled.30min|MetabolicLabelled.60min|ProCap"

use rule FilterBAM_WaspTags as ReadLengthMapExperiment_FilterBAM_WaspTags with:
    input:
        bam = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam"
    output:
        bam = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam",
        bai = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam.bai",
    log:
        "logs/ReadLengthMapExperiment_FilterBAM_WaspTags/{Phenotype}/{IndID}.{Rep}.log"

rule ExtractAndAnnotateJuncs:
    input:
        bam = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam",
        bai = "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam.bai",
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        "ReadLengthMapExperimentSpliceCounts/juncfiles/{MinAnchor}/{Phenotype}/{IndID}.{Rep}.junc.gz"
    log:
        "logs/ExtractJuncs/{MinAnchor}.{Phenotype}.{IndID}.{Rep}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions extract -m 20 -a {wildcards.MinAnchor} -s 0 {input.bam} | regtools junctions annotate /dev/stdin {input.fa} {input.gtf} | gzip - > {output} ) &> {log}
        """

rule ReadLengthMapExperiment_featureCounts:
    input:
        bam = expand(
            "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        bai = expand(
            "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam.bai",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        annotations = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf"
    output:
        "ReadLengthMapExperimentResults/featureCounts/Counts.txt"
    threads:
        8
    resources:
        mem = 12000,
        cpus_per_node = 9,
    params:
        extraParams = ""
    log:
        "logs/ReadLengthMapExperiment_featureCounts.log"
    shell:
        """
        featureCounts {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

rule ReadLengthMapExperiment_ProcessJuncFiles:
    input:
        juncs = expand(
            "ReadLengthMapExperimentSpliceCounts/juncfiles/8/{Phenotype}/{IndID}.{Rep}.junc.gz",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        annotations ="/project2/yangili1/yangili/chRNA/annotation_leaf_JAN28.txt.gz" 
    output:
        SummedByDataset = "ReadLengthMapExperimentResults/JuncCountsTidy/SummarisedByJuncAndDataset.txt.gz",
        longtable = "ReadLengthMapExperimentResults/JuncCountsTidy/LongTable.txt.gz",
        SummedByAnnotation = "ReadLengthMapExperimentResults/JuncCountsTidy/SummarisedBySampleAndAnnotation.txt.gz"
    conda:
        "../envs/r_2.yaml"
    log:
        "logs/ReadLengthMapExperiment_ProcessJuncFiles.log"
    shell:
        """
        Rscript scripts/ReadLengthMapExperiment.ProcessJuncFiles.R {output.SummedByDataset} {output.longtable} {output.SummedByAnnotation} &> {log}
        """

rule ReadLengthMapExperiment_CollectBams:
    input:
        expand(
            "ReadLengthMapExperiment/{Phenotype}/{IndID}/{Rep}/Filtered.bam",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        expand(
            "ReadLengthMapExperimentSpliceCounts/juncfiles/1/{Phenotype}/{IndID}.{Rep}.junc.gz",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        expand(
            "ReadLengthMapExperimentSpliceCounts/juncfiles/8/{Phenotype}/{IndID}.{Rep}.junc.gz",
            zip,
            Phenotype=ReadLengthMapExperiment_samples["Phenotype"],
            IndID=ReadLengthMapExperiment_samples["IndID"],
            Rep=ReadLengthMapExperiment_samples["RepNumber"],
        ),
        "ReadLengthMapExperimentResults/featureCounts/Counts.txt",
        "ReadLengthMapExperimentResults/JuncCountsTidy/SummarisedByJuncAndDataset.txt.gz" 
