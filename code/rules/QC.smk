
rule CountReadsPerSample:
    input:
        bam = AllBams,
        bai = AllBais
    output:
        "../output/QC/ReadCountsAndJunctionsPerSamples.tsv"
    log:
        "logs/CountReadsPerSample.log"
    shell:
        """
        # exec > {log} 2>&1
        # set -x
        for f in {input.bam}
        do
           printf "%s\\t%s\\n" $f $(samtools idxstats $f | awk -F'\\t' '$1~"^chr[1-9]" {{sum+=$3}} END {{print sum}}') >> {output}
        done
        """

# rule CountJunctionsPerSample:
#     input:

rule QualimapRnaseq:
    input:
        gtf="ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        bam=GetBamForBigwig,
        bai=GetBaiForBigwig
    output:
        "QC/QualimapRnaseq/{Phenotype}.{IndID}.{Rep}/rnaseq_qc_results.txt"
    log:
        "logs/QualimapRnaseq/{Phenotype}.{IndID}.{Rep}.log"
    conda:
        "../envs/qualimap.yml"
    params:
        libtype = GetQualimapLibtype
    resources:
        mem_mb = 16000
    shell:
        """
        unset DISPLAY
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} {params.libtype} --java-mem-size=12G -outdir QC/QualimapRnaseq/{wildcards.Phenotype}.{wildcards.IndID}.{wildcards.Rep}/ &> {log}
        """

rule GetCommonSnpsForMbv:
    input:
        "Genotypes/1KG_GRCh38/22.vcf.gz"
    output:
        vcf = "QC/mbv/chr22.vcf.gz",
        tbi = "QC/mbv/chr22.vcf.gz.tbi"
    log:
        "logs/GetCommonSnpsForMbv.log"
    shell:
        """
        bcftools view -O z -q 0.01:minor {input} > {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

rule PlotMultiPhenotypeHeatmap:
    input:
        Genes = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        Counts = expand("featureCounts/{PhenotypeRegion}/Counts.txt", PhenotypeRegion=["polyA.Expression", "chRNA.Expression", "AtTSS/H3K27AC", "AtTSS/H3K4ME3", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
    output:
        "QC/MultiphenotypeSpearmanHeatmap.pdf"
    log:
        "logs/PlotMultiPhenotypeHeatmap.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PlotSequentialPhenotypeHeatmap.R {input.Genes} {output} {input.Counts} &> {log}
        """


rule mbv:
    input:
        bam = GetBamForBigwig,
        vcf = "QC/mbv/chr22.vcf.gz",
        tbi = "QC/mbv/chr22.vcf.gz.tbi"
    output:
        text ="QC/mbv/data/{Phenotype}/{IndID}.{Rep}.txt",
        plot = "QC/mbv/plots/{Phenotype}/{IndID}.{Rep}.pdf"
    log:
        "logs/mbv/{Phenotype}/{IndID}.{Rep}.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        QTLtools_1.2_CentOS7.8_x86_64 mbv --vcf {input.vcf} --bam {input.bam} --out {output.text}  --reg chr22 &> {log}
        Rscript scripts/Plot_mbv.R {output.text} {output.plot} {wildcards.IndID} &>> {log}
        """


rule MultiQC:
    input:
        expand("QC/QualimapRnaseq/{Phenotype}.{IndID}.{Rep}/rnaseq_qc_results.txt",  zip, Phenotype=RNASeqSamples_df['Phenotype'], IndID=RNASeqSamples_df['IndID'], Rep=RNASeqSamples_df['RepNumber']),
        expand("Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Log.final.out", zip, Phenotype=RNASeqSamples_df['Phenotype'], IndID=RNASeqSamples_df['IndID'], Rep=RNASeqSamples_df['RepNumber']),
        expand("FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.json", zip, Phenotype=Fastq_samples['Phenotype'], IndID=Fastq_samples['IndID'], Rep=Fastq_samples['RepNumber']),
    log: "logs/Multiqc.log"
    output:
        "Multiqc/multiqc_report.html"
    shell:
        """
        multiqc -f -dd 3  -o Multiqc/ QC/QualimapRnaseq/ FastqFastp/ Alignments/STAR_Align/ &> {log}
        """
