rule idxstats:
    input:
        bam = GetBamForBigwig,
        bai=GetBaiForBigwig
    output:
        "QC/idxstats/{Phenotype}/{IndID}.{Rep}.txt"
    log:
        "logs/idxstats/{Phenotype}/{IndID}.{Rep}.log"
    shell:
        """
        samtools idxstats {input.bam} > {output} 2> {log}
        """

rule CountAutosomalReadsPerSample:
    input:
        expand("QC/idxstats/{Phenotype}/{IndID}.{Rep}.txt", zip, Phenotype=Fastq_samples['Phenotype'], IndID=Fastq_samples['IndID'], Rep=Fastq_samples['RepNumber'])
    output:
        "../output/QC/AutosomeCountsPerSamples.tsv"
    shell:
        """
        awk -F'\\t' '{{ print $1, $2, $3, FILENAME }}' {input} > {output}
        """

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
