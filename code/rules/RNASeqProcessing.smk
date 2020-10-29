# rule CutAdapters:
#     input:
#     output:

rule STAR_make_index:
    """
    did not work on bigmem2. Never figured out why (the log file didn't
    indicate anything). Ran on login node with success.
    """
    input:
        fasta = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
    log:
        "logs/STAR_make_index.log"
    params:
        genomeDir = "ReferenceGenome/STARIndex/"
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """


rule STAR_alignment:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        vcf = "Genotypes/GEUVADIS_1000G/All.vcf.gz",
        R1_fastq = lambda wildcards: "Fastq/GEUVADIS_RNAseq/" + GEUVADIS_line_fastq_dict[wildcards.RNASeqSample] + "1.fastq.gz",
        R2_fastq = lambda wildcards: "Fastq/GEUVADIS_RNAseq/" + GEUVADIS_line_fastq_dict[wildcards.RNASeqSample] + "2.fastq.gz"
    log:
        "logs/GEUVADIS_RNAseq/STAR/{RNASeqSample}.log"
    threads: 12
    params:
        AdapterClip = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        ReadNum = "-1"
    output:
        "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Aligned.sortedByCoord.out.bam",
        "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Log.final.out"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir ReferenceGenome/STARIndex/ --readFilesIn {input.R1_fastq} {input.R2_fastq} --outSAMtype BAM SortedByCoordinate --clip3pAdapterSeq {params.AdapterClip} --readMapNumber {params.ReadNum} --outSAMmultNmax 1 --alignEndsType EndToEnd --twopassMode Basic --readFilesCommand zcat --outFileNamePrefix Alignments/GEUVADIS_RNAseq/{wildcards.RNASeqSample}/ &> {log}
        # STAR --runThreadN {threads} --genomeDir ReferenceGenome/STARIndex/ --varVCFfile <(zcat {input.vcf}) --readFilesIn {input.R1_fastq} {input.R2_fastq} --waspOutputMode SAMtag --outSAMattributes vA vG vW  --outSAMtype BAM SortedByCoordinate --clip3pAdapterSeq {params.AdapterClip} --readMapNumber {params.ReadNum} --outSAMmultNmax 1 --alignEndsType EndToEnd --twopassMode Basic --readFilesCommand zcat --outFileNamePrefix Alignments/GEUVADIS_RNAseq/{wildcards.RNASeqSample}/ &> {log}
        """

rule index_RNA_seq_bams:
    input:
        bam = "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Aligned.sortedByCoord.out.bam"
    output:
        bai= "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule RNASeq_BamToBigwig:
    input:
        bam = "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Aligned.sortedByCoord.out.bam",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
    output:
        bigwig = "Bigwigs/GEUVADIS_RNAseq/{RNASeqSample}.bw"
    log:
        "logs/GEUVADIS_RNAseq/BamToBig/{RNASeqSample}.log"
    shell:
        """
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bigwig} -split &> {log}
        """

rule STAR_to_leafcutter_junc:
    input:
        "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/SJ.out.tab"
    output:
        "Phenotypes/GEUVADIS_RNAseq/BySample/{RNASeqSample}.junc"
    log:
        "logs/STAR_to_leafcutter_junc/{RNASeqSample}.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4==1 && $1!="MT" {{ print $1,$2,$3,".",$7,"+" }} $4==2&& $1!="MT" {{ print $1,$2,$3,".",$7,"-" }}' {input} > {output}
        """

#TODO: regtools to make juncfiles.

rule make_leafcutter_juncfile:
    input:
        expand ("Phenotypes/GEUVADIS_RNAseq/BySample/{RNASeqSample}.junc", RNASeqSample=GEUVADIS_line_fastq_dict.keys()),
    output:
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/juncfilelist.txt"
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
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/juncfilelist.txt",
    output:
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz",
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind_numers.counts.gz"
    log:
        "logs/Phenotypes/GEUVADIS_RNAseq/leafcutter_cluster.log"
    shell:
        """
        leafcutter_cluster.py -s -j {input} -r Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/ &> {log}
        """

rule MakeLeafcutterBlacklistChromsFile:
    input:
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz",
    output:
        "ReferenceGenome/Fasta/GRCh38.primary_assembly.NonAutosomes.list",
    shell:
        """
        zcat {input} | awk -F: 'NR>1 {{print $1}}' | sort | uniq | grep -v -P "^chr\d+" > {output}
        """

rule FilterLeafcutterClustersForAutosomes:
    input:
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz",
    output:
        "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.autosomes.gz",
    shell:
        """
        cat <(zcat {input} | head -1) <(zcat {input} | awk 'NR>1' | grep -P "chr\d+" | sed 's/^chr//') | gzip - > {output}
        """

rule leafcutter_prepare_phenotype_table:
    input:
        counts ="Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.autosomes.gz",
    output:
        phenotypes = "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm",
        PCs = "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz.PCs"
    log:
        "logs/Phenotypes/GEUVADIS_RNAseq/leafcutter_prepare_phenotype_table.log"
    conda:
        "../envs/py27.yaml"
    shell:
        """
        scripts/prepare_phenotype_table.py -p 15 {input.counts}
        """

# rule leafcutter_prepare_cluster_groups:
    
