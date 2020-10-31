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

juncfile_source = "regtools"

if juncfile_source == "STAR":
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

elif juncfile_source == "regtools":
    rule regtools_bam_to_junc:
        input:
            bam = "Alignments/GEUVADIS_RNAseq/{RNASeqSample}/Aligned.sortedByCoord.out.bam",
            fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
            gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
        output:
            junc = "Phenotypes/GEUVADIS_RNAseq/BySample/{RNASeqSample}.junc",
            regtools = "Phenotypes/GEUVADIS_RNAseq/BySample/{RNASeqSample}.regtools.junc"
        log:
            "logs/regtools_to_leafcutter_junc/{RNASeqSample}.log"
        params:
            regtools_junction_strand = 0
        shell:
            """
            (regtools junctions extract -s {params.regtools_junction_strand} {input.bam} | regtools junctions annotate -S - {input.fa} {input.gtf} | tee {output.regtools} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{ print $1, $2, $3, $4, $5, $6 }}' > {output.junc}) &> {log}
            """

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
    params:
        ExtraParams = "-s True"
    shell:
        """
        python2.7 scripts/leafcutter_cluster.py {params.ExtraParams} -j {input} -r Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/ &> {log}
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
    """
    phenotype table compatible with QTLtools (6 column bed format, bgzip
    compressed with tbi index file). groupid for each phenotype is the
    leafcutter clusterid, so permutation testing will be for phenotypes grouped
    by cluster
    """
    input:
        counts ="Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.autosomes.gz",
    output:
        phenotypes_perchrom = expand("Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr{chrom}", chrom=range(1, 23)),
        phenotypes = "Phenotypes/GEUVADIS_RNAseq/leafcutter/SplicingPhenotypesQQNormed.bed.gz",
        phenotypes_tbi = "Phenotypes/GEUVADIS_RNAseq/leafcutter/SplicingPhenotypesQQNormed.bed.gz.tbi",
        PCs = "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz.PCs"
    log:
        "logs/Phenotypes/GEUVADIS_RNAseq/leafcutter_prepare_phenotype_table.log"
    params:
        NumberPCs = 15
    conda:
        "../envs/py27.yaml"
    shell:
        """
        scripts/prepare_phenotype_table.py -p {params.NumberPCs} {input.counts}
        cat {output.phenotypes_perchrom} | awk -F'\\t' -v OFS='\\t' 'NR==1 {{print $1,$2,$3,"pid","gid","strand", $0}} NR>1 {{split($4,a,":"); split(a[4],b,"_"); print $1,$2,$3,$4[4],b[3],$0}}' | sed -r 's/(^(\S+\s+){{6}})(\S+\s+){{4}}/\\1/' | bedtools sort -i - | bgzip -c /dev/stdin > {output.phenotypes}
        tabix -p bed {output.phenotypes}
        """

rule sQTL_QTLtools_cis_permutation_pass:
    input:
        phenotypes = "Phenotypes/GEUVADIS_RNAseq/leafcutter/SplicingPhenotypesQQNormed.bed.gz",
        phenotypes_tbi = "Phenotypes/GEUVADIS_RNAseq/leafcutter/SplicingPhenotypesQQNormed.bed.gz.tbi",
        genotypes = "Genotypes/GEUVADIS_Lappalainnen.vcf",
        genotypes_tbi = "Genotypes/GEUVADIS_Lappalainnen.vcf",

        PCs = "Phenotypes/GEUVADIS_RNAseq/leafcutter/clustering/leafcutter_perind.counts.gz.PCs"
    params:
        CisWindow = 100000
    output:
        chunk = "QTLs/sQTLs/permutation_pass_chunks/{chunk_num}.txt"
    shell:
        """
        QTLtools_1.2_CentOS7.8_x86_64 cis --permute 10000 --vcf {input.genotypes} --bed {phenotypes} --cov {input.PCs} --out {output.chunk} --grp-best
        """

# rule sQTL_QTLtools_cis_nominal_pass:
#     """
#     Do a nominal cis pass to get all phenotype:SNP pairs for phenotypes that are significant in permutation pass. This might be necessary for colocalization methods.
#     """
