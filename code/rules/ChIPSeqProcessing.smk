rule MakeHisat2Index:
    input:
        "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2"
    log:
        "logs/MakeIndex.log"
    threads: 8
    shell:
        """
        hisat2-build {input} {input} -p {threads} &> {log}
        """

rule Align_From_SRA:
    """
    align reads from SRA accession number. If errors occur, check out if this
    solves the problem:
    https://standage.github.io/that-darn-cache-configuring-the-sra-toolkit.html
    """
    input:
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2"
    output:
        bam = "Alignments/Grubert_ChIPSeq/{Antibody}/{Cell_Line}.sorted.bam"
    params:
        SRA = lambda wildcards: Grubert_ChIP_seq_dict[(wildcards.Antibody, wildcards.Cell_Line)]
    log:
        "logs/ChIPSeq/Align_From_SRA/{Antibody}.{Cell_Line}.log"
    threads: 8
    shell:
        """
        (hisat2 -x ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa --sra-acc {params.SRA} --no-spliced-alignment --no-discordant -p {threads} | samtools sort -n | samtools fixmate -m - - |   samtools sort | samtools markdup - -  > {output.bam}) &> {log}
        """

rule ChIPSeq_BamToBigwig:
    input:
        bam = "Alignments/Grubert_ChIPSeq/{Antibody}/{Cell_Line}.sorted.bam",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
    output:
        bigwig = "Bigwigs/Grubert_ChIPSeq/{Antibody}.{Cell_Line}.bw"
    log:
        "logs/ChIPSeq/BamToBig/{Antibody}.{Cell_Line}.log"
    shell:
        """
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bigwig} -pc &> {log}
        """
