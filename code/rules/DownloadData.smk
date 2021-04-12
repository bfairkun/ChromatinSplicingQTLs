rule Download_GEUVADIS_RNAseq:
    output:
        "Fastq/GEUVADIS_RNAseq/{ScanName}"
    params:
        lambda wildcards: GEUVADIS_sample_links_dict[wildcards.ScanName]
    shell:
        """
        wget -O {output} {params}
        """

rule DownloadHg38Ref:
    output:
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        primary_gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        chr_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    shell:
        """
        wget -O- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz | zcat > {output.fa}
        wget -O- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz | zcat > {output.primary_gtf}
        wget -O- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz | zcat > {output.chr_gtf}
        """

rule indexHg38Ref:
    input:
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
    output:
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
    shell:
        """
        samtools faidx {input.fa}
        """

rule Download_GEUVADIS_Genotypes_official:
    """
    The official vcf from 1000Genomes project downloaded, filtered for GEUVADIS
    samples, and maybe some other things.
    """
    input:
        GEU_samples = "../data/E-GEUV-1.sdrf.txt",
    output:
        vcf = "Genotypes/GEUVADIS_1000G/{Chrom}.vcf.gz",
        tbi = "Genotypes/GEUVADIS_1000G/{Chrom}.vcf.gz.tbi",
    shell:
        """
        wget -O - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.Chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools view --force-samples -S <(awk -F'\\t' 'NR>1 {{print $1}}' {input.GEU_samples} | sort | uniq) -O z -m 2 -M 2 - > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule Download_Grubert_ChIP_seq:
    """
    Rule is actually obselete, since Hisat2 aligner can accept SRA accesion as
    input, so no need to download raw fastq
    """
    input:
        Grubert_samples = "../data/PRJNA268086_SraRunTable.GrubertEtAl.csv"
    output:
        "Fastq/Grubert_{antibody}/{SRR}_1.fastq.gz",
        "Fastq/Grubert_{antibody}/{SRR}_2.fastq.gz",
    log:
        "logs/Download_Grubert_ChIP_seq/{antibody}.{SRR}.log"
    shell:
        """
        fasterq-dump {wildcards.SRR} --split-files -O Fastq/Grubert_{wildcards.antibody}  &> {log}
        gzip Fastq/Grubert_{wildcards.antibody}/{wildcards.SRR}_1.fastq
        gzip Fastq/Grubert_{wildcards.antibody}/{wildcards.SRR}_2.fastq
        """

# rule Download
