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

rule Download_hg19_to_hg38_chain:
    output:
        "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    shell:
        """
        wget -O {output} "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
        """

rule DownloadHg38Gencode_basic:
    """
    only transcripts flagged as basic
    """
    output:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf"
    shell:
        """
        wget -O- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.basic.annotation.gtf.gz | zcat > {output}
        """

rule Download1KG_GRCh38:
    output:
        vcf = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/{chrom}.vcf.gz.tbi",
    log:
        "logs/Download1KG_GRCh38/{chrom}.log"
    params:
        aspera_key = config['aspera_key']
    shell:
        """
        #Use aspera if aspera key file and aspera links parameters are defined (not empty strings)
        if [[ ! -z "{params.aspera_key}" ]]
        then
            ascp -i {params.aspera_key} -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz {output.vcf} &> {log}
            ascp -i {params.aspera_key} -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi {output.tbi} &>> {log}
        else
            curl -f -o {output.vcf} "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz" &> {log}
            curl -f -o {output.tbi} "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi" &>> {log}
        fi
        """

rule Gather1KGData:
    input:
        expand("Genotypes/1KG_GRCh38/{chrom}.vcf.gz", chrom=autosomes),

rule CopyFastqFromLocal:
    input:
        R1 = GetFastqLocalFuncs('R1_local'),
        R2 = GetFastqLocalFuncs('R2_local'),
    output:
        R1 = temp("Fastq/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz"),
        R2 = temp("Fastq/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz")
    log:
        "logs/MoveFastqFromLTS/{Phenotype}.{IndID}.{Rep}.log"
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log}
        cat {input.R2} > {output.R2} 2>> {log}
        """


rule DownloadFastqFromLink:
    """
    Download w aspera if possible, then try ftp link, then try fasterq-dump
    """
    input:
        aspera_key = config['aspera_key']
    output:
        fastq = temp("Fastq/{Phenotype}/{IndID}/{Rep}.{Read}.fastq.gz"),
    log:
        "logs/DownloadFastqFromAsperaLink/{Phenotype}.{IndID}.{Rep}.{Read}.log"
    wildcard_constraints:
        Read = "R1|R2"
    shadow: "shallow"
    params:
        ftp_link = GetDownloadLinkFuncs('ftp'),
        aspera_link  = GetDownloadLinkFuncs('aspera'),
    shell:
        """
        #Use aspera if aspera key file and aspera links parameters are defined (not empty strings)
        if [[ ! -z "{input.aspera_key}" && ! -z "{params.aspera_link}" ]]; then
            for link in {params.aspera_link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                ascp -v -QT -l 300m -P33001 -i {input.aspera_key} era-fasp@${{link}} $tmpfile &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        else
            for link in {params.ftp_link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                wget -O $tmpfile ${{link}} &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        fi
        """

use rule DownloadFastqFromLink as DownloadFastqFromLink_SE with:
    output:
        fastq = "FastqSE/{Phenotype}/{IndID}/{Rep}.{Read}.fastq.gz",
    log:
        "logs/DownloadFastqFromAsperaLink_SE/{Phenotype}.{IndID}.{Rep}.{Read}.log"
    wildcard_constraints:
        Read = "SE"


# TODO: make Download per fastq, and make inherited rule for SE

# def GetFastqR2(wildcards):
#     """
#     return empty list single ended
#     """
#     df_subset = Fastq_samples.loc[
#             (Fastq_samples['IndID'] == wildcards.IndID) &
#             (Fastq_samples['Phenotype'] == wildcards.Phenotype) &
#             (Fastq_samples['RepNumber'] == wildcards.Rep)]
#     if sum(df_subset['PairedEnd']) > 0:
#         return

SingleEnd_df = Fastq_samples.loc[ Fastq_samples['PairedEnd']==False, ['Phenotype', 'IndID', 'RepNumber'] ].drop_duplicates()
rule GatherAllFasterqDump_SE:
    input:
        expand("FastqSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz", zip, Phenotype=SingleEnd_df['Phenotype'], IndID=SingleEnd_df['IndID'], Rep=SingleEnd_df['RepNumber'])




rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "Fastq/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz"
    output:
        R1 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz",
        html = "FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.html",
        json = "FastqFastp/{Phenotype}/{IndID}/{Rep}.fastp.json"
    params:
        umi = GetFastpParamsUmi,
        I = "-I",
        O = "-O"
    resources:
        mem_mb = 8000
    log:
        "logs/fastp/{Phenotype}.{IndID}.{Rep}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} {params.I} {input.R2} -o {output.R1} {params.O} {output.R2} --html {output.html} --json {output.json} {params.umi} &> {log}
        """

use rule fastp as fastp_SE with:
    input:
        R1 = "FastqSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz",
        R2 = []
    output:
        R1 = "FastqFastpSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz",
        R2 = [],
        html = "FastqFastpSE/{Phenotype}/{IndID}/{Rep}.fastp.html",
        json = "FastqFastpSE/{Phenotype}/{IndID}/{Rep}.fastp.json"
    log:
        "logs/fastpSE/{Phenotype}.{IndID}.{Rep}.log"
    params:
        umi = GetFastpParamsUmi,
        I = "",
        O = ""
