
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
            ascp -i {params.aspera_key} -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz {output.vcf} &> {log}
            ascp -i {params.aspera_key} -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi {output.tbi} &>> {log}
        else
            curl -f -o {output.vcf} "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz" &> {log}
            curl -f -o {output.tbi} "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{wildcards.chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi" &>> {log}
        fi
        """

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
    input:
        aspera_key = config['aspera_key']
    output:
        R1 = temp("Fastq/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz"),
        R2 = temp("Fastq/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz")
    log:
        "logs/DownloadFastqFromAsperaLink/{Phenotype}.{IndID}.{Rep}.log"
    params:
        R1_ftp = GetDownloadLinkFuncs('R1_ftp'),
        R2_ftp = GetDownloadLinkFuncs('R2_ftp'),
        R1_aspera = GetDownloadLinkFuncs('R1_aspera'),
        R2_aspera = GetDownloadLinkFuncs('R2_aspera'),
    shell:
        """
        #Use aspera if aspera key file and aspera links parameters are defined (not empty strings)
        if [[ ! -z "{input.aspera_key}" && ! -z "{params.R1_aspera}" && ! -z "{params.R2_aspera}" ]]
        then
            ascp -v -QT -l 300m -P33001 -i {input.aspera_key} era-fasp@{params.R1_aspera} {output.R1} &> {log}
            ascp -v -QT -l 300m -P33001 -i {input.aspera_key} era-fasp@{params.R2_aspera} {output.R2} &> {log}
        else
            wget -O {output.R1} {params.R1_ftp} &> {log}
            wget -O {output.R2} {params.R2_ftp} &>> {log}
        fi
        """

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
        GetFastpParams
    resources:
        mem_mb = 8000
    log:
        "logs/fastp/{Phenotype}.{IndID}.{Rep}.log"
    conda:
        "envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} {params} &> {log}
        """

