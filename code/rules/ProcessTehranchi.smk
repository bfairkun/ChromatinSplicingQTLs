
rule DownloadTehranchiTableS2:
    output:
        temp("Tehranchi/TableS2.xlsx")
    shell:
        """
        wget -O {output} 'https://www.cell.com/cms/10.1016/j.cell.2016.03.041/attachment/f7554f28-b82d-4b06-99b0-cbcfa7698d4c/mmc2.xlsx'
        """

rule TehranchiTableS2_to_bed:
    input:
        "Tehranchi/TableS2.xlsx"
    output:
        temp("Tehranchi/TableS2.hg19.bed")
    conda:
        "../envs/r_slopes.yml"
    shell:
        """
        Rscript scripts/TehranchiTableS2_to_bed.R {input} {output}
        """

rule TehranchiBedToHg38:
    input:
        bed = "Tehranchi/TableS2.hg19.bed",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        "Tehranchi/TableS2.hg38.bed"
    conda:
        "../envs/crossmap.yml"
    log:
        "logs/TehranchiBedToHg38.log"
    shell:
        """
        CrossMap.py bed {input.chain} {input.bed} {output} &> {log}
        """

