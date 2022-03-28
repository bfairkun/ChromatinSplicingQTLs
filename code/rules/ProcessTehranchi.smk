
rule DownloadTehranchiTableS2:
    output:
        "Tehranchi/TableS2.xlsx"
    shell:
        """
        wget -O {output} 'https://www.cell.com/cms/10.1016/j.cell.2016.03.041/attachment/f7554f28-b82d-4b06-99b0-cbcfa7698d4c/mmc2.xlsx'
        """

rule TehranchiTableS2_to_bed:
    input:
        "Tehranchi/TableS2.xlsx"
    output:
    conda:
        "envs/r_essentials.yml"
    shell:
        """

        """
