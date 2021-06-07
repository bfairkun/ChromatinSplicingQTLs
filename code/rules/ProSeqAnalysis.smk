rule CollectProSeqData:
    input:
        "ProCapAnalysis/SupplementDat2.csv",
        expand("ProCapAnalysis/RefSeqClassification_{tTRE_type}.bed", tTRE_type=["distal_enhancer", "promoter", "proximal_enhancer"])

rule DownloadKristjansdottirSupplementData2:
    output:
        "ProCapAnalysis/SupplementDat2.csv"
    shell:
        """
        wget -O {output} https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19829-z/MediaObjects/41467_2020_19829_MOESM6_ESM.csv
        """

rule GetRefSeqtTRE_beds:
    input:
        "ProCapAnalysis/SupplementDat2.csv"
    output:
        "ProCapAnalysis/RefSeqClassification_{tTRE_type}.bed"
    shell:
        """
        awk -F, -v OFS='\\t' 'NR>1 && $9==0 && $5=="{wildcards.tTRE_type}" {{print $1, $2, $2+1, ".", ".", $4}}' {input} > {output}
        """

rule MakeProCapPhenotypeTable:
    input:
        "ProCapAnalysis/SupplementDat2.csv"
    output:
        "QTLs/QTLTools/ProCap/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/MakeProCapPhenotypeTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTable_ProCap.R {input} {output} &> {log}
        """

