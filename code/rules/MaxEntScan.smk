rule IntronAnnotationToBED_MES:
    input:
        "/project2/yangili1/bjf79/ChromatinSplicingQTLs/data/IntronAnnotationsFromYang.tsv.gz"
    output:
        "SplicingAnalysis/MaxEntScan/Annotation/IntronAnnotationsFromYang.bed.gz"
    log:
        'logs/MaxEntScan/Annotation.log'
    resources:
        mem_mb = 8000
    shell:
        """
        (zcat {input} | awk '{{
        print $1, $2, $3, $1":"$2":"$3":"$4, $7, $4, $5, $6, $8
        }}' FS='\\t' OFS='\\t' - | tail -n+2 - | bedtools sort -i - | gzip - > {output}) &> {log}
        """
        
rule GetIntronFivePrime_MES:
    input:
        "SplicingAnalysis/MaxEntScan/Annotation/IntronAnnotationsFromYang.bed.gz"
    output:
        "SplicingAnalysis/MaxEntScan/Annotation/FivePrime.bed.gz"
    log:
        "logs/MaxEntScan/GetFivePrimeFromAnnotation.log"
    resources:
        mem_mb = 8000
    shell:
        """
        (zcat {input} | awk '{{
        if ($6=="+")
            print $1, $2-3, $2+6, $4, $8, $6
        else
            print $1, $3-7, $3+2, $4, $8, $6
        }}' FS='\\t' OFS='\\t' - | gzip - > {output}) &> {log}
        """

rule GetIntronThreePrime_MES:
    input:
        "SplicingAnalysis/MaxEntScan/Annotation/IntronAnnotationsFromYang.bed.gz"
    output:
        "SplicingAnalysis/MaxEntScan/Annotation/ThreePrime.bed.gz"
    log:
        "logs/MaxEntScan/GetThreePrimeFromAnnotation.log"
    resources:
        mem_mb = 8000
    shell:
        """
        (zcat {input} | awk '{{
        if ($6=="+")
            print $1, $3-21, $3+2, $4, $8, $6
        else
            print $1, $2-3, $2+20, $4, $8, $6
        }}' FS='\\t' OFS='\\t' - | gzip - > {output}) &> {log}
        """

    
rule GetSpliceSiteFasta_MES:
    input:
        fasta = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        bed = "SplicingAnalysis/MaxEntScan/Annotation/{Junction}.bed.gz"
    output:
        "SplicingAnalysis/MaxEntScan/Annotation/{Junction}.tab"
    log:
        "logs/MaxEntScan/GetSpliceSiteFasta.{Junction}.tab.log"
    wildcard_constraints:
        Junction = "FivePrime|ThreePrime",
    resources:
        mem_mb = 12000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        bedtools getfasta -name -s -tab -fi {input.fasta} -bed {input.bed} -fo {output} &> {log}
        """

rule ScoreMaxEntScan_MES:
    input:
        "SplicingAnalysis/MaxEntScan/Annotation/{Junction}.tab" 
    output:
        "SplicingAnalysis/MaxEntScan/Annotation/{Junction}.score.tab.gz",
    log:
        "logs/MaxEntScan/{Junction}.MaxEntScan.log"
    wildcard_constraints:
        Junction = "FivePrime|ThreePrime",
    conda:
        "../envs/maxentpy.yml"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/ScoreMaxEntScan.py {input} {output} &> {log}
        """

rule CollectMaxtEntScan:
    input:
        FivePrimeMES = "SplicingAnalysis/MaxEntScan/Annotation/FivePrime.score.tab.gz",
        ThreePrimeMES = "SplicingAnalysis/MaxEntScan/Annotation/ThreePrime.score.tab.gz",