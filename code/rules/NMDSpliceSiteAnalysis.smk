rule LeafCutterCountsPerJunction:
    input:
        'SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.{Phenotype}.gz'
    output:
        'SplicingAnalysis/NMDJunctions/tables/leafcutter_perind.counts_only.{Phenotype}.bed.gz'
    log:
        'logs/NMDSpliceSite/LeafCutterCountsPerJunction.{Phenotype}.log'
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/LeafCutterCountsPerJunction.py --input {input} --output {output} &> {log}
        """
        
rule IntronAnnotationToBED:
    input:
        "/project2/yangili1/yangili/chRNA/annotation_leaf_JAN28.txt.gz"
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed.gz"
    log:
        'logs/NMDSpliceSite/Annotation.log'
    resources:
        mem_mb = 8000
    shell:
        """
        (awk '{{
        print $1, $2, $3, $1":"$2":"$3":"$4, $7, $4, $5, $6, $8
        }}' FS='\\t' OFS='\\t' {input}  | bedtools sort -i - | gzip - > {output}) &> {log}
        """
        
rule GetIntronFivePrime:
    input:
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed.gz"
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/FivePrime.bed.gz"
    log:
        "logs/NMDSpliceSite/GetFivePrimeFromAnnotation.log"
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
  
def GetExtensionForAnnot(wildcards):
    if wildcards.Annot == '.ForMaxEntScan':
        return "14"
    else:
        return "0"

rule GetIntronThreePrime:
    input:
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed.gz"
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/ThreePrime{Annot}.bed.gz"
    log:
        "logs/NMDSpliceSite/GetThreePrimeFromAnnotation{Annot}.log"
    resources:
        mem_mb = 8000
    wildcard_constraints:
        Annot = "|.ForMaxEntScan"
    params:
        GetExtensionForAnnot
    shell:
        """
        (zcat {input} | awk '{{
        if ($6=="+")
            print $1, $3-7-{params}, $3+2, $4, $8, $6
        else
            print $1, $2-3, $2+6+{params}, $4, $8, $6
        }}' FS='\\t' OFS='\\t' - | gzip - > {output}) &> {log}
        """
    
def GetFastaOtherParams(wildcards):
    if wildcards.ext == 'tab':
        return "-tab"
    else:
        return ""
    
rule GetSpliceSiteFasta:
    input:
        fasta = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        bed = "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.bed.gz"
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.{ext}"
    log:
        "logs/NMDSpliceSite/GetSpliceSiteFasta.{Junction}.{ext}.log"
    wildcard_constraints:
        Junction = "FivePrime|ThreePrime|ThreePrime.ForMaxEntScan",
        ext = "tab|fa"
    params:
        GetFastaOtherParams
    resources:
        mem_mb = 12000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        bedtools getfasta -name -s {params} -fi {input.fasta} -bed {input.bed} -fo {output} &> {log}
        """

rule ScorePWMSpliceSite:
    input:
        "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.tab" 
    output:
        Tab = "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.scored.tab.gz",
        WeblogoImg = "SplicingAnalysis/NMDJunctions/Plots/{Junction}.PWM.png"
    log:
        "logs/NMDSpliceSite/{Junction}.PWM.log"
    wildcard_constraints:
        Junction = "FivePrime|ThreePrime",
    conda:
        "../envs/biopython.yml"
    shell:
        """
        python scripts/ScorePWM.py {input} {output.Tab} {output.WeblogoImg} &> {log}
        """
        
rule ScoreMaxEntScan:
    input:
        "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.tab" 
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/{Junction}.MaxEntScan_score.tab.gz",
    log:
        "logs/NMDSpliceSite/{Junction}.MaxEntScan.log"
    wildcard_constraints:
        Junction = "FivePrime|ThreePrime.ForMaxEntScan",
    conda:
        "../envs/maxentpy.yml"
    resources:
        mem_mb = 24000
    shell:
        """
        python scripts/ScoreMaxEntScan.py {input} {output} &> {log}
        """
    
    
    
    
    
    

        