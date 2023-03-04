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
        
#rule CollectScoreOutput:
#    input:
#        annotation = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed.gz",
#        FivePrimePWM = "SplicingAnalysis/NMDJunctions/Annotation/FivePrime.scored.tab.gz",
#        ThreePrimePWM = "SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.scored.tab.gz",
#        FivePrimeMES = "SplicingAnalysis/NMDJunctions/Annotation/FivePrime.MaxEntScan_score.tab.gz",
#        ThreePrimeMES = #"SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.ForMaxEntScan.MaxEntScan_score.tab.gz",
#    output:
#        annotation_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed"),
#        FivePrimePWM_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/FivePrime.scored.tab"),
#        ThreePrimePWM_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.scored.tab"),
#        FivePrimeMES_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/FivePrime.MaxEntScan_score.tab"),
#        ThreePrimeMES_tmp = #temp("SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.ForMaxEntScan.MaxEntScan_score.tab"),
#        annotation_scored_1 = #temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_1.bed"),
#        annotation_scored_2 = #temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_2.bed"),
#        annotation_scored_3 = #temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_3.bed"),
#        annotation_scored = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored.bed.gz"
#    log:
#        "logs/NMDSpliceSite/CollectScores.log"
#    resources:
#        mem_mb = 62000
#    shell:
#        """
#        echo "Uncompressing..." > {log};
#        (zcat {input.annotation} > {output.annotation_tmp}) &>> {log};
#        (zcat {input.FivePrimePWM} > {output.FivePrimePWM_tmp}) &>> {log};
#        (zcat {input.ThreePrimePWM} > {output.ThreePrimePWM_tmp}) &>> {log};
#        (zcat {input.FivePrimeMES} > {output.FivePrimeMES_tmp}) &>> {log};
#        (zcat {input.ThreePrimeMES} > {output.ThreePrimeMES_tmp}) &>> {log};
#        echo "Done uncompressing. Merging 5'SS PWM..." >> {log};
#        (paste {output.annotation_tmp} {output.FivePrimePWM_tmp} | cut --complement -f 10 - > #{output.annotation_scored_1}) &>> {log};
#        echo "Merging 3'SS PWM..." >> {log};
#        (paste {output.annotation_scored_1} {output.ThreePrimePWM_tmp} | cut --complement -f 12 - > #{output.annotation_scored_2}) &>> {log};
#        echo "Merging 5'SS MaxEntScan..." >> {log};
#        (paste {output.annotation_scored_2} {output.FivePrimeMES_tmp} | cut --complement -f 14 - > #{output.annotation_scored_3}) &>> {log};
#        echo "Merging 3'SS MaxEntScan..." >> {log};
#        (paste {output.annotation_scored_3} {output.ThreePrimeMES_tmp} | cut --complement -f 16 - | gzip - > #{output.annotation_scored}) &>> {log};
#        echo "Done!" >> {log};
#        """
        
        
rule CollectScoreOutput_1:
    input:
        annotation = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed.gz",
        score = "SplicingAnalysis/NMDJunctions/Annotation/FivePrime.scored.tab.gz",
    output:
        annotation_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.bed"),
        score_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/FivePrime.scored.tab"),
        annotation_scored = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_1.bed.gz"),
    log:
        "logs/NMDSpliceSite/CollectScores_1.log"
    resources:
        mem_mb = 32000
    shell:
        """
        echo "Uncompressing..." > {log};
        (zcat {input.annotation} > {output.annotation_tmp}) &>> {log};
        (zcat {input.score} > {output.score_tmp}) &>> {log};
        echo "Done uncompressing. Merging score..." >> {log};
        (paste {output.annotation_tmp} {output.score_tmp} | gzip - > {output.annotation_scored}) &>> {log};
        echo "Done!" >> {log};
        """
        
use rule CollectScoreOutput_1 as CollectScoreOutput_2 with:
    input:
        annotation = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_1.bed.gz",
        score = "SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.scored.tab.gz",
    output:
        annotation_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_1.bed"),
        score_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.scored.tab"),
        annotation_scored = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_2.bed.gz"),
    log:
        "logs/NMDSpliceSite/CollectScores_2.log"
        
use rule CollectScoreOutput_1 as CollectScoreOutput_3 with:
    input:
        annotation = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_2.bed.gz",
        score = "SplicingAnalysis/NMDJunctions/Annotation/FivePrime.MaxEntScan_score.tab.gz",
    output:
        annotation_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_1.bed"),
        score_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/FivePrime.MaxEntScan_score.tab"),
        annotation_scored = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_3.bed.gz"),
    log:
        "logs/NMDSpliceSite/CollectScores_3.log"
        
rule CollectScoreOutput_4:
    input:
        annotation = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_3.bed.gz",
        score = "SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.ForMaxEntScan.MaxEntScan_score.tab.gz",
    output:
        annotation_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored_3.bed"),
        score_tmp = temp("SplicingAnalysis/NMDJunctions/Annotation/ThreePrime.ForMaxEntScan.MaxEntScan_score.tab"),
        annotation_scored = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored.bed.gz",
    log:
        "logs/NMDSpliceSite/CollectScores_4.log"
    resources:
        mem_mb = 48000
    shell:
        """
        echo "Uncompressing..." > {log};
        (zcat {input.annotation} > {output.annotation_tmp}) &>> {log};
        (zcat {input.score} > {output.score_tmp}) &>> {log};
        echo "Done uncompressing. Merging score..." >> {log};
        (paste {output.annotation_tmp} {output.score_tmp} | cut --complement -f 10,13,16,19 - | gzip - > {output.annotation_scored}) &>> {log};
        echo "Done!" >> {log};
        """
        
rule AddIDToScoredAnnotation:
    input:
        annot_junc = "SplicingAnalysis/NMDJunctions/Annotation/annotation_leaf_JAN28.scored.bed.gz",
        clusters = "/project2/yangili1/yangili/chRNA/intron_clusters.txt"
    output:
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_junctions.bed.gz",
        "SplicingAnalysis/NMDJunctions/Annotation/annotation_clusters.bed.gz",
    log:
        "logs/NMDSpliceSite/add_info_to_annotation.log"
    shell:
        """
        python scripts/add_additional_annotation_to_junctions.py {input.annot_junc} {input.clusters} &> {log}
        """
        
rule GetHistone_peaks_close_to_gene:
    input:
        counts = "featureCounts/{histone}/Counts.txt",
        tss = "ReferenceGenome/Annotations/GTFTools_BasicAnnotations/gencode.v34.chromasomal.tss.bed"
    output:
        temp_counts = temp('QTLs/QTLTools/{histone}/featureCounts.bed'),
        peaks_tss = 'QTLs/QTLTools/{histone}/CountsPeaksAtTSS.bed.gz'
    wildcard_constraints:
        histone = 'H3K4ME3|H3K27AC'
    log:
        "logs/Get{histone}_peaks_close_to_gene.log"
    resources:
        mem_mb = 62000
    shell:
        """
        (tail -n+3 featureCounts/{wildcards.histone}/Counts.txt | awk -v n=7 '{{print $2, $3, $4, $1, $6, $5, $0}}' FS='\\t' OFS='\\t' - | cut -f 7,8,9,10,11,12 --complement - > {output.temp_counts}) &> {log};
        (sed -e 's/^/chr/' ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed | awk '{{print $1, $2, $3, $5, $6, $4, $7}}' OFS='\\t' FS='\\t' - | bedtools sort -i - | bedtools closest -b {output.temp_counts} -a - -d | awk '$92>=0 && $92<=2000' OFS='\\t' FS='\\t' - | gzip - > {output.peaks_tss}) &>> {log}
        """
            

        