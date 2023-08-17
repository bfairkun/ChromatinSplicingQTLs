rule UnzipBed12ForMetaplot:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.bed12.bed.gz"
    output:
        "Metaplots/AssayProfiles/References/ExpressedGeneList.bed12",
    log:
        "logs/Metaplots/unzip_bed12.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (zcat {input} > {output}) &> {log}
        """
        
        
rule SplitBed12ByLengthQuartile:
    input:
        "Metaplots/AssayProfiles/References/ExpressedGeneList.bed12"
    output:
        "Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q1.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q2.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q3.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.CDS_length.Q4.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q1.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q2.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q3.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.gene_length.Q4.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.expression.Q1.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.expression.Q2.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.expression.Q3.bed12",
        "Metaplots/AssayProfiles/References/ExpressedGeneList.expression.Q4.bed12"
    log:
        "logs/Metaplots/split_bed12_into_quartiles.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/split_bed12.py &> {log}
        """
        
        
rule splitBed12:
    input:
        "Metaplots/AssayProfiles/References/ExpressedGeneList{bed12}.bed12"
    output:
        plus_bed12 = "Metaplots/AssayProfiles/References/ExpressedGeneList{bed12}.plus.bed12",
        minus_bed12 = "Metaplots/AssayProfiles/References/ExpressedGeneList{bed12}.minus.bed12"
    log:
        "logs/Metaplots/split_bed.{bed12}.log"
    wildcard_constraints:
        bed12="|".join(["", ".CDS_length.Q1", ".CDS_length.Q2", ".CDS_length.Q3", ".CDS_length.Q4",
                        ".gene_length.Q1", ".gene_length.Q2", ".gene_length.Q3", ".gene_length.Q4",
                        ".expression.Q1", ".expression.Q2", ".expression.Q3", ".expression.Q4"])
    resources:
        mem_mb = 12000
    shell:
        """
        (awk '$6=="+"' OFS='\\t' FS='\\t' {input} > {output.plus_bed12}) &> {log};
        (awk '$6=="-"' OFS='\\t' FS='\\t' {input} > {output.minus_bed12}) &> {log}
        """
        
gene_assays = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "H3K36ME3", "chRNA.Expression.Splicing"]
tss_assays = ["H3K27AC", "H3K4ME3", "H3K4ME1"]
        
def MatrixType(wildcards):
    if wildcards.Phenotype in gene_assays:
        return "scale-regions"
    else:
        return "reference-point"
        
#def IsMetaplotMetagene(wildcards):
#    if wildcards.metaplot == "CDS":
#        return "--metagene"
#    else:
#        return ""
        
def GetMetaplotParams(wildcards):
    if wildcards.Phenotype in gene_assays:
        return '-m 5000 -b 3000 -a 3000'
    else:
        return "-b 3000 -a 3000"
        
        
rule ComputeMatrixForCoverageMetaplots:
    input:
        bigwigs = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}/{IndID}.1.bw",
        bed12 = "Metaplots/AssayProfiles/References/ExpressedGeneList.bed12",
    output:
        mat = "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.mat"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene"
    conda:
        "../envs/deeptools.yml"
    params:
        MatrixType = MatrixType,
        MetaplotParams = GetMetaplotParams 
    log:
        "logs/Metaplots/ComputeMatrix.{Phenotype}.{IndID}.{metaplot}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix {params.MatrixType} -S {input.bigwigs} -R {input.bed12} {params.MetaplotParams} -o {output} &> {log}
        """
        
def SwitchStrandBW(wildcards):
    if wildcards.strand == "plus":
        return "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}_stranded/{IndID}.1.minus.bw"
    else:
        return "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}_stranded/{IndID}.1.plus.bw"
        
use rule ComputeMatrixForCoverageMetaplots as ComputeMatrixForCoverageMetaplots_stranded with:
    input:
        bigwigs = SwitchStrandBW,
        bed12 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{strand}.bed12",
    output:
        mat = "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.{strand}.mat"
    log:
        "logs/Metaplots/ComputeMatrix.{Phenotype}.{IndID}.{metaplot}.{strand}.log"
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene",
        strand = "plus|minus"
        
        
rule plotHeatmapCoverageMetaplots:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.{stat}.png"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene",
        stat = 'median|mean'
    params:
        "{metaplot}"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/Metaplots/PlotHeatmap.{Phenotype}.{IndID}.{metaplot}.{stat}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot {wildcards.stat} --regionsLabel {params} --heatmapHeight 10 &>> {log}
        """
        
use rule plotHeatmapCoverageMetaplots as plotHeatmapCoverageMetaplots_stranded with:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.{strand}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.{strand}.{stat}.png"
    wildcard_constraints:
        Phenotype =  "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene",
        strand = "plus|minus",
        stat = 'median|mean'
    log:
        "logs/Metaplots/PlotHeatmap.{Phenotype}.{IndID}.{metaplot}.{strand}.{stat}.log"
    params:
        "{metaplot}.{strand}"
    
        
def MetageneParams(wildcards):
    if wildcards.quartiles == 'CDS_length':
        return "--metagene"
    else:
        return ""
        
rule ComputeMatrixForCoverageMetaplots_quartiles:
    input:
        bigwigs = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}/{IndID}.1.bw",
        Q1 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q1.bed12",
        Q2 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q2.bed12",
        Q3 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q3.bed12",
        Q4 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q4.bed12",
    output:
        mat = "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.mat"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
            "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length|expression"
    log:
        "logs/Metaplots/ComputeMatrixQuartiles.{Phenotype}.{IndID}.{quartiles}.log"
    params:
        MatrixType = MatrixType,
        MetaplotParams = GetMetaplotParams,
        MetageneParams = MetageneParams
    conda:
        "../envs/deeptools.yml"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix {params.MatrixType} -S {input.bigwigs} -R {input.Q4} {input.Q3} {input.Q2} {input.Q1} {params.MetaplotParams} {params.MetageneParams} -o {output} &> {log}
        """

use rule ComputeMatrixForCoverageMetaplots_quartiles as ComputeMatrixForCoverageMetaplots_quartiles_stranded with:
    input:
        bigwigs = SwitchStrandBW,
        Q1 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q1.{strand}.bed12",
        Q2 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q2.{strand}.bed12",
        Q3 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q3.{strand}.bed12",
        Q4 = "Metaplots/AssayProfiles/References/ExpressedGeneList.{quartiles}.Q4.{strand}.bed12",
    output:
        mat = "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.{strand}.mat"
    log:
        "logs/Metaplots/ComputeMatrix.{Phenotype}.{IndID}.{quartiles}.{strand}.log"
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        strand = "plus|minus",
        quartiles = "CDS_length|gene_length|expression"
        
rule plotHeatmapCoverageMetaplotsQuartiles:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.{stat}.png"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3", "H3K27AC", "H3K4ME3", "H3K4ME1"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length|expression",
        stat = 'median|mean'
    conda:
        "../envs/deeptools.yml"
    params:
        '"Q4" "Q3" "Q2" "Q1"'
    log:
        "logs/Metaplots/PlotHeatmap.Quartiles.{Phenotype}.{IndID}.{quartiles}.{stat}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot {wildcards.stat} --regionsLabel {params} --heatmapHeight 10 &>> {log}
        """
        
use rule plotHeatmapCoverageMetaplotsQuartiles as plotHeatmapCoverageMetaplotsQuartiles_stranded with:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.{strand}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.{strand}.{stat}.png"
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length|expression",
        strand = "plus|minus",
        stat = 'median|mean'
    log:
        "logs/Metaplots/PlotHeatmap.Quartiles.{Phenotype}.{IndID}.{quartiles}.{strand}.{stat}.log"
        

def GetStrandParams(wildcards):
    if wildcards.strand == 'plus':
        return '"+"'
    elif wildcards.strand == 'minus':
        return '"-"'
        
def GetIntronLengAwkParams(wildcards):
    if wildcards.intron_length == ".long_introns":
        return " && (($3-$2)>10000)"
    else:
        return ""

rule GetIntronsBed:
    input:
        "QTLs/QTLTools/chRNA.IER/OnlyFirstReps.qqnorm.bed.gz"
    output:
        "Metaplots/AssayProfiles/References/spliceq.introns.{strand}{intron_length}.bed"
    log:
        "logs/Metaplots/spliceq.introns.{strand}{intron_length}.log"
    resources:
        mem_mb = 12000
    params:
        strand = GetStrandParams,
        intron_len_params = GetIntronLengAwkParams
    wildcard_constraints:
        strand = 'plus|minus',
        intron_length = '|.long_introns'
    shell:
        """
        (zcat {input} | awk '(NR>1 && $6=={params.strand}{params.intron_length})  {{print $1, $2, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' | bedtools sort -i - > {output}) &> {log}
        """
        
#rule GetIntronsBed_long:
#    input:
#        "QTLs/QTLTools/chRNA.IER/OnlyFirstReps.qqnorm.bed.gz"
#    output:
#        "Metaplots/AssayProfiles/References/spliceq.introns.{strand}.long_introns.bed"
#    log:
#        "logs/Metaplots/spliceq.introns.{strand}.log"
#    resources:
#        mem_mb = 12000
#    wildcard_constraints:
#        strand = 'plus|minus'
#    params:
#        GetStrandParams
#    shell:
#        """
#        (zcat {input} | awk '(NR>1 && $6=={params} && (($3-$2)>10000))  {{print $1, $2, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' | bedtools sort -i - > {output}) &> {log}
#        """
    
def GetChRNASample(wildcards):
    if wildcards.strand == 'plus':
        return "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/chRNA.Expression.Splicing_stranded/{IndID}.1.minus.bw",
    elif wildcards.strand == 'minus':
        return "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/chRNA.Expression.Splicing_stranded/{IndID}.1.plus.bw",

rule ComputeMatrixForIntrons_chRNA:
    input:
        chRNA = GetChRNASample,
        ml30 = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/MetabolicLabelled.30min/{IndID}.1.bw",
        ml60 = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/MetabolicLabelled.60min/{IndID}.1.bw",
        polyA = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/Expression.Splicing/{IndID}.1.bw",
        bed = "Metaplots/AssayProfiles/References/spliceq.introns.{strand}{intron_length}.bed"
    output:
        "Metaplots/AssayProfiles/Matrix/Introns.{IndID}.{strand}{intron_length}.mat"
    log:
        "logs/Metaplots/ComputeMatrix.splicing.chRNA.Expression.Splicing.{IndID}.{strand}{intron_length}.log"
    conda:
        "../envs/deeptools.yml"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        strand = 'plus|minus',
        intron_length = '|.long_introns'
    shell:
        """
        computeMatrix scale-regions -S {input.chRNA} {input.ml30} {input.ml60} {input.polyA} -R {input.bed} -m 600 -b 200 -a 200 -o {output} &> {log}
        """
        
#rule ComputeMatrixForIntrons_long:
#    input:
#        chRNA = GetChRNASample,
#        ml30 = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/MetabolicLabelled.30min/{IndID}.1.bw",
#        ml60 = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/MetabolicLabelled.60min/{IndID}.1.bw",
#        polyA = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/Expression.Splicing/{IndID}.1.bw",
#        bed = "Metaplots/AssayProfiles/References/spliceq.introns.{strand}.long_introns.bed"
#    output:
#        "Metaplots/AssayProfiles/Matrix/Introns.{IndID}.{strand}.long_introns.mat"
#    log:
#        "logs/Metaplots/ComputeMatrix.splicing.chRNA.Expression.Splicing.{IndID}.{strand}.long_introns.log"
#    conda:
#        "../envs/deeptools.yml"
#    resources:
#        mem_mb = 12000
#    wildcard_constraints:
#        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
#        strand = 'plus|minus'
#    shell:
#        """
#        computeMatrix scale-regions -S {input.chRNA} {input.ml30} {input.ml60} {input.polyA} -R {input.bed} -m 600 -b 200 -a 200 -o {output} &> {log}
#        """

#rule ComputeMatrixForIntrons:
#    input:
#        bigwigs = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}/{IndID}.1.bw",
#        bed = "Metaplots/AssayProfiles/References/spliceq.introns.{strand}.bed"
#    output:
#        "Metaplots/AssayProfiles/Matrix/Introns.{Phenotype}.{IndID}.{strand}.mat"
#    log:
#        "logs/Metaplots/ComputeMatrix.splicing.{Phenotype}.{IndID}.{strand}.log"
#    conda:
#        "../envs/deeptools.yml"
#    wildcard_constraints:
#        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min"]),
#        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
#    resources:
#        mem_mb = 12000
#    shell:
#        """
#        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed} -m 500 -b 100 -a 100 -o {output} &> {log}
#        """
        
rule plotHeatmapMetaplots_Introns:
    input:
        "Metaplots/AssayProfiles/Matrix/Introns.{IndID}.{strand}{intron_length}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/Introns.{IndID}.{strand}{intron_length}.png"
    wildcard_constraints:
        strand = 'plus|minus',
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        intron_length = '|.long_introns'
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/Metaplots/PlotHeatmap.Introns.{IndID}.{strand}{intron_length}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot mean --heatmapHeight 14 &>> {log}
        """
        
#rule plotHeatmapMetaplots_Introns_long:
#    input:
#        "Metaplots/AssayProfiles/Matrix/Introns.{IndID}.{strand}.long_introns.mat"
#    output:
#        "Metaplots/AssayProfiles/Plots/Introns.{IndID}.{strand}.long_introns.png"
#    wildcard_constraints:
#        strand = 'plus|minus',
#        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
#    conda:
#        "../envs/deeptools.yml"
#    log:
#        "logs/Metaplots/PlotHeatmap.Introns.{IndID}.{strand}.log"
#    shell:
#        """
#        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot mean --heatmapHeight 14 &>> {log}
#        """


