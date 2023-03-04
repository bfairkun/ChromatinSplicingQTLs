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
    log:
        "logs/Metaplots/split_bed12_into_quartiles.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/split_bed12_by_size.py &> {log}
        """
        
def GetMetaplotParamsMetagene(wildcards):
    if wildcards.metaplot == "CDS":
        return "--metagene"
    else:
        return ""
        
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
                        ".gene_length.Q1", ".gene_length.Q2", ".gene_length.Q3", ".gene_length.Q4"])
    resources:
        mem_mb = 12000
    shell:
        """
        (awk '$6=="+"' OFS='\\t' FS='\\t' {input} > {output.plus_bed12}) &> {log};
        (awk '$6=="-"' OFS='\\t' FS='\\t' {input} > {output.minus_bed12}) &> {log}
        """
        
rule ComputeMatrixForCoverageMetaplots:
    input:
        bigwigs = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}/{IndID}.1.bw",
        bed12 = "Metaplots/AssayProfiles/References/ExpressedGeneList.bed12",
    output:
        mat = "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.mat"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene"
    conda:
        "../envs/deeptools.yml"
    params:
        GetMetaplotParamsMetagene
    log:
        "logs/Metaplots/ComputeMatrix.{Phenotype}.{IndID}.{metaplot}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.bed12} -m 5000 -b 3000 -a 3000 -o {output} {params} &> {log}
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
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.png"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene"
    params:
        "{metaplot}"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/Metaplots/PlotHeatmap.{Phenotype}.{IndID}.{metaplot}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot mean --regionsLabel {params} --heatmapHeight 14 &>> {log}
        """
        
use rule plotHeatmapCoverageMetaplots as plotHeatmapCoverageMetaplots_stranded with:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{metaplot}.{strand}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{metaplot}.{strand}.png"
    wildcard_constraints:
        Phenotype =  "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        metaplot = "CDS|gene",
        strand = "plus|minus"
    log:
        "logs/Metaplots/PlotHeatmap.{Phenotype}.{IndID}.{metaplot}.{strand}.log"
    params:
        "{metaplot}.{strand}"
    
        
def GetMetaplotMatrixParams(wildcards):
    if wildcards.quartiles == 'CDS_length':
        return "--metagene"
    elif wildcards.quartiles == 'gene_length':
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
            "MetabolicLabelled.60min", "H3K36ME3"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length"
    log:
        "logs/Metaplots/ComputeMatrixQuartiles.{Phenotype}.{IndID}.{quartiles}.log"
    params:
        GetMetaplotMatrixParams
    conda:
        "../envs/deeptools.yml"
    resources:
        mem_mb = 12000
    shell:
        """
        computeMatrix scale-regions -S {input.bigwigs} -R {input.Q1} {input.Q2} {input.Q3} {input.Q4} -m 5000 -b 3000 -a 3000 -o {output} {params} &> {log}
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
        quartiles = "CDS_length|gene_length"
        
rule plotHeatmapCoverageMetaplotsQuartiles:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.png"
    wildcard_constraints:
        Phenotype = "|".join(["Expression.Splicing", "MetabolicLabelled.30min",
        "MetabolicLabelled.60min", "H3K36ME3"]),
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length"
    conda:
        "../envs/deeptools.yml"
    params:
        '"Q1" "Q2" "Q3" "Q4"'
    log:
        "logs/Metaplots/PlotHeatmap.Quartiles.{Phenotype}.{IndID}.{quartiles}.log"
    shell:
        """
        plotHeatmap -m {input} -o {output} --averageTypeSummaryPlot mean --regionsLabel {params} --heatmapHeight 14 &>> {log}
        """
        
use rule plotHeatmapCoverageMetaplotsQuartiles as plotHeatmapCoverageMetaplotsQuartiles_stranded with:
    input:
        "Metaplots/AssayProfiles/Matrix/{Phenotype}.{IndID}.{quartiles}.{strand}.mat"
    output:
        "Metaplots/AssayProfiles/Plots/{Phenotype}.{IndID}.{quartiles}.{strand}.png"
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing",
        IndID = "|".join(['NA18486', 'NA19137', 'NA19152', 'NA19153']),
        quartiles = "CDS_length|gene_length",
        strand = "plus|minus"
    log:
        "logs/Metaplots/PlotHeatmap.Quartiles.{Phenotype}.{IndID}.{quartiles}.{strand}.log"
        
        