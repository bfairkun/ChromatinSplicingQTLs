rule DownloadChromHMM:
    output:
        "ChromHMM/hg19.wgEncodeBroadHmmGm12878HMM.bed.gz"
    shell:
        """
        wget -O {output} http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
        """

rule liftOverChromHMM:
    input:
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz",
        bed = "ChromHMM/hg19.wgEncodeBroadHmmGm12878HMM.bed.gz"
    output:
        bed = "ChromHMM/hg38.wgEncodeBroadHmmGm12878HMM.bed",
        bedgz = "ChromHMM/hg38.wgEncodeBroadHmmGm12878HMM.bed.gz",
        tbi = "ChromHMM/hg38.wgEncodeBroadHmmGm12878HMM.bed.gz.tbi"
    shadow: "shallow"
    shell:
        """
        CrossMap.py bed {input.chain} {input.bed} {output.bed}
        cat {output.bed} | bedtools sort -i - | bgzip /dev/stdin -c > {output.bedgz}
        tabix -p bed {output.bedgz}
        """

rule Download_miRNA_BindingSites:
    output:
        "Misc/miRNA_binding_sites.hg19.bed"
    shell:
        """
        wget -O- --no-check-certificate "https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip" | zcat | awk -v OFS='\\t' -F'\\t' '{{print $1,$2,$3,"miRNA_BS", $4, $6}}' > {output}
        """

use rule liftOverChromHMM as liftOver_miRNA with:
    input:
        bed = "Misc/miRNA_binding_sites.hg19.bed",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz",
    output:
        bed = "Misc/miRNA_binding_sites.hg38.bed",
        bedgz = "Misc/miRNA_binding_sites.hg38.bed.gz",
        tbi = "Misc/miRNA_binding_sites.hg38.bed.gz.tbi"

rule CreatePAS_RegionAnnotation:
    input:
        bed = expand("QTLs/QTLTools/APA_{Fraction}/OnlyFirstReps.sorted.qqnorm.bed.gz", Fraction=["Nuclear", "Total"]),
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        "APA_Processing/hg38.TestFeature.SetWidthRegions.bed.gz"
    shell:
        """
        zcat {input.bed} | awk -v OFS='\\t' '$1 !~ "^#" {{print $1,$2,$3,$4,$5,$6}}' | sort | uniq | grep -P '^chr[0-9]+\\t' |  bedtools sort -i - | bedtools flank -s -l 0  -r 50 -g {input.fai} | bedtools shift -i - -m 50 -p -50 -g {input.fai} | bedtools slop -s -l 0 -r 10 -g {input.fai} | awk -v OFS='\\t' '{{print $1,$2,$3, "PAS_Region", $4, $6}}' | bedtools sort -i - | gzip - > {output}
        """


rule ConcatAllAnnotations:
    """
    annotation type in column 4. optional feature name in col 5.
    """
    input:
        "NonCodingRNA_annotation/annotation/ncRNA.categorized.bed.gz",
        "ChromHMM/hg38.wgEncodeBroadHmmGm12878HMM.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.3ss.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.5ss.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.bptregion.bed.gz",
        "APA_Processing/hg38.TestFeature.SetWidthRegions.bed.gz",
        "Misc/miRNA_binding_sites.hg38.bed.gz"
    output:
        bed = "QTL_SNP_Enrichment/Annotations.bed.gz",
        tbi = "QTL_SNP_Enrichment/Annotations.bed.gz.tbi"
    shell:
        """
        zcat {input} | awk -F'\\t' -v OFS='\\t' '{{ print $1,$2,$3,$4,$5,$6 }}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule IntersectFinemapSNPsWithAnnotations:
    """
    Bed should be a bed with all genomic region annotations, with the annotation type in column4. Optional feature name in column5. Finemap should be the hyprcoloc finemap output.
    """
    input:
        bed = "QTL_SNP_Enrichment/Annotations.bed.gz",
        Finemap = "hyprcoloc/Results/ForColoc/{ColocRun}/snpscores.txt.gz"
    output:
        "QTL_SNP_Enrichment/FinemapIntersections/{ColocRun}.bed.gz"
    log:
        "logs/IntersectFinemapSNPsWithAnnotations/{ColocRun}.log"
    shell:
        """
        zcat {input.Finemap} | awk -v OFS='\\t' -F'\\t' 'NR>1 {{split($1,snp,":"); print "chr"snp[1], snp[2], snp[2]+1, $1"_"$2"_"$4, $3}}' | bedtools sort -i - | bedtools intersect -a - -b <(zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $3, $4 }}')  -wao | gzip - > {output}
        """ 

def GetAwkCommandForParseQTLToolsPermutationPass(wildcards):
    if wildcards.Phenotype in ["polyA.Splicing", "chRNA.Splicing", "polyA.Splicing.Subset_YRI", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing", "chRNA.RNA.Editing", "chRNA.Splicing.Order"]:
        return '$11, $12, $13, $10, $1, ".", $6";"$3";"$4";"$20";"$23";"$24'
    else:
        return '$9, $10, $11, $8, ".", ".", $1";"$3";"$4";"$18";"$21";"$22'

rule IntersectTopSNPsFromPermutationPassWithAnnotations:
    input:
        bed = "QTL_SNP_Enrichment/Annotations.bed.gz",
        PermutationPass = "QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz",
    output:
        "QTL_SNP_Enrichment/TopSNPIntersections/{Phenotype}.bed.gz"
    log:
        "logs/IntersectTopSNPsFromPermutationPassWithAnnotations/{Phenotype}.log"
    params:
        GetAwkCommandForParseQTLToolsPermutationPass
    shell:
        """
        (zcat {input.PermutationPass} | awk -F' ' -v OFS='\\t' 'NR>1 {{print {params} }}' |  grep -v '^NA' | bedtools sort -i - | bedtools intersect -a - -b <(zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,$4}}' ) -wao | gzip - > {output} ) &> {log}
        """

rule IntersectTopSNPsFromPermutationPass_ColocWindow_WithAnnotations:
    input:
        bed = "QTL_SNP_Enrichment/Annotations.bed.gz",
        PermutationPass = "QTLs/QTLTools/{Phenotype}/PermutationPassForColoc.txt.gz",
    output:
        "QTL_SNP_Enrichment/TopSNPIntersections_ForColocWindows/{Phenotype}.bed.gz"
    log:
        "logs/IntersectTopSNPsFromPermutationPass_ColocWindow_WithAnnotations/{Phenotype}.log"
    params:
        GetAwkCommandForParseQTLToolsPermutationPass = '$9, $10, $11, $8, ".", ".", $1";"$3";"$4";"$18";"$21";"$22'
    shell:
        """
        (zcat {input.PermutationPass} | awk -F' ' -v OFS='\\t' 'NR>1 {{print {params} }}' |  grep -v '^NA' | bedtools sort -i - | bedtools intersect -a - -b <(zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,$4}}' ) -wao | gzip - > {output} ) &> {log}
        """

# rule Classify_eQTLs_AsTranscriptionalOrPosttranscriptional:
#     input:
#         pi1 = "pi1/PairwisePi1Traits.P.all.txt.gz",
#         TSS = expand("Misc/PeaksClosestToTSS/{Phenotype}_assigned.tsv.gz", Phenotype=["H3K27AC", "H3K4ME3"])
#     output:
#         "QTL_SNP_Enrichment/PostTranscriptionalVsTranscriptional_eQTLs/List.tsv.gz"
#     conda:
#         "../envs/r_2.yaml"
#     shell:
#         """
#         Rscript scripts/ClassifyTxnVsPostTxnEqtls.R {input.pi1}
#         """

rule IntersectFinemapSNPsWithAnnotations_susie:
    """
    Bed should be a bed with all genomic region annotations, with the annotation type in column4. Optional feature name in column5. Finemap should be the hyprcoloc finemap output.
    """
    input:
        bed = "QTL_SNP_Enrichment/Annotations.bed.gz",
        Finemap = "FineMapping/susie_runs_{Pop}/susie_output.tab"
    output:
        "QTL_SNP_Enrichment/FinemapIntersections_Susie/{Pop}.bed.gz"
    log:
        "logs/IntersectFinemapSNPsWithAnnotations_Susie/{Pop}.log"
    shell:
        """
        awk -v OFS='\\t' -F'\\t' 'NR>1 && $2!="NA" {{split($2,snp,":"); print "chr"snp[1], snp[2], snp[2]+1, $2"_"$1"_"$3, $4}}' {input.Finemap} | bedtools sort -i - | bedtools intersect -a - -b <(zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $3, $4 }}')  -wao | gzip - > {output}
        """ 

rule GatherFinemapSNPAnnotationIntersections_susie:
    input:
        expand("QTL_SNP_Enrichment/FinemapIntersections_Susie/{Pop}.bed.gz", Pop=["YRI", "Geuvadis"])

rule GatherFinemapSNPAnnotationIntersections:
    input:
        expand("QTL_SNP_Enrichment/FinemapIntersections/{ColocRun}.bed.gz", ColocRun = colocs_genewise.index)

rule GatherIntersectTopSNPsFromPermutationPassWithAnnotations:
    input:
        expand("QTL_SNP_Enrichment/TopSNPIntersections_ForColocWindows/{Phenotype}.bed.gz", Phenotype = PhenotypesToColoc)
