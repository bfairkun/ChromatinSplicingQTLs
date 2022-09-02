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
    input:
        "ChromHMM/hg38.wgEncodeBroadHmmGm12878HMM.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.3ss.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.5ss.bed.gz",
        "SplicingAnalysis/regtools_annotate_combined/comprehensive.bptregion.bed.gz",
        "APA_Processing/hg38.TestFeature.SetWidthRegions.bed.gz"
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

rule GatherFinemapSNPAnnotationIntersections:
    input:
        expand("QTL_SNP_Enrichment/FinemapIntersections/{ColocRun}.bed.gz", ColocRun = colocs_genewise.index)
