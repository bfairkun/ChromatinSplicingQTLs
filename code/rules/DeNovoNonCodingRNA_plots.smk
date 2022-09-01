def GetBigWigFile(wildcards):
    location = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/chRNA.Expression.Splicing_stranded/" + wildcards.IndID

    if wildcards.strand == 'plus':
        bigwig_file = location + ".1.minus.bw"
    elif wildcards.strand == 'minus':
        bigwig_file = location + ".1.plus.bw"
        
    return bigwig_file
    

rule LogTransformBigWig:
    input:
        chrom_sizes = "../data/Chrome.sizes", 
        bigwig = GetBigWigFile,
    output:
        bedgraph = temp("NonCodingRNA_annotation/bigwig/chRNA.{IndID}.{strand}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA_annotation/bigwig/chRNA.{IndID}.{strand}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA_annotation/bigwig/chRNA.{IndID}.{strand}.log1p.bw"
    wildcard_constraints:
        IndID = "NA18486|NA19201|NA19210",
        strand = 'plus|minus'
    log:
        "logs/NonCodingRNA_annotation/chRNA.log1p.{IndID}.{strand}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bigWigToBedGraph {input.bigwig} {output.bedgraph}) &> {log};
        python scripts/LogTransformBigWig.py --input {output.bedgraph} --output {output.bedgraph_log} &>> {log};
        (bedGraphToBigWig {output.bedgraph_log} {input.chrom_sizes} {output.bigwig_log}) &>> {log};
        """

rule LogTransformPolyA:
    input:
        chrom_sizes = "../data/Chrome.sizes", 
        bigwig = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/Expression.Splicing/{IndID}.1.bw"
    output:
        bedgraph = temp("NonCodingRNA_annotation/bigwig/polyA.{IndID}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA_annotation/bigwig/polyA.{IndID}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA_annotation/bigwig/polyA.{IndID}.log1p.bw"
    wildcard_constraints:
        IndID = "NA18486|NA19201|NA19210",
    log:
        "logs/NonCodingRNA_annotation/polyA.log1p.{IndID}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bigWigToBedGraph {input.bigwig} {output.bedgraph}) &> {log};
        python scripts/LogTransformBigWig.py --input {output.bedgraph} --output {output.bedgraph_log} &>> {log};
        (bedGraphToBigWig {output.bedgraph_log} {input.chrom_sizes} {output.bigwig_log}) &>> {log};
        """




rule deepTools_chRNA_plotHeatmap:
    input:
        chRNA_plus = "NonCodingRNA_annotation/bigwig/chRNA.{IndID}.plus.log1p.bw",
        chRNA_minus = "NonCodingRNA_annotation/bigwig/chRNA.{IndID}.minus.log1p.bw",
        polyA = "NonCodingRNA_annotation/bigwig/polyA.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA_annotation/deeptools/matrix/{IndID}.RNA{strictness}.matrix.gz",
        png = "NonCodingRNA_annotation/deeptools/plots/{IndID}.RNA{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA_annotation/plots.{IndID}{strictness}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = "NA18486|NA19201|NA19210",
        strictness = "|.strict1|.strict2|.strict3"
    shell:
        """
        computeMatrix reference-point --regionBodyLength 5000 -a 3000 -b 3000 -R {input.plus_bed} {input.minus_bed} -S {input.chRNA_plus} {input.chRNA_minus} {input.polyA} -o {output.mat} &> {log}
        plotHeatmap -m {output.mat} -o {output.png} --regionsLabel "ncRNA.plus{wildcards.strictness}" "ncRNA.minus{wildcards.strictness}" --samplesLabel "chRNA plus {wildcards.IndID}" "chRNA minus {wildcards.IndID}" "polyA {wildcards.IndID}" --averageTypeSummaryPlot mean &>> {log}
        """
        

rule StrictFilter:
    input:
        ncPlus = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus.bed",
        ncMinus = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus.bed",
        genes = "NonCodingRNA_annotation/annotation/allGenes.bed.gz",
        chrom = "../data/Chrome.sizes"
    output:
        genes1 = "NonCodingRNA_annotation/deeptools/bed/allGenes_expanded1.bed.gz",
        genes2 = "NonCodingRNA_annotation/deeptools/bed/allGenes_expanded2.bed.gz",
        ncPlus1 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus.strict1.bed",
        ncMinus1 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus.strict1.bed",
        ncPlus2 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus.strict2.bed",
        ncMinus2 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus.strict2.bed",
        ncPlus3 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus.strict3.bed",
        ncMinus3 = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus.strict3.bed",
    log:
        "logs/NonCodingRNA_annotation/strict_filtering.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools slop -i {input.genes} -g {input.chrom} -b 5000 | awk -F'\t' '{{print $1"\\t"$2"\\t"$3}}' | gzip - > {output.genes1}) &> {log};
        (bedtools slop -i {input.genes} -g {input.chrom} -b 20000 | awk -F'\t' '{{print $1"\\t"$2"\\t"$3}}' | gzip - > {output.genes2}) &>> {log};
        (bedtools subtract -A -a {input.ncPlus} -b {output.genes1} -f 0.5 > {output.ncPlus1}) &>> {log};
        (bedtools subtract -A -a {input.ncMinus} -b {output.genes1} -f 0.5 > {output.ncMinus1}) &>> {log};
        (bedtools subtract -A -a {input.ncPlus} -b {output.genes2} -f 0.5 > {output.ncPlus2}) &>> {log};
        (bedtools subtract -A -a {input.ncMinus} -b {output.genes2} -f 0.5 > {output.ncMinus2}) &>> {log};
        (bedtools subtract -A -a {input.ncPlus} -b {output.genes1} > {output.ncPlus3}) &>> {log};
        (bedtools subtract -A -a {input.ncMinus} -b {output.genes1} > {output.ncMinus3}) &>> {log};
        """
        
        
################



use rule LogTransformPolyA as LogTransformHistone with:
    input:
        chrom_sizes = "../data/Chrome.sizes", 
        bigwig = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/bigwigs/{Phenotype}/{IndID}.1.bw"
    output:
        bedgraph = temp("NonCodingRNA_annotation/bigwig/{Phenotype}.{IndID}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA_annotation/bigwig/{Phenotype}.{IndID}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA_annotation/bigwig/{Phenotype}.{IndID}.log1p.bw"
    wildcard_constraints:
        IndID = "NA19201|NA19210",
        Phenotype = "H3K27AC|H3K4ME1|H3K4ME3|H3K36ME3|ProCap"
    log:
        "logs/NonCodingRNA_annotation/{Phenotype}.log1p.{IndID}.log"
        

rule deepTools_histone_plotHeatmap:
    input:
        procap = "NonCodingRNA_annotation/bigwig/ProCap.{IndID}.log1p.bw",
        h3k4me1 = "NonCodingRNA_annotation/bigwig/H3K4ME1.{IndID}.log1p.bw",
        h3k4me3 = "NonCodingRNA_annotation/bigwig/H3K4ME3.{IndID}.log1p.bw",
        h3k27ac = "NonCodingRNA_annotation/bigwig/H3K27AC.{IndID}.log1p.bw",
        h3k36me3 = "NonCodingRNA_annotation/bigwig/H3K36ME3.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA_annotation/deeptools/matrix/{IndID}.histone{strictness}.matrix.gz",
        png = "NonCodingRNA_annotation/deeptools/plots/{IndID}.histone{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA_annotation/histone_plots.{IndID}{strictness}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = "NA19201|NA19210",
        strictness = "|.strict1|.strict2|.strict3"
    shell:
        """
        computeMatrix reference-point --regionBodyLength 5000 -a 3000 -b 3000 -R {input.plus_bed} {input.minus_bed} -S {input.procap} {input.h3k27ac} {input.h3k4me1} {input.h3k4me3} {input.h3k36me3} -o {output.mat} &> {log}
        plotHeatmap -m {output.mat} -o {output.png} --regionsLabel "ncRNA.plus{wildcards.strictness}" "ncRNA.minus{wildcards.strictness}" --samplesLabel "ProCap {wildcards.IndID}" "H3K27ac {wildcards.IndID}" "H3K4me1 {wildcards.IndID}" "H3K4me3 {wildcards.IndID}" "H3K36me3 {wildcards.IndID}" --averageTypeSummaryPlot mean &>> {log}
        """
        
rule plotHeatmapAll:
    input:
        chRNA_plus = "NonCodingRNA_annotation/bigwig/chRNA.{IndID}.plus.log1p.bw",
        chRNA_minus = "NonCodingRNA_annotation/bigwig/chRNA.{IndID}.minus.log1p.bw",
        polyA = "NonCodingRNA_annotation/bigwig/polyA.{IndID}.log1p.bw",
        procap = "NonCodingRNA_annotation/bigwig/ProCap.{IndID}.log1p.bw",
        h3k4me1 = "NonCodingRNA_annotation/bigwig/H3K4ME1.{IndID}.log1p.bw",
        h3k4me3 = "NonCodingRNA_annotation/bigwig/H3K4ME3.{IndID}.log1p.bw",
        h3k27ac = "NonCodingRNA_annotation/bigwig/H3K27AC.{IndID}.log1p.bw",
        h3k36me3 = "NonCodingRNA_annotation/bigwig/H3K36ME3.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA_annotation/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA_annotation/deeptools/matrix/{IndID}.combined{strictness}.matrix.gz",
        png = "NonCodingRNA_annotation/deeptools/plots/{IndID}.combined{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA_annotation/combined_plots.{IndID}{strictness}.log"
    resources:
        mem_mb = 12000
    wildcard_constraints:
        IndID = "NA19201|NA19210",
        strictness = "|.strict1|.strict2|.strict3"
    shell:
        """
        computeMatrix reference-point --regionBodyLength 5000 -a 3000 -b 3000 -R {input.plus_bed} {input.minus_bed} -S {input.chRNA_plus} {input.chRNA_minus} {input.polyA} {input.procap} {input.h3k27ac} {input.h3k4me1} {input.h3k4me3} {input.h3k36me3} -o {output.mat} &> {log}
        plotHeatmap -m {output.mat} -o {output.png} --regionsLabel "ncRNA.plus{wildcards.strictness}" "ncRNA.minus{wildcards.strictness}" --samplesLabel "chRNA plus {wildcards.IndID}" "chRNA minus {wildcards.IndID}" "polyA {wildcards.IndID}" "ProCap {wildcards.IndID}" "H3K27ac {wildcards.IndID}" "H3K4me1 {wildcards.IndID}" "H3K4me3 {wildcards.IndID}" "H3K36me3 {wildcards.IndID}" --averageTypeSummaryPlot mean &>> {log}
        """
 
rule BamToBWWithIntronCoverage:
    input:
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = "Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam",
        bai = "Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam.bai"
    output:
        bw_plus = "NonCodingRNA_annotation/bigwig/{IndID}.plus.bw",
        bw_minus = "NonCodingRNA_annotation/bigwig/{IndID}.minus.bw",
    wildcard_constraints:
        IndID = "NA18486|NA18497|NA18498|NA18499|NA19201|NA19210"
    resources:
        mem_mb = 12000
    log:
        "logs/NonCodingRNA_annotation/{IndID}.bed2bigwig.log"
    shell:
        """
        bash scripts/BamToBigwig_PrefilterBySize.sh {input.fai} {input.bam} {output.bw_minus} bw_minus={output.bw_plus} SORT_ARGS='-T /scratch/midway2/cnajar/' MKTEMP_ARGS='-p /scratch/midway2/cnajar/'
        """
 