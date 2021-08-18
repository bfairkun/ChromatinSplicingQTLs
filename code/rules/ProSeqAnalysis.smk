rule CollectProSeqData:
    input:
        "ProCapAnalysis/SupplementDat2.csv",
        expand("ProCapAnalysis/RefSeqClassification_{tTRE_type}.bed", tTRE_type=["distal_enhancer", "promoter", "proximal_enhancer"]),
        # "ProCapAnalysis/GM18505_PRO-seq.{strand}.bw",
        expand("ProCapAnalysis/deeptools.{libtype}.plot.png", libtype=["Pro-seq", "chRNA", "polyA-RNA"])


rule DownloadKristjansdottirSupplementData2:
    output:
        "ProCapAnalysis/SupplementDat2.csv"
    shell:
        """
        wget -O {output} https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19829-z/MediaObjects/41467_2020_19829_MOESM6_ESM.csv
        """

rule ProSeqCountTable:
    input:
        SupDat = "ProCapAnalysis/SupplementDat2.csv",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        "ProCapAnalysis/ProCap.CountTable.hg38.bed"
    shadow: "shallow"
    shell:
        """
        awk -F, -v OFS=',' 'NR>1 && $8=="TRUE" {{print $1"_"$2, $0}}' {input.SupDat} | cut -d, --complement -f2-9 | sed 's/,/\\t/g' | sort > ProCapAnalysis/CountTable.values.tsv
        awk -F, -v OFS='\\t' 'NR>1 && $8=="TRUE" {{print $1, $2, $2+1, $1"_"$2, $5, $4}}' {input.SupDat} | sed 's/,/\\t/g' | CrossMap.py bed {input.chain} /dev/stdin ProCapAnalysis/CountTable.hg38.features.bed
        join --nocheck-order <(awk -F'\\t' -v OFS='\\t' '{{print $4, $0}}' ProCapAnalysis/CountTable.hg38.features.bed | sort) <(cat ProCapAnalysis/CountTable.values.tsv | sort) > ProCapAnalysis/Output.headerless
        cat <(cat {input.SupDat} | awk -F, -v OFS="," 'NR==1 {{ print "pid", "#Chr", "start", "end", "pid", "gid", "strand", $0 }}' | sed 's/,/\\t/g' | cut --complement -f8-15) <(cat ProCapAnalysis/Output.headerless | sed 's/ /\\t/g' ) | cut --complement -f1 | awk -F'\\t' -v OFS='\\t' 'NR==1 {{print}} NR>1 {{$4=$5"_"$1"_"$2; $5=$4;  print $0}}' > {output}
        """


rule GetRefSeqtTRE_beds:
    input:
        SupDat = "ProCapAnalysis/SupplementDat2.csv",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        "ProCapAnalysis/RefSeqClassification_{tTRE_type}.bed"
    shell:
        """
        awk -F, -v OFS='\\t' 'NR>1 && $9==0 && $5=="{wildcards.tTRE_type}" {{print $1, $2, $2+1, ".", ".", $4}}' {input.SupDat} | CrossMap.py bed {input.chain} /dev/stdin {output}
        """

def GetProSeqExample_DownloadLink(wildcards):
    if wildcards.strand == "minus":
        return "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004714/suppl/GSM3004714%5F48651%5FR1%5Fcap%2Eclip%2ETTAGGC%2Ehs37d5%2Ebwa%2EuniqueUMI%2Emn%2Ebedgraph%2Egz"
    elif wildcards.strand == "plus":
        return "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004714/suppl/GSM3004714%5F48651%5FR1%5Fcap%2Eclip%2ETTAGGC%2Ehs37d5%2Ebwa%2EuniqueUMI%2Epl%2Ebedgraph%2Egz"
rule GetProSeq_bigwig:
    """
    hg18 bedgraphs on GEO, download and convert to hg38 bedgraph
    """
    input:
        SupDat = "ProCapAnalysis/SupplementDat2.csv",
        chromSizes = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        bg = "ProCapAnalysis/GM18505_PRO-seq.{strand}.bg",
    shadow: "shallow"
    wildcard_constraints:
        strand = "plus|minus"
    params:
        Link = GetProSeqExample_DownloadLink
    shell:
        """
        wget -O hg19.bedgraph.gz {params.Link}
        CrossMap.py bed {input.chain} hg19.bedgraph.gz hg38.bg
        bedtools sort -i hg38.bg | grep -P '^chr\\d+\\t' > {output.bg}
        """

rule MergeProSeqIntoUnstrandedBigwig:
    input:
        bg = expand("ProCapAnalysis/GM18505_PRO-seq.{strand}.bg", strand=["plus", "minus"]),
        chromSizes = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
    output:
        bw = "ProCapAnalysis/GM18505_Pro-seq.unstranded.bw"
    shadow: "shallow"
    shell:
        """
        bedtools unionbedg -i {input.bg} -empty -g {input.chromSizes} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,$4+$5*-1}}' > union.bedgraph
        bedGraphToBigWig union.bedgraph {input.chromSizes} {output.bw}
        """

rule NormalizeProCapPhenotypeTable:
    input:
        "ProCapAnalysis/ProCap.CountTable.hg38.bed"
    output:
        "QTLs/QTLTools/ProCap/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/NormalizeProCapPhenotypeTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTable_ProCap.R {input} {output} &> {log}
        """

def GetBigwigForDeeptToolseRNA(wildcards):
    if wildcards.libtype == "Pro-seq":
        return "ProCapAnalysis/GM18505_Pro-seq.unstranded.bw"
    elif wildcards.libtype == "chRNA":
        return "bigwigs/chRNA.Expression.Splicing/NA18497.1.bw"
    elif wildcards.libtype == "polyA-RNA":
        return "bigwigs/Expression.Splicing/NA18505.1.bw"

rule DownloadEnhancerAtlas:
    input:
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        "ProCapAnalysis/eRNA_18505.bed"
    shell:
        """
        wget -O-  http://www.enhanceratlas.org/data/download/enhancer/hs/GM18505.bed | CrossMap.py bed {input.chain} /dev/stdin {output}
        """

rule DownloadSlideBaseBCellEnhancers:
    output:
        "ProCapAnalysis/enhancers_LCL.bed"
    shell:
        "wget -O {output} https://slidebase.binf.ku.dk/human_enhancers/presets/serve/lymphocyte_of_B_lineage"

rule deepTools_eRNA_heatmap:
    input:
        bigwigs = GetBigwigForDeeptToolseRNA,
        bedsFromProCapTable = expand("ProCapAnalysis/RefSeqClassification_{tTRE_type}.bed", tTRE_type=["distal_enhancer", "promoter", "proximal_enhancer"]),
        expressed_genes_bed = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        ehancer_map_bed = "ProCapAnalysis/enhancers_LCL.bed"
    output:
        mat = "ProCapAnalysis/deeptools.{libtype}.matrix.gz",
        png = "ProCapAnalysis/deeptools.{libtype}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/deepTools_eRNA_heatmap.{libtype}.log"
    shell:
        """
        computeMatrix reference-point -a 1500 -b 1500 -R {input.bedsFromProCapTable} {input.expressed_genes_bed} {input.ehancer_map_bed} -S {input.bigwigs} -o {output.mat} &> {log}
        plotHeatmap -m {output.mat} -o {output.png} --averageTypeSummaryPlot median
        """

# rule Plot_deepTools_eRNA_heatmap:
