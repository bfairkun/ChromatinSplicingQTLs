rule featureCountsExons:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForPhenotype
    output:
        "featureCounts/{Phenotype}/CountsExons.txt"
    params:
        extraParams = GetFeatureCountsParams,
    threads:
        8
    wildcard_constraints:
        Phenotype = "chRNA.Expression.Splicing"
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.log"
    shell:
        """
        featureCounts -p {params.extraParams} -f -t exon -g exon_id -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

rule GetGenomeElements:
    input: 
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        "Misc/genome_elements.bed"
    shell:
        """
        grep -v protein_coding {input.gtf} | tail -n+6 | awk -F'\t' '{{print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7}}' > {output} && grep protein_coding {input.gtf} | awk -F'\t' '$3=="exon"' | awk -F'\t' '{{print $1"\t"$4"\t"$5"\tprotein_coding_"$3"\t"$6"\t"$7}}' >> {output}
        """
    
rule GetLogRPKM:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        CountTable = "featureCounts/chRNA.Expression.Splicing/Counts.txt",
        bed = "../data/Introns.GencodeV34.hg38.UCSC.bed.gz",
    output:
        "featureCounts/chRNA.Expression.Splicing/CountTable.MeanLogRPKM.txt.gz",
        "Misc/GencodeHg38_all_introns.expressedHostGenes.bed.gz"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/chRNA-seq_FilterForExpressedGenes.R
        """

rule GetUniqIntrons:
    input:
        bed = "Misc/GencodeHg38_all_introns.expressedHostGenes.bed.gz",
        elements = "Misc/genome_elements.bed",
        faidx = "ReferenceGenome/Annotations/Chrome.sizes"
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed"
    shell:
        """
        zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '$1 !~ "_" {{ print $1, $2, $3,$7"_"$1"_"$2"_"$3"_"$6,".", $6 }}' | sort | uniq | bedtools sort -i - -faidx {input.faidx} | bedtools intersect -v -a - -b {input.elements} > {output}
        """        
        
rule MakeWindowsForIntrons:
    """
    For calculating downward slope over introns. Like in Jonkers & Lis. Split
    introns into equal number of bins
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = "ReferenceGenome/Annotations/Chrome.sizes"
    params:
        MinIntronLength = config["MinIntronLength"],
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed.IntronWindows.bed"
    shell:
        """
        set +o pipefail;
        bedtools intersect -a {input.bed} -b {input.bed} -c -s -sorted | awk -F '\\t' '$NF==1 && ($3-$2)>={params.MinIntronLength}' | bedtools makewindows -b - -n 100 -i srcwinnum | bedtools sort -i - -faidx {input.faidx} | awk -F'\\t' -v OFS='\\t' '{{ split($4,a,"_") }} a[5]=="+" {{print $1,$2,$3,$4,".",a[5]}} a[5]=="-" {{print $1,$2,$3, a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"101-a[6],".", a[5] }}' > {output}
        """


rule MakeWindowsForIntrons_equalSized:
    """
    For calculating downward slope over introns. Like in Jonkers & Lis. Split
    introns into equal length bins
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = "ReferenceGenome/Annotations/Chrome.sizes"
    params:
        WinLen = config["WinLen"],
        MinIntronLength = config["MinIntronLength"],
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed.IntronWindows_equalLength.bed"
    shell:
        """
        set +o pipefail;
        bedtools intersect -a {input.bed} -b {input.bed} -c -s -sorted | awk -F '\\t' '$NF==1 && ($3-$2)>={params.MinIntronLength}' | bedtools makewindows -b - -w {params.WinLen} -i srcwinnum | bedtools sort -i - -faidx {input.faidx} | awk -F'\\t' -v OFS='\\t' '{{ split($4,a,"_"); print $1,$2,$3,$4,".",a[5] }}'  > {output}
        """


rule CountReadsInIntronWindows:
    """
    Ignore alignments that are secondary or second read pair. Intersect reads
    in anti-stranded fashion (R1 should map to antisense strand with this
    library prep method)
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed.{windowStyle}.bed",
        faidx = "ReferenceGenome/Annotations/Chrome.sizes",
        bam = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam'
    log:
        "logs/CountReadsInIntronWindows/{IndID}.{windowStyle}.log"
    output:
        bed = "IntronWindowCounts/{IndID}.{windowStyle}.bed.gz"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples_df['IndID'].unique()),
        windowStyle="|".join(["IntronWindows_equalLength", "IntronWindows"])
    shell:
        """
        echo {input.bam};
        echo {input.faidx};
        echo {input.bed};
        echo {output};
        set +o pipefail;
        (samtools view -bh -F 256 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.bed} -b - -c -split -F 0.5 | gzip - > {output}) &> {log}
        """
        

rule GetSlopes:
    input:
        bed = "IntronWindowCounts/{IndID}.IntronWindows_equalLength.bed.gz"
    params:
        WinLen = config["WinLen"],
        minIntronCounts = config["minIntronCounts"],
        minCoverageCounts = config["minCoverageCounts"],
        minCoverage = config["minCoverage"]
    output:
        "slopes/{IndID}.tab.gz",
        "slopes/{IndID}_glm.nb.tab.gz"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples_df['IndID'].unique()),
    log:
        "logs/slopes/{IndID}.slope.log"
    shell:
        """
        (Rscript scripts/GetSlopes.R {input.bed} {params.minIntronCounts} {params.minCoverageCounts} {params.minCoverage} {params.WinLen}) &> {log}
        """
