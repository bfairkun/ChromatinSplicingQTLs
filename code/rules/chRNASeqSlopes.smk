rule GetGenomeElements:
    input: 
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        "IntronSlopes/Annotation/genome_exons.bed"
    shell:
        """
        awk -F'\t' '$3=="exon"' {input.gtf} | awk -F'\t' '{{print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7}}' >> {output}
        """

#grep -v protein_coding {input.gtf} | tail -n+6 | awk -F'\t' '{{print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7}}' > {output} && grep protein_coding {input.gtf} | awk -F'\t' '$3=="exon"' | awk -F'\t' '{{print $1"\t"$4"\t"$5"\tprotein_coding_"$3"\t"$6"\t"$7}}' >> {output}
    
    
rule GetIntronsInExpressedGenes:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        bed = "../data/Introns.GencodeV34.hg38.UCSC.bed.gz",
        genelist = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        "IntronSlopes/Annotation/GencodeHg38_all_introns.expressedHostGenes.bed.gz"
    log:
        "logs/chRNA-seq_FilterForExpressedGenes.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/chRNA-seq_FilterForExpressedGenes.R {input} {output} &> {log}
        """

rule GetUniqIntrons:
    input:
        bed = "IntronSlopes/Annotation/GencodeHg38_all_introns.expressedHostGenes.bed.gz",
        elements = "IntronSlopes/Annotation/genome_exons.bed",
        faidx = "../data/Chrome.sizes"
    output:
        "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed"
    params:
        max_overlap = 0.01,
    shell:
        """
        zcat {input.bed} | awk -F'\\t' -v OFS='\\t' '$1 !~ "_" {{ print $1, $2, $3,$7"_"$1"_"$2"_"$3"_"$6,".", $6 }}' | sort | uniq | bedtools sort -i - -faidx {input.faidx} | bedtools intersect -s -v -f {params.max_overlap} -r -a - -b {input.elements} > {output}
        """        
        
rule MakeWindowsForIntrons:
    """
    For calculating downward slope over introns. Like in Jonkers & Lis. Split
    introns into equal number of bins
    """
    input:
        bed = "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = "../data/Chrome.sizes"
    params:
        MinIntronLength = 500,
    output:
        "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed.IntronWindows.bed"
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
        bed = "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = "../data/Chrome.sizes"
    params:
        WinLen = 200,
        MinIntronLength = 500,
    output:
        "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed.IntronWindows_equalLength.bed"
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
        bed = "IntronSlopes/Annotation/GencodeHg38_all_introns.corrected.uniq.bed.{windowStyle}.bed",
        faidx = "../data/Chrome.sizes",
        bam = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam'
    log:
        "logs/CountReadsInIntronWindows/{IndID}.{windowStyle}.log"
    output:
        bed = "IntronSlopes/IntronWindowCounts/{IndID}.{windowStyle}.bed.gz"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        windowStyle="IntronWindows_equalLength|IntronWindows"
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
        bed = "IntronSlopes/IntronWindowCounts/{IndID}.{windowStyle}.bed.gz"
        #bed = "IntronSlopes/IntronWindowCounts/{IndID}.IntronWindows_equalLength.bed.gz"
    params:
        WinLen = "200",
        minIntronCounts = "1000",
        minCoverageCounts = "20",
        minCoverage = "0.9",
        minIntronLen = "1000"
    output:
        "IntronSlopes/slopes/{IndID}.{windowStyle}.tab.gz",
        #"IntronSlopes/slopes/{IndID}.{windowStyle}.glm_nb.tab.gz"
    resources:
        mem_mb = 16000
    conda:
        "../envs/r_slopes.yml"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        windowStyle = "IntronWindows_equalLength|IntronWindows"
    log:
        "logs/slopes/{IndID}.{windowStyle}.slope.log"
    shell:
        """
        (Rscript scripts/GetSlopes.R {input.bed} {params.minIntronCounts} {params.minCoverageCounts} {params.minCoverage} {params.WinLen} {params.minIntronLen}) &> {log}
        """
        
rule SlopesPreparePhenotypes:
    input: 
        expand("IntronSlopes/IntronWindowCounts/{IndID}.IntronWindows_equalLength.bed.gz", IndID=chRNASeqSamples)
    output:
        "QTLs/QTLTools/chRNA.Slopes/OnlyFirstReps.qqnorm.bed.gz"
    conda:
        "../envs/py_tools.yml"
    params:
        skip_samples = 'NA18855',
        skip_introns = ':'.join(['ENSG00000153944.11_chr17_57596950_57615969_+',
                                'ENSG00000187231.14_chr2_179191866_179264498_-',
                                'ENSG00000170832.13_chr17_60301704_60345480_-',
                                'ENSG00000062725.10_chr17_60500487_60525793_-',
                                'ENSG00000204217.15_chr2_202377550_202464808_+',
                                'ENSG00000182670.13_chr21_37172744_37182773_+',
                                'ENSG00000131374.14_chr3_17428519_17508473_-',
                                'ENSG00000188033.10_chr19_12583556_12609157_-',
                                'ENSG00000152492.15_chr3_191329723_191357087_+',
                                'ENSG00000198874.13_chr7_67195337_67238307_+',
                                'ENSG00000170776.22_chr15_85485753_85521427_+',
                                'ENSG00000181722.16_chr3_114900342_114974365_-',
                                'ENSG00000088930.8_chr20_21303473_21326278_+',
                                'ENSG00000136813.14_chr9_111397153_111408570_-',
                                'ENSG00000073282.13_chr3_189631577_189737739_+']),
        max_missing = '0.1',
        top_introns = '10000',
        #FDR = "0.25"
    log:
        "logs/slopes/OnlyFirstReps.log"
    shell:
        """
        python scripts/PreparePhenotypeTablesSlopes.py --windowStyle IntronWindows_equalLength --skip_samples {params.skip_samples} --skip_introns {params.skip_introns} --max_missing {params.max_missing} --top_introns {params.top_introns} &> {log}
        """
        





