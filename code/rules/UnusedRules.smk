localrules: DownloadFasterqDump_SE

rule DownloadFasterqDump_SE:
    """
    FasterqDump for single end read datasets. can't be parallelized on midway
    compute nodes. is painfully slow relative to aspera. Using aspera links are
    better...
    This rule is probably not needed because aspera from ENA (which is mostly snychronized with SRA) is faster and more reliable. There is a different rule for using aspera to download fastq files
    """
    output:
        "FastqSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz"
    log:
        "logs/DownloadFasterqDump_SE/{Phenotype}/{IndID}.{Rep}.log"
    params:
        SRA_Run = GetDownloadLinkFuncs('SRA_Run'),
    shadow: "shallow"
    threads: 1
    resources:
        tasks_per_node = 7
    shell:
        """
        fasterq-dump -S -v -e {threads} {params.SRA_Run} &> {log}
        for accession in {params.SRA_Run};
        do
            gzip -c $accession.fastq >> {output}
            rm $accession.fastq
        done
        """



##############################

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
        (samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.bed} -b - -c -split | gzip - > {output}) &> {log}
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
        mem_mb = 32000
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
        
rule SlopesQQnorm:
    input: 
        expand("IntronSlopes/slopes/{IndID}.IntronWindows_equalLength.tab.gz", IndID=chRNASeqSamples)
        #expand("IntronSlopes/IntronWindowCounts/{IndID}.IntronWindows_equalLength.bed.gz", IndID=chRNASeqSamples)
    output:
        "QTLs/QTLTools/chRNA.Slopes/OnlyFirstReps.qqnorm.bed.gz",
        'QTLs/QTLTools/chRNA.Slopes/OnlyFirstReps.Slopes.bed.gz',
        'QTLs/QTLTools/chRNA.Slopes/OnlyFirstReps.Intercept.bed.gz'
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
                                'ENSG00000073282.13_chr3_189631577_189737739_+',
                                'ENSG00000135945.10_chr2_99464985_99489816_-',
                                'ENSG00000154447.15_chr4_169136620_169155479_-',
                                'ENSG00000105708.9_chr19_19714487_19732955_-',
                                'ENSG00000153561.13_chr2_86720809_86740926_+',
                                'ENSG00000148459.16_chr10_26709768_26720217_+',
                                'ENSG00000166478.10_chr11_9479546_9494645_+',
                                'ENSG00000111371.16_chr12_46229639_46239678_-',
                                'ENSG00000250312.8_chr4_131505_160911_+',
                                'ENSG00000166478.10_chr11_9479546_9494645_+',
                                'ENSG00000166348.18_chr10_73575675_73625566_-',
                                'ENSG00000078674.17_chr8_17993619_18006262_+',
                                'ENSG00000035499.13_chr5_60647533_60686961_-',
                                'ENSG00000150403.18_chr13_113503588_113510237_+',
                                'ENSG00000204149.12_chr10_49994429_50001995_+',
                                'ENSG00000164236.12_chr5_10638168_10649265_+',
                                'ENSG00000173273.16_chr8_9615677_9679950_+',
                                'ENSG00000107036.12_chr9_5629453_5656582_+',
                                'ENSG00000054118.15_chr1_36259484_36282532_+',
                                'ENSG00000114999.8_chr2_112503181_112520281_+',
                                'ENSG00000155849.15_chr7_37271882_37314849_-',
                                'ENSG00000180370.10_chr3_196740157_196782625_+'
                                'ENSG00000166575.17_chr11_87236684_87295781_+',
                                'ENSG00000142599.19_chr1_8466023_8495062_-',
                                'ENSG00000114126.17_chr3_142101841_142149182_-',
                                'ENSG00000113761.12_chr5_177064601_177079381_+',
                                'ENSG00000135842.17_chr1_184831962_184884632_-',
                                'ENSG00000119487.17_chr9_125585727_125657650_-',
                                'ENSG00000205209.7_chr19_34596594_34675629_-',
                                'ENSG00000121988.18_chr2_135315530_135345549_-',
                                'ENSG00000164182.11_chr5_61099032_61152703_+',
                                'ENSG00000152127.9_chr2_134120291_134254261_+',
                                'ENSG00000170946.15_chr11_31370859_31414810_+',
                                'ENSG00000275778.2_chr12_11047188_11171421_-',
                                'ENSG00000014164.7_chr8_143475585_143507745_-',
                                'ENSG00000111877.17_chr6_118829250_118856370_-']),
        max_missing = '0.1',
        top_introns = '10000',
        #FDR = "0.25"
    log:
        "logs/slopes/OnlyFirstReps.log"
    shell:
        """
        python scripts/PrepareSlopesQQnorm.py --windowStyle IntronWindows_equalLength --skip_samples {params.skip_samples} --skip_introns {params.skip_introns} --max_missing {params.max_missing} --top_introns {params.top_introns} &> {log}
        """
        
##############################

        
rule CollectSpliceq:
    input:
        expand("QTLs/QTLTools/{Phenotype}.IER/OnlyFirstReps.ratios.bed.gz",
        Phenotype = ['chRNA', 'polyA', 'MetabolicLabelled.30min', 'MetabolicLabelled.60min']),
        expand("QTLs/QTLTools/{Phenotype}.IER/OnlyFirstReps.qqnorm.bed.gz",
        Phenotype = ['chRNA', 'polyA', 'MetabolicLabelled.30min', 'MetabolicLabelled.60min']),
        "QTLs/QTLTools/polyA.IER.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"