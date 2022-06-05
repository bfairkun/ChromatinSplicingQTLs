

rule DownloadProSeq:
    output:
        temp(expand("ProSeq/bedgraph/{IndID}.{strand}.hg19.bedgraph.gz",
                IndID = proseq_samples,
                strand = ['plus', 'minus'])
        )
    shell:
        """
        wget -O ProSeq/bedgraph/GM18505.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004714/suppl/GSM3004714_48651_R1_cap.clip.TTAGGC.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18505.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004714/suppl/GSM3004714_48651_R1_cap.clip.TTAGGC.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19239.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004723/suppl/GSM3004723_48636_R1_cap.clip.TGACCA.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19239.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004723/suppl/GSM3004723_48636_R1_cap.clip.TGACCA.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19238.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004722/suppl/GSM3004722_48636_R1_cap.clip.TTAGGC.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19238.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004722/suppl/GSM3004722_48636_R1_cap.clip.TTAGGC.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19222.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004721/suppl/GSM3004721_48636_R1_cap.clip.CGATGT.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19222.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004721/suppl/GSM3004721_48636_R1_cap.clip.CGATGT.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19193.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004720/suppl/GSM3004720_48636_R1_cap.clip.ATCACG.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19193.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004720/suppl/GSM3004720_48636_R1_cap.clip.ATCACG.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19131.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004719/suppl/GSM3004719_48636_R1_cap.clip.ACTTGA.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19131.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004719/suppl/GSM3004719_48636_R1_cap.clip.ACTTGA.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19099.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004718/suppl/GSM3004718_51634_R1_cap.clip.GCCAAT.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM19099.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004718/suppl/GSM3004718_51634_R1_cap.clip.GCCAAT.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18522.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004717/suppl/GSM3004717_48651_R1_cap.clip.GCCAAT.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18522.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004717/suppl/GSM3004717_48651_R1_cap.clip.GCCAAT.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18520.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004716/suppl/GSM3004716_48651_R1_cap.clip.ACAGTG.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18520.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004716/suppl/GSM3004716_48651_R1_cap.clip.ACAGTG.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18517.minus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004715/suppl/GSM3004715_48651_R1_cap.clip.TGACCA.hs37d5.bwa.uniqueUMI.mn.bedgraph.gz
        wget -O ProSeq/bedgraph/GM18517.plus.hg19.bedgraph.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3004nnn/GSM3004715/suppl/GSM3004715_48651_R1_cap.clip.TGACCA.hs37d5.bwa.uniqueUMI.pl.bedgraph.gz
        """

rule LiftOverProSeq:
    input:
        bedgraph = "ProSeq/bedgraph/{IndID}.{strand}.hg19.bedgraph.gz",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        "ProSeq/bedgraph/{IndID}.{strand}.hg38.bedgraph.gz",
        temp("ProSeq/unmapped.{IndID}.{strand}")
    wildcard_constraints:
        IndID = "|".join(proseq_samples),
        strand = 'minus|plus'
    resources:
        mem_mb = 48000
    log:
        "logs/liftover.{IndID}.{strand}.log"
    shell:
        """
        liftOver {input.bedgraph} {input.chain} ProSeq/bedgraph/{wildcards.IndID}.{wildcards.strand}.hg38.bedgraph ProSeq/unmapped &> {log};
        gzip ProSeq/bedgraph/{wildcards.IndID}.{wildcards.strand}.hg38.bedgraph
        """


#rule DefineGenomeWindows:
#    output:
#        expand("NonCodingRNA/bed/{chrom}.{strand}.{segment}.bed.gz",
#               chrom=chrom_list,
#               strand = ['plus', 'minus'],
#               segment = GenomeWindowsByChrom)
#    log:
#        "logs/makegenomewindows.log"
#    resources:
#        mem_mb = 24000
#    params:
#        length = GetChromLength
#    shell:
#        """
#        python scripts/MakeGenomeWindows.py --chrom {wildcards.chrom} --length {params.length} &> {log}
#        """

rule MakeWindows:
    input:
        "../data/Chrome.sizes"
    output:
        plus = "NonCodingRNA/bed/{chrom}.windows.plus.bed.gz",
        minus = "NonCodingRNA/bed/{chrom}.windows.minus.bed.gz",
    params:
        win_size = 200
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    shell:
        """
        bedtools makewindows -g {input} -w {params} | awk '$1=="{wildcards.chrom}" {{print $1"\\t"$2"\\t"$3"\\t.\\t.\\t+"}}' > {output.plus};
        bedtools makewindows -g {input} -w {params} | awk '$1=="{wildcards.chrom}" {{print $1"\\t"$2"\\t"$3"\\t.\\t.\\t-"}}' > {output.minus}
        """

rule SplitBamByChromosome:
    input:
        bam = "Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam",
    output:
        temp("NonCodingRNA/bam/{IndID}.{chrom}.bam")
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        chrom = "|".join(chrom_list),
    resources:
        mem_mb = 58000
    log:
        "logs/NonCodingRNA/split_bam.{IndID}.{chrom}.log"
    shell:
        """
        (samtools view -bh -F 256 -f 64 {input.bam} {wildcards.chrom} > {output}) &> {log}
        """

rule WindowCountsChRNA:
    input:
        faidx = "../data/Chrome.sizes",
        bam = "NonCodingRNA/bam/{IndID}.{chrom}.bam",
        chrom_bed = "NonCodingRNA/bed/{chrom}.windows.{strand}.bed.gz"
    output:
        "NonCodingRNA/counts/chRNA/{chrom}.{IndID}.{strand}.bed.gz"
    wildcard_constraints:
        IndID = "|".join([x for x in chRNASeqSamples if x != "NA18855"]),
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 58000
    log:
        "logs/NonCodingRNA/counts_chRNA.{IndID}.{chrom}.{strand}.log"
    shell:
        """
        (samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.chrom_bed} -b - -c -split | gzip - > {output}) &> {log}
        """

rule WindowCountsProSeq:
    input:
        windows_bed = "NonCodingRNA/bed/{chrom}.windows.{strand}.bed.gz",
        proseq_bed = "ProSeq/bedgraph/{IndID}.{strand}.hg38.bedgraph.gz"
    output:
        "NonCodingRNA/counts/ProSeq/{chrom}.{IndID}.{strand}.bed.gz"
    wildcard_constraints:
        IndID = "|".join(proseq_samples),
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 24000
    log:
        "logs/NonCodingRNA/counts_proseq.{IndID}.{chrom}.{strand}.log"
    shell:
       """
       (bedtools intersect -wao -a {input.windows_bed} -b {input.proseq_bed} | awk '{{print $1"\\t"$2"\\t"$3"\\t"$10}}' - | bedtools merge -d -1 -c 4 -o sum -i - | gzip - > {output}) &> {log}
       """

rule MakeInputForHMM:
    input:
        expand("NonCodingRNA/counts/chRNA/{{chrom}}.{IndID}.{{strand}}.bed.gz",
               IndID = [x for x in chRNASeqSamples if x != "NA18855"]),
        expand("NonCodingRNA/counts/ProSeq/{{chrom}}.{IndID}.{{strand}}.bed.gz",
               IndID = proseq_samples)
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    log:
        "logs/NonCodingRNA/make_input_for_hmm.{chrom}.{strand}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        python scripts/MakeInputForHMM.py --chrom {wildcards.chrom} --counts_dir NonCodingRNA/counts/ --strand {wildcards.strand} --output {output} > {log}
        """

rule RunHMM:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.predicted.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    log:
        "logs/NonCodingRNA/run_hmm.{chrom}.{strand}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        Rscript scripts/runHMM.R {wildcards.chrom} {wildcards.strand} > {log};
        gzip NonCodingRNA/tables/{wildcards.chrom}.{wildcards.strand}.predicted.tab
        """

#rule ProcessAnnotationsHMM:

rule deepTools_chRNA_TSS_heatmap:
    input:
        bigwigs = GetBigwigForDeeptToolscheRNA,
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
