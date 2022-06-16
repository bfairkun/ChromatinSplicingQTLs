

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
        
rule FilterLongJunc:
    input:
        "Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam",
    output:
        temp("NonCodingRNA/bam/{IndID}.filtered.bam")
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
    resources:
        mem_mb = 12000
    log:
        "logs/NonCodingRNA/filter_bam.{IndID}.log"
    conda:
        "../envs/py_tools.yml"
    params:
        max_junc_length = 10000
    shell:
        """
        python scripts/filter_long_junc.py --input {input} --output {output} --max_junc_length {params.max_junc_length} &> {log}
        """
        
rule IndexBam:
    input:
        "NonCodingRNA/bam/{IndID}.filtered.bam"
    output:
        "NonCodingRNA/bam/{IndID}.filtered.bam.bai"
    log:
        "logs/index.{IndID}_filtered.bam.log"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
    resources:
        mem_mb = 12000
    shell:
        """
        (samtools index -b {input}) &> {log}
        """

rule SplitBamByChromosome:
    input:
        bam = "NonCodingRNA/bam/{IndID}.filtered.bam",
        bai = "NonCodingRNA/bam/{IndID}.filtered.bam.bai"
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
        (samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.chrom_bed} -b - -c | gzip - > {output}) &> {log}
        """
        #(samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.chrom_bed} -b - -c -split | gzip - > {output}) &> {log}

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
        #expand("NonCodingRNA/counts/ProSeq/{{chrom}}.{IndID}.{{strand}}.bed.gz",
        #       IndID = proseq_samples)
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
        python scripts/MakeInputForHMM.py --chrom {wildcards.chrom} --counts_dir NonCodingRNA/counts/ --strand {wildcards.strand} --output {output} &> {log}
        """

rule RunHMM:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.predicted_{states}states.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        states = '2|3'
    log:
        "logs/NonCodingRNA/run_hmm.{chrom}.{strand}.{states}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        Rscript scripts/runHMM.R {wildcards.chrom} {wildcards.strand} {wildcards.states} &> {log};
        gzip NonCodingRNA/tables/{wildcards.chrom}.{wildcards.strand}.predicted_{wildcards.states}states.tab
        """

def getStrandString(wildcards):
    if wildcards.strand == 'plus':
        return '+'
    else:
        return '-'
        
        
rule GetAllGenesForOverlap:
    input:
        "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
    output:
        "NonCodingRNA/annotation/allGenes.bed.gz"
    shell:
        """
        awk -F'\\t' '{{print "chr"$1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$4}}' {input} | bedtools sort -i - | gzip -  > {output}
        """
    
def GetMergedDistance(wildcards):
    if wildcards.nstates == '3':
        return 1000
    elif wildcards.nstates == '2':
        return 2000

rule ProcessAnnotationsHMM:
    input:
        pred = "NonCodingRNA/tables/{chrom}.{strand}.predicted_{nstates}states.tab.gz",
        genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz"
    output:
        TU_bed = "NonCodingRNA/tables/{chrom}.{strand}.TU_{nstates}states.tab.gz",
        merged_bed = "NonCodingRNA/tables/{chrom}.{strand}.merged_{nstates}states.tab.gz",
        ncRNA_bed = "NonCodingRNA/tables/{chrom}.{strand}.{nstates}states.ncRNA.bed.gz",
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        nstates = '2|3'
    params:
        strand = getStrandString,
        merge_distance = 1000,#GetMergedDistance,
        max_overlap = 0.1
    log:
        "logs/NonCodingRNA/merge_hmm.{chrom}.{strand}.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/processHMM.py --input {input.pred} --output {output.TU_bed} --nstates {wildcards.nstates} &> {log};
        (bedtools merge -i {output.TU_bed} -d {params.merge_distance} | awk -F'\\t' '{{print $0"\\t"$1"_"$2"_"$3"_{params.strand}\\t"$1"_"$2"_"$3"_{params.strand}\\t{params.strand}"}}' - | sort -u | gzip - > {output.merged_bed}) &>> {log};
        (bedtools subtract -A -s -a {output.merged_bed} -b {input.genes_bed} -f {params.max_overlap} | sort -u | gzip - > {output.ncRNA_bed}) &>> {log}
        """
        
rule MergeNonCodingRNA:
    input:
        expand("NonCodingRNA/tables/{chrom}.{strand}.{{nstates}}states.ncRNA.bed.gz",
                chrom = chrom_list, strand = ['plus', 'minus'])
    output:
        tmp_all_bed = temp("NonCodingRNA/annotation/ncRNA.{nstates}states.bed.gz"),
        tmp_filtered_bed = temp("NonCodingRNA/annotation/ncRNA_filtered.{nstates}states.bed.gz"),
        tmp_hmm_bed = temp("NonCodingRNA/annotation/allHMM.{nstates}states.bed.gz"),
        all_bed = "NonCodingRNA/annotation/ncRNA.{nstates}states.sorted.bed.gz",
        filtered_bed = "NonCodingRNA/annotation/ncRNA_filtered.{nstates}states.sorted.bed.gz",
        hmm_bed = "NonCodingRNA/annotation/allHMM.{nstates}states.sorted.bed.gz"
    params:
        length = 1000,
    wildcard_constraints:
        nstates = '2|3'
    log:
        "logs/NonCodingRNA/hmm_annotation.{nstates}states.log"
    resources:
        mem_mb = 58000
    shell:
        """
        python scripts/merge_HMM_annotation.py --length {params.length} --nstates {wildcards.nstates} &> {log};
        bedtools sort -i {output.tmp_all_bed} | gzip - > {output.all_bed};
        bedtools sort -i {output.tmp_filtered_bed} | gzip - > {output.filtered_bed};
        bedtools sort -i {output.tmp_hmm_bed} | gzip - > {output.hmm_bed};
        """
   
rule ProteinCodingOverlap:
    input:
        gene_list = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        hmm_bed = "NonCodingRNA/annotation/allHMM.{nstates}states.sorted.bed.gz"
    output:
        protein_coding = "NonCodingRNA/annotation/AnnotatedProteinCoding_overlap.{nstates}states.bed.gz",
        de_novo = "NonCodingRNA/annotation/DeNovoProteinCoding_overlap.{nstates}states.bed.gz"
    params:
        nstates = '2|3'
    log:
        "logs/NonCodingRNA/hmm_protein_coding.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        bedtools intersect -s -a {input.gene_list} -b {input.hmm_bed} -f 0.5 -wa | sort -u | gzip - > {output.protein_coding};
        bedtools intersect -s -b {input.gene_list} -a {input.hmm_bed} -f 0.5 -wa | sort -u | gzip - > {output.de_novo}
        """


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
        bedgraph = temp("NonCodingRNA/bigwig/chRNA.{IndID}.{strand}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA/bigwig/chRNA.{IndID}.{strand}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA/bigwig/chRNA.{IndID}.{strand}.log1p.bw"
    wildcard_constraints:
        IndID = "NA18486|NA19201|NA19210",
        strand = 'plus|minus'
    log:
        "logs/NonCodingRNA/chRNA.log1p.{IndID}.{strand}.log"
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
        bedgraph = temp("NonCodingRNA/bigwig/polyA.{IndID}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA/bigwig/polyA.{IndID}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA/bigwig/polyA.{IndID}.log1p.bw"
    wildcard_constraints:
        IndID = "NA18486|NA19201|NA19210",
    log:
        "logs/NonCodingRNA/polyA.log1p.{IndID}.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bigWigToBedGraph {input.bigwig} {output.bedgraph}) &> {log};
        python scripts/LogTransformBigWig.py --input {output.bedgraph} --output {output.bedgraph_log} &>> {log};
        (bedGraphToBigWig {output.bedgraph_log} {input.chrom_sizes} {output.bigwig_log}) &>> {log};
        """
    
rule SplitAnnotationByStrand:
    input:
        "NonCodingRNA/annotation/ncRNA_filtered.2states.sorted.bed.gz"
    output:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
    log:
        "logs/NonCodingRNA/split_annotation.log"
    shell:
        """
        (zcat NonCodingRNA/annotation/ncRNA_filtered.2states.sorted.bed.gz | awk -F'\\t' '$6=="+"' -  > {output.plus_bed}) &> {log};
        (zcat NonCodingRNA/annotation/ncRNA_filtered.2states.sorted.bed.gz | awk -F'\\t' '$6=="-"' -  > {output.minus_bed}) &>> {log};
        """
        
        
rule GetTSSAnnotations:
    input:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz",
    output:
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA_filtered.2states.TSS.bed.gz",
        genes_plus = temp("NonCodingRNA/annotation/allGenes.plus.bed"),
        genes_minus = temp("NonCodingRNA/annotation/allGenes.minus.bed"),
        genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    log:
        "logs/TSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 12000
    shell:
        """
        awk -F'\\t' '{{print $1"\\t"$2-{params.window}"\\t"$2+{params.window}"\\t"$4"\\t"$5"\\t"$6}}' {input.plus_bed} > NonCodingRNA/annotation/ncRNA_filtered.2states.TSS.bed;
        awk -F'\\t' '{{print $1"\\t"$3-{params.window}"\\t"$3+{params.window}"\\t"$4"\\t"$5"\\t"$6}}' {input.minus_bed} >> NonCodingRNA/annotation/ncRNA_filtered.2states.TSS.bed;
        gzip NonCodingRNA/annotation/ncRNA_filtered.2states.TSS.bed;
        zcat {input.genes_bed} | awk -F'\\t' '$6=="+"' -  > {output.genes_plus};
        zcat {input.genes_bed} | awk -F'\\t' '$6=="-"' -  > {output.genes_minus};
        awk -F'\\t' '{{print $1"\\t"$2-{params.window}"\\t"$2+{params.window}"\\t"$4"\\t"$5"\\t"$6}}' {output.genes_plus} > NonCodingRNA/annotation/allGenes.TSS.bed;
        awk -F'\\t' '{{print $1"\\t"$3-{params.window}"\\t"$3+{params.window}"\\t"$4"\\t"$5"\\t"$6}}' {output.genes_minus} >> NonCodingRNA/annotation/allGenes.TSS.bed;
        gzip NonCodingRNA/annotation/allGenes.TSS.bed;
        """
        

rule deepTools_chRNA_plotHeatmap:
    input:
        chRNA_plus = "NonCodingRNA/bigwig/chRNA.{IndID}.plus.log1p.bw",
        chRNA_minus = "NonCodingRNA/bigwig/chRNA.{IndID}.minus.log1p.bw",
        polyA = "NonCodingRNA/bigwig/polyA.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA/deeptools/matrix/{IndID}.RNA{strictness}.matrix.gz",
        png = "NonCodingRNA/deeptools/plots/{IndID}.RNA{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA/plots.{IndID}{strictness}.log"
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
        ncPlus = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        ncMinus = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        genes = "NonCodingRNA/annotation/allGenes.bed.gz",
        chrom = "../data/Chrome.sizes"
    output:
        genes1 = "NonCodingRNA/deeptools/bed/allGenes_expanded1.bed.gz",
        genes2 = "NonCodingRNA/deeptools/bed/allGenes_expanded2.bed.gz",
        ncPlus1 = "NonCodingRNA/deeptools/bed/ncRNA.plus.strict1.bed",
        ncMinus1 = "NonCodingRNA/deeptools/bed/ncRNA.minus.strict1.bed",
        ncPlus2 = "NonCodingRNA/deeptools/bed/ncRNA.plus.strict2.bed",
        ncMinus2 = "NonCodingRNA/deeptools/bed/ncRNA.minus.strict2.bed",
        ncPlus3 = "NonCodingRNA/deeptools/bed/ncRNA.plus.strict3.bed",
        ncMinus3 = "NonCodingRNA/deeptools/bed/ncRNA.minus.strict3.bed",
    log:
        "logs/NonCodingRNA/strict_filtering.log"
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
        bedgraph = temp("NonCodingRNA/bigwig/{Phenotype}.{IndID}.bedgraph"),
        bedgraph_log = temp("NonCodingRNA/bigwig/{Phenotype}.{IndID}.log1p.bedgraph"),
        bigwig_log = "NonCodingRNA/bigwig/{Phenotype}.{IndID}.log1p.bw"
    wildcard_constraints:
        IndID = "NA19201|NA19210",
        Phenotype = "H3K27AC|H3K4ME1|H3K4ME3|H3K36ME3|ProCap"
    log:
        "logs/NonCodingRNA/{Phenotype}.log1p.{IndID}.log"
        

rule deepTools_histone_plotHeatmap:
    input:
        procap = "NonCodingRNA/bigwig/ProCap.{IndID}.log1p.bw",
        h3k4me1 = "NonCodingRNA/bigwig/H3K4ME1.{IndID}.log1p.bw",
        h3k4me3 = "NonCodingRNA/bigwig/H3K4ME3.{IndID}.log1p.bw",
        h3k27ac = "NonCodingRNA/bigwig/H3K27AC.{IndID}.log1p.bw",
        h3k36me3 = "NonCodingRNA/bigwig/H3K36ME3.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA/deeptools/matrix/{IndID}.histone{strictness}.matrix.gz",
        png = "NonCodingRNA/deeptools/plots/{IndID}.histone{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA/histone_plots.{IndID}{strictness}.log"
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
        chRNA_plus = "NonCodingRNA/bigwig/chRNA.{IndID}.plus.log1p.bw",
        chRNA_minus = "NonCodingRNA/bigwig/chRNA.{IndID}.minus.log1p.bw",
        polyA = "NonCodingRNA/bigwig/polyA.{IndID}.log1p.bw",
        procap = "NonCodingRNA/bigwig/ProCap.{IndID}.log1p.bw",
        h3k4me1 = "NonCodingRNA/bigwig/H3K4ME1.{IndID}.log1p.bw",
        h3k4me3 = "NonCodingRNA/bigwig/H3K4ME3.{IndID}.log1p.bw",
        h3k27ac = "NonCodingRNA/bigwig/H3K27AC.{IndID}.log1p.bw",
        h3k36me3 = "NonCodingRNA/bigwig/H3K36ME3.{IndID}.log1p.bw",
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus{strictness}.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus{strictness}.bed",
    output:
        mat = "NonCodingRNA/deeptools/matrix/{IndID}.combined{strictness}.matrix.gz",
        png = "NonCodingRNA/deeptools/plots/{IndID}.combined{strictness}.plot.png"
    conda:
        "../envs/deeptools.yml"
    log:
        "logs/NonCodingRNA/combined_plots.{IndID}{strictness}.log"
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
        
 
