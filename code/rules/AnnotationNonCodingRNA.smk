#######################################
#    Prepare input for HMM and run    #
#######################################

rule GetAllGenesForOverlap:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
        #"ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
    output:
        "NonCodingRNA/annotation/allGenes.bed.gz"
    shell:
        """
        awk -F'\\t' '$3=="gene" {{print $1"\\t"$4"\\t"$5"\\t"$7"\\t"$9}}' {input} | awk -F'"' '{{print $1"\\t"$6"\\t"$4}}' | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4}}' | bedtools sort -i - | gzip - > {output}
        """

rule MakeWindows:
    input:
        "../data/Chrome.sizes"
    output:
        plus = "NonCodingRNA/bed/{chrom}.windows.plus.bed.gz",
        minus = "NonCodingRNA/bed/{chrom}.windows.minus.bed.gz",
    params:
        win_size = 50,
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
        HMM_pass = '_annotation|_merged'
    resources:
        mem_mb = 12000
    log:
        "logs/NonCodingRNA/NonCodingRNA/filter_bam.{IndID}.log"
    conda:
        "../envs/py_tools.yml"
    params:
        max_junc_length = 20000
    shell:
        """
        (samtools view -h -F256 {input} | perl -lane 'print if (( $F[0] =~ /^@/ ) || ( abs($F[8]) <= {params.max_junc_length}  ))' | samtools view -bh > {output}) &> {log}
        """
        
        
rule IndexBam:
    input:
        "NonCodingRNA/bam/{IndID}.filtered.bam"
    output:
        temp("NonCodingRNA/bam/{IndID}.filtered.bam.bai")
    log:
        "logs/NonCodingRNA/index.{IndID}_filtered.bam.log"
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        HMM_pass = '_annotation|_merged'
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
        HMM_pass = '_annotation|_merged'
    resources:
        mem_mb = 58000
    log:
        "logs/NonCodingRNA/NonCodingRNA/split_bam.{IndID}.{chrom}.log"
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
        "logs/NonCodingRNA/NonCodingRNA/counts_chRNA.{IndID}.{chrom}.{strand}.log"
    shell:
        """
        (samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.chrom_bed} -b - -c | gzip - > {output}) &> {log}
        """

rule MakeInputForHMM:
    input:
        expand("NonCodingRNA/counts/chRNA/{{chrom}}.{IndID}.{{strand}}.bed.gz",
               IndID = [x for x in chRNASeqSamples if x != "NA18855"]),
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz",
        "NonCodingRNA/tables/{chrom}.{strand}.counts.all_samples.tab.gz",
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    log:
        "logs/NonCodingRNA/NonCodingRNA/make_input_for_hmm.{chrom}.{strand}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        python scripts/MakeInputForHMM.py --chrom {wildcards.chrom} --counts_dir NonCodingRNA/counts/ --strand {wildcards.strand} --output NonCodingRNA/tables/{wildcards.chrom}.{wildcards.strand} &> {log}
        """

rule RunHMM:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.counts.tab.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.predicted_states.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    log:
        "logs/NonCodingRNA/NonCodingRNA/run_hmm.{chrom}.{strand}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        Rscript scripts/runHMM.R {wildcards.chrom} {wildcards.strand} NonCodingRNA/tables/ &> {log};
        gzip NonCodingRNA/tables/{wildcards.chrom}.{wildcards.strand}.predicted_states.tab
        """



########################################
#          Process HMM output          #
########################################

rule ProcessHMMOutput:
    input:
        "NonCodingRNA/tables/{chrom}.plus.counts.tab.gz",
        "NonCodingRNA/tables/{chrom}.minus.counts.tab.gz",
        "NonCodingRNA/tables/{chrom}.{strand}.predicted_states.tab.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.low_expression.bed.gz",
        "NonCodingRNA/tables/{chrom}.{strand}.expression.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 42000
    log:
        "logs/NonCodingRNA/procesHMM.{chrom}.{strand}.log"
    shell:
        """
        python scripts/process_hmm_output.py --chrom {wildcards.chrom} --strand {wildcards.strand} &> {log}
        """
def GetStrandParam(wildcards):
    if wildcards.strand == 'plus':
        return "+"
    else:
        return "-"
       
rule GetTranscriptsFromHMM:
    input:
        counts = "NonCodingRNA/tables/{chrom}.{strand}.counts.all_samples.tab.gz",
        low_expression = "NonCodingRNA/tables/{chrom}.{strand}.low_expression.bed.gz",
        expression = "NonCodingRNA/tables/{chrom}.{strand}.expression.bed.gz",
        genes = "NonCodingRNA/annotation/allGenes.bed.gz",
        chrom_sizes = "../data/Chrome.sizes"
    output:
        sample_counts = "NonCodingRNA/tables/{chrom}.{strand}.sample.counts.bed.gz",
        expression_all = "NonCodingRNA/tables/{chrom}.{strand}.expression.all.bed.gz",
        low_expression_all = "NonCodingRNA/tables/{chrom}.{strand}.low_expression.all.bed.gz",
        pc1 = "NonCodingRNA/tables/{chrom}.{strand}.pc1.bed",
        pc = "NonCodingRNA/tables/{chrom}.{strand}.pc.bed.gz",
        expression_merged = "NonCodingRNA/tables/{chrom}.{strand}.expression.merged.bed.gz",
        low_expression_merged = "NonCodingRNA/tables/{chrom}.{strand}.low_expression.merged.bed.gz",
        expression_filtered = "NonCodingRNA/tables/{chrom}.{strand}.expression.filtered.bed.gz",
        low_expression_filtered = "NonCodingRNA/tables/{chrom}.{strand}.low_expression.filtered.bed.gz",
        filtered_all = "NonCodingRNA/tables/{chrom}.{strand}.filtered.bed",
        overlap_counts = "NonCodingRNA/tables/{chrom}.{strand}.overlap_counts.bed.gz",
        overlap_expression = "NonCodingRNA/tables/{chrom}.{strand}.overlap_expression.bed.gz",
        overlap_low_expression = "NonCodingRNA/tables/{chrom}.{strand}.overlap_low_expression.bed.gz",
        ncRNA_overlaps = "NonCodingRNA/tables/{chrom}.{strand}.ncRNA_overlaps.bed.gz",
        ncRNA_candidates = "NonCodingRNA/tables/{chrom}.{strand}.ncRNA_candidates.bed.gz",
        hmm_unsorted = "NonCodingRNA/tables/{chrom}.{strand}.hmm.bed",
        hmm_annotation = "NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 58000
    params:
        init_counts = ','.join(['0']*86),
        merge_c = ','.join([str(x) for x in range(4,93)]),
        merge_o = ','.join(['sum']*86),
        merge_hash = ','.join(['$' + str(x) for x in range(5,93)]),
        strand = GetStrandParam
    log:
        "logs/NonCodingRNA/GetTranscriptHMM.{chrom}.{strand}.log"
    shell:
        """
        (zcat {input.counts} | tail -n+2 | awk -F'_' '{{print $1"\\t"$2"\\t"$3"\\t.\\t.\\t{params.strand}\\t"$4}}' - | awk -F'\\t{wildcards.strand}\\t' '{{print $1"\\t"$2}}' - | gzip - > {output.sample_counts}) &> {log};
        (bedtools intersect -s -a {output.sample_counts} -b {input.expression} -u | awk '$5=2' FS='\\t' OFS='\\t' - | gzip - > {output.expression_all}) &>> {log};
        (bedtools intersect -s -a {output.sample_counts} -b {input.low_expression} -u | awk '$5=1' FS='\\t' OFS='\\t' - | gzip - > {output.low_expression_all}) &>> {log};
        (zcat {input.genes} | grep protein_coding | grep {wildcards.chrom} - | awk '$6=="{params.strand}" {{print $0, {params.init_counts} }}' OFS='\\t' -  > {output.pc1}) &>> {log};
        (zcat {input.genes} | grep snoRNA | grep {wildcards.chrom} - | awk '$6=="{params.strand}" {{print $0, {params.init_counts} }}' OFS='\\t' - >> {output.pc1}) &>> {log};
        (zcat {output.expression_all} | bedtools merge -i - -s -c {params.merge_c} -o first,collapse,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.expression_merged}) &>> {log};
        (zcat {output.low_expression_all} | bedtools merge -i - -s -c {params.merge_c} -o first,collapse,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.low_expression_merged}) &>> {log};
        (bedtools intersect -s -a {output.sample_counts} -b {output.pc1} -u | bedtools merge -i - -s -c {params.merge_c} -o first,first,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.overlap_counts}) &>> {log};
        (bedtools intersect -s -a {output.expression_merged} -b {output.pc1} -u | bedtools merge -i - -s -c {params.merge_c} -o first,first,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.overlap_expression}) &>> {log};
        (bedtools intersect -s -a {output.low_expression_merged} -b {output.pc1} -u | bedtools merge -i - -s -c {params.merge_c} -o first,first,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.overlap_low_expression}) &>> {log};
        (zcat {output.overlap_counts} >> {output.pc1}) &>> {log};
        #(zcat {output.overlap_expression} >> {output.pc1}) &>> {log};
        (zcat {output.overlap_low_expression} >> {output.pc1}) &>> {log};
        (bedtools sort -i {output.pc1} | bedtools merge -s -i - -c {params.merge_c} -o first,first,first,{params.merge_o} -d 1000 | awk '{{print sep='\\t' $1, $2, $3, "pc_" NR, {params.merge_hash} }}' FS='\\t' OFS='\\t' | gzip - > {output.pc}) &>> {log};
        (bedtools slop -i {output.pc} -b 51 -g {input.chrom_sizes} | bedtools subtract -a {output.low_expression_merged} -b - | gzip - > {output.low_expression_filtered}) &>> {log};
        (bedtools subtract -s -a {output.expression_merged} -b {output.pc} | gzip - > {output.expression_filtered}) &>> {log};
        (zcat {output.low_expression_filtered} > {output.filtered_all}) &>> {log};
        (zcat {output.expression_filtered} >> {output.filtered_all}) &>> {log};
        (bedtools sort -i {output.filtered_all} | bedtools merge -s -i - -c {params.merge_c} -o first,collapse,first,{params.merge_o} -d 1000 | bedtools sort -i - | awk '{{print sep='\\t' $1, $2, $3, "ncRNA_" NR, {params.merge_hash} }}' FS='\\t' OFS='\\t' | gzip - > {output.ncRNA_overlaps}) &>> {log};
        (bedtools intersect -s -a {output.ncRNA_overlaps} -b {output.expression_filtered} -u | gzip - > {output.ncRNA_candidates}) &>> {log};
        (zcat {output.ncRNA_candidates} > {output.hmm_unsorted}) &>> {log};
        (zcat {output.pc} >> {output.hmm_unsorted}) &>> {log};
        (bedtools sort -i {output.hmm_unsorted} | gzip - > {output.hmm_annotation}) &>> {log};
        """
        
        
        
        
        
        
        
        
        

