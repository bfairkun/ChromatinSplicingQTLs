#######################################
#      Annotated non-coding RNAs      #
#######################################

def FeatureCountsNonCodingStrandParams(wildcards):
    if wildcards.Phenotype == 'chRNA.Expression':
        return "-s 2"
    else:
        return ""
        
rule GetAdditionalNonCodingRNAFromFeatureCounts:
    input:
        fCRNA = "featureCounts/{Phenotype}/Counts.txt",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf"
    output:
        "featureCounts/{Phenotype}_annotated_ncRNA/Counts.txt",
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"]),
    log:
        "logs/NonCodingRNA/{Phenotype}.get_annotated_ncRNA.log",
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/NonCodingRNA/GetNonCodingRNAFromFeatureCounts.py --phenotype {wildcards.Phenotype} &> {log};
        """

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
        HMM_pass = '_annotation|'
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
        HMM_pass = '_annotation|'
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
        HMM_pass = '_annotation|'
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
        python scripts/NonCodingRNA/MakeInputForHMM.py --chrom {wildcards.chrom} --counts_dir NonCodingRNA/counts/ --strand {wildcards.strand} --output NonCodingRNA/tables/{wildcards.chrom}.{wildcards.strand} &> {log}
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
        Rscript scripts/NonCodingRNA/runHMM.R {wildcards.chrom} {wildcards.strand} NonCodingRNA/tables/ &> {log};
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
        python scripts/NonCodingRNA/process_hmm_output.py --chrom {wildcards.chrom} --strand {wildcards.strand} &> {log}
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
        expression_all = temp("NonCodingRNA/tables/{chrom}.{strand}.expression.all.bed.gz"),
        low_expression_all = temp("NonCodingRNA/tables/{chrom}.{strand}.low_expression.all.bed.gz"),
        pc1 = temp("NonCodingRNA/tables/{chrom}.{strand}.pc1.bed"),
        pc = temp("NonCodingRNA/tables/{chrom}.{strand}.pc.bed.gz"),
        expression_filtered = temp("NonCodingRNA/tables/{chrom}.{strand}.expression.filtered.bed.gz"),
        low_expression_filtered = temp("NonCodingRNA/tables/{chrom}.{strand}.low_expression.filtered.bed.gz"),
        filtered_all = temp("NonCodingRNA/tables/{chrom}.{strand}.filtered.bed"),
        overlap_counts = temp("NonCodingRNA/tables/{chrom}.{strand}.overlap_counts.bed.gz"),
        #overlap_expression = temp("NonCodingRNA/tables/{chrom}.{strand}.overlap_expression.bed.gz"),
        #overlap_low_expression = temp("NonCodingRNA/tables/{chrom}.{strand}.overlap_low_expression.bed.gz"),
        ncRNA_overlaps = temp("NonCodingRNA/tables/{chrom}.{strand}.ncRNA_overlaps.bed.gz"),
        ncRNA_candidates = temp("NonCodingRNA/tables/{chrom}.{strand}.ncRNA_candidates.bed.gz"),
        hmm_unsorted = temp("NonCodingRNA/tables/{chrom}.{strand}.hmm.bed"),
        hmm_annotation = "NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz",
        hmm_igv = "NonCodingRNA/tables/{chrom}.{strand}.hmm.igv.bed.gz",
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
        intersect_hash = ','.join(['$' + str(x) for x in range(7,93)]),
        strand = GetStrandParam
    log:
        "logs/NonCodingRNA/GetTranscriptHMM.{chrom}.{strand}.log"
    shell:
        """
        (zcat {input.counts} | tail -n+2 | awk -F'_' '{{print $1"\\t"$2"\\t"$3"\\t.\\t0\\t{params.strand}\\t"$4}}' - | awk -F'\\t{wildcards.strand}\\t' '{{print $1"\\t"$2}}' - | gzip - > {output.sample_counts}) &> {log};
        echo 'step 1' >> {log};
        (bedtools intersect -s -a {output.sample_counts} -b {input.expression} -wb | awk '{{print $1, $2, $3, $96, $97, $6, {params.intersect_hash} }}' OFS='\\t' - | gzip - > {output.expression_all}) &>> {log};
        echo 'step 2' >> {log};
        (bedtools intersect -s -a {output.sample_counts} -b {input.low_expression} -wb | awk '{{print $1, $2, $3, $96, $97, $6, {params.intersect_hash} }}' OFS='\\t' - | gzip - > {output.low_expression_all}) &>> {log};
        echo 'step 3' >> {log};
        #(bedtools intersect -s -a {output.sample_counts} -b {input.expression} -u | gzip - > {output.expression_all}) &>> {log};
        echo 'step 4' >> {log};
        #(bedtools intersect -s -a {output.sample_counts} -b {input.low_expression} -u | gzip - > {output.low_expression_all}) &>> {log};
        echo 'step 5' >> {log};
        (zcat {input.genes} | grep protein_coding | awk '$1=="{wildcards.chrom}" && $6=="{params.strand}" {{print $0, {params.init_counts} }}' OFS='\\t' -  > {output.pc1}) &>> {log};
        echo 'step 6' >> {log};
        (zcat {input.genes} | grep "snoRNA\\|snRNA" - | grep {wildcards.chrom} - | awk '$1=="{wildcards.chrom}" && $6=="{params.strand}" {{print $0, {params.init_counts} }}' OFS='\\t' - >> {output.pc1}) &>> {log};
        echo 'step 7' >> {log};
        ########
        (bedtools slop -i {output.pc1} -b 2000 -g {input.chrom_sizes} | bedtools intersect -s -a {output.sample_counts} -b - -u | bedtools merge -i - -s -c {params.merge_c} -o first,first,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.overlap_counts}) &>> {log};
        echo 'step 8' >> {log};
        (zcat {output.overlap_counts} >> {output.pc1}) &>> {log};
        echo 'step 9' >> {log};
        (bedtools sort -i {output.pc1} | bedtools merge -s -i - -c {params.merge_c} -o first,first,first,{params.merge_o} -d 1000 | awk '{{print sep='\\t' $1, $2, $3, "pc_" NR, {params.merge_hash} }}' FS='\\t' OFS='\\t' | gzip - > {output.pc}) &>> {log};
        echo 'step 10' >> {log};
        #########
        (bedtools subtract -s -a {output.expression_all} -b {output.pc} | bedtools merge -i - -s -c {params.merge_c} -o first,collapse,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.expression_filtered} )&>> {log};
        echo 'step 11' >> {log};
        (bedtools subtract -s -a {output.low_expression_all} -b {output.pc} | bedtools merge -i - -s -c {params.merge_c} -o first,collapse,first,{params.merge_o} | bedtools sort -i - | gzip - > {output.low_expression_filtered} )&>> {log};
        echo 'step 12' >> {log};
        ##########
        (zcat {output.low_expression_filtered} > {output.filtered_all}) &>> {log};
        echo 'step 13' >> {log};
        (zcat {output.expression_filtered} >> {output.filtered_all}) &>> {log};
        echo 'step 14' >> {log};
        (bedtools sort -i {output.filtered_all} | bedtools merge -s -i - -c {params.merge_c} -o first,collapse,first,{params.merge_o} | bedtools sort -i - | awk '{{print sep='\\t' $1, $2, $3, "ncRNA_" NR, {params.merge_hash} }}' FS='\\t' OFS='\\t' | gzip - > {output.ncRNA_overlaps}) &>> {log};
        echo 'step 15' >> {log};
        (bedtools intersect -s -a {output.ncRNA_overlaps} -b {output.expression_filtered} -u | gzip - > {output.ncRNA_candidates}) &>> {log};
        echo 'step 16' >> {log};
        (zcat {output.ncRNA_candidates} > {output.hmm_unsorted}) &>> {log};
        echo 'step 17' >> {log};
        (zcat {output.pc} >> {output.hmm_unsorted}) &>> {log};
        echo 'step 18' >> {log};
        (bedtools sort -i {output.hmm_unsorted} | gzip - > {output.hmm_annotation}) &>> {log};
        echo 'step 19' >> {log};
        (zcat {output.hmm_annotation} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - > {output.hmm_igv}) &>> {log};
        echo 'finished!' >> {log};
        """
       
rule trim_ncRNAs:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.hmm_trimmed.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 12000
    log:
        "logs/NonCodingRNA/TrimHMM.{chrom}.{strand}.log"
    shell:
        """
        python scripts/NonCodingRNA/trim_ncRNAs.py --chrom {wildcards.chrom} --strand {wildcards.strand} &> {log}
        """
        
        
rule AddCountsToTrimmedNcRNAs:
    input:
        sample_counts = "NonCodingRNA/tables/{chrom}.{strand}.sample.counts.bed.gz",
        hmm = "NonCodingRNA/tables/{chrom}.{strand}.hmm_trimmed.bed.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.hmm_w_counts.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 24000
    log:
        "logs/NonCodingRNA/AddCountsToTrimmedNcRNAs.{chrom}.{strand}.log"
    params:
        intersect_hash = ','.join(['$' + str(x) for x in range(7,93)]),
        merge_c = ','.join([str(x) for x in range(4,93)]),
        merge_o = ','.join(['sum']*86),
    shell:
        """
        (bedtools sort -i {input.hmm} | bedtools merge -s -i - -c 4,5,6 -o first,first,first | bedtools intersect -s -a {input.sample_counts} -b - -wb | awk '{{print $1, $2, $3, $96, $97, $6, {params.intersect_hash} }}' OFS='\\t' - | bedtools merge -s -i - -c {params.merge_c} -o first,first,first,{params.merge_o} | bedtools sort -i - | awk '{{print $1, $2, $3, $4, $5, $6, {params.intersect_hash} }}' OFS='\\t' - | gzip - > {output}) &> {log}
        """
        
def GetReverseCounts(wildcards):
    if wildcards.strand == 'plus':
        rev_strand = 'minus'
    else:
        rev_strand = 'plus'
    rev_counts = "NonCodingRNA/tables/{chrom}.{strand}.sample.counts.bed.gz".format(
    chrom=wildcards.chrom, strand=rev_strand)
    
    return rev_counts

rule GetReverseCoverage:
    input:
        rev_counts = GetReverseCounts,
        hmm = "NonCodingRNA/tables/{chrom}.{strand}.hmm_w_counts.bed.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.reverse_expression.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 24000
    log:
        "logs/NonCodingRNA/GetReverseCoverage.{chrom}.{strand}.log"
    params:
        merge_c = '4,5,6,' + ','.join([str(x) for x in range(13,99)]),
        merge_o = 'first,first,first,' + ','.join(['sum']*86),
    shell:
        """
        zcat {input.hmm} | awk -F'\\t' '{{ print $1, $2, $3, $4, $5, $6}}' OFS='\\t' - | bedtools intersect -S -a - -b {input.rev_counts} -wb | bedtools sort -i - | bedtools merge -s -i - -c {params.merge_c} -o {params.merge_o} | gzip - > {output}
        """

rule FilterAndMergeNcRNAs:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.hmm.sorted.bed.gz",
        "NonCodingRNA/tables/{chrom}.{strand}.hmm_w_counts.bed.gz",
        "NonCodingRNA/tables/{chrom}.{strand}.reverse_expression.bed.gz"
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.hmm.filtered.bed.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
    resources:
        mem_mb = 12000
    conda:
        "../envs/py_tools.yml"
    log:
        "logs/NonCodingRNA/FilterAndMerge.{chrom}.{strand}.log"
    shell:
        """
        python scripts/NonCodingRNA/filter_and_merge_ncRNAs.py --chrom {wildcards.chrom} --strand {wildcards.strand} &> {log}
        """
        
rule MergeNonCodingAnnotation:
    input:
        expand("NonCodingRNA/tables/{chrom}.{strand}.hmm.filtered.bed.gz",
        chrom = chrom_list, strand = ['plus', 'minus'])
    output:
        tmp_bed = temp("NonCodingRNA/tables/NonCodingRNA.merged.bed"),
        bed = "NonCodingRNA/tables/NonCodingRNA.merged.bed.gz"
    log:
        "logs/NonCodingRNA/MergeNonCodingRNA.log"
    params:
        chrom = ' '.join(chrom_list),
        strand = 'minus plus',
    shell:
        """
        echo 'init merge' > {log}
        for CHROM in {params.chrom}
        do
          for STRAND in {params.strand}
          do
            echo $CHROM >> {log}
            echo $STRAND >> {log}
            (zcat NonCodingRNA/tables/${{CHROM}}.${{STRAND}}.hmm.filtered.bed.gz >> {output.tmp_bed}) &>> {log} ;
          done;
        done;
        echo "done merging" >> {log}
        (bedtools sort -i {output.tmp_bed} | gzip - > {output.bed}) &>> {log};
        """
    
        
rule RemovePseudogenesAndLncRNA:
    input:
        genes = "NonCodingRNA/annotation/allGenes.bed.gz",
        bed = "NonCodingRNA/tables/NonCodingRNA.merged.bed.gz",
    output:
        "NonCodingRNA/annotation/ncRNA.bed.gz"
    log:
        "logs/NonCodingRNA/RemovePseudoAndLncRNANonCodingRNA.log"
    shell:
        """
        zcat {input.bed} | awk '{{print sep='\\t' $1, $2, $3, "ncRNA_" NR, $5, $6 }}' FS='\\t' OFS='\\t' - | gzip - > {output};
        #zcat {input.genes} | grep 'pseudogene\\|lncRNA' - | bedtools intersect -v -a {input.bed} -b - -s -f 0.5 | awk '{{print sep='\\t' $1, $2, $3, "ncRNA_" NR, $5, $6 }}' FS='\\t' OFS='\\t' - | gzip - > {output}
        """
        
#########################
# Get other annotations #
#########################

rule GetAllGenesGencode:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    output:
        "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    shell:
        """
        awk -F'\\t' '$3=="gene" {{print $1"\\t"$4"\\t"$5"\\t"$7"\\t"$9}}' {input} | awk -F'"' '{{print $1"\\t"$2"\\t"$4}}' | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4}}' | bedtools sort -i - | gzip - > {output}
        """
        
rule GetTSSAnnotations:
    input:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        #chrom_sizes = "../data/Chrome.sizes",
        #tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        ncRNA_temp = temp("NonCodingRNA/annotation/ncRNA.TSS_temp.bed"),
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
    log:
        "logs/NonCodingRNA/TSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 24000
    shell:
        """
        (awk '{{if ((($3-$2)/10)>1000) a=int(($3-$2)/10); else a=1000}} {{print $1, $2-a, $2+a, $4, $4, $6}}' OFS='\\t' {input.plus_bed} > {output.ncRNA_temp}) &> {log};
        (awk '{{if ((($3-$2)/10)>1000) a=int(($3-$2)/10); else a=1000}} {{print $1, $3-a, $3+a, $4, $4, $6}}' OFS='\\t' {input.minus_bed} >> {output.ncRNA_temp}) &>> {log};
        (bedtools sort -i {output.ncRNA_temp} | gzip - > {output.ncRNA_tss}) &>> {log};
        """

rule GetAllTSSAnnotations:
    input:
        chrom_sizes = "../data/Chrome.sizes",
        tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        all_genes_tss_bp = "NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz",
        all_genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    log:
        "logs/NonCodingRNA/AllTSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 58000
    output:
        all_genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    shell:
        """
        (awk '{{print "chr"$1, $2, $3, $7, $6":_"NR, $4}}' FS='\\t' OFS='\\t' {input.tss_bed} | gzip - > {output.all_genes_tss_bp}) &> {log};
        (bedtools slop -b 999 -i {output.all_genes_tss_bp} -g {input.chrom_sizes} | gzip - > {output.all_genes_tss}) & >> {log};
        """
     
rule Get3PrimeEndAnnotations:
    input:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        utr_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.utr.bed",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
        chrom_sizes = "../data/Chrome.sizes"
    output:
        ncRNA_temp = temp("NonCodingRNA/annotation/ncRNA.3PrimeEnd_temp.bed"),
        genes_temp = temp("NonCodingRNA/annotation/allGenes.3PrimeEnd_temp.bed"),
        genes_last_temp = temp("NonCodingRNA/annotation/allGenes.Last3PrimeEnd_temp.bed"),
        ncRNA_3prime_end = "NonCodingRNA/annotation/ncRNA.3PrimeEnd.bed.gz",
        genes_3prime_end = "NonCodingRNA/annotation/allGenes.3PrimeEnd.bed.gz",
        genes_3prime_last = "NonCodingRNA/annotation/allGenes.Last3PrimeEnd.bed.gz",
    log:
        "logs/NonCodingRNA/EndAnnotations.log",
    params:
        window = 200
    resources:
        mem_mb = 58000
    shell:
        """
        (awk '{{print $1, $3, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' {input.plus_bed} > {output.ncRNA_temp}) &> {log};
        (awk '{{print $1, $2, $2, $4, $5, $6}}' FS='\\t' OFS='\\t' {input.minus_bed} >> {output.ncRNA_temp}) &>> {log};
        (bedtools sort -i {output.ncRNA_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.ncRNA_3prime_end}) &>> {log};
        (awk '$6=="3UTR" && $4=="+" {{print "chr"$1, $3, $3, $5, $7, $4}}' FS='\\t' OFS='\\t' {input.utr_bed} > {output.genes_temp}) &>> {log};
        (awk '$6=="3UTR" && $4=="-" {{print "chr"$1, $2, $2, $5, $7, $4}}' FS='\\t' OFS='\\t' {input.utr_bed} >> {output.genes_temp}) &>> {log};
        (bedtools sort -i {output.genes_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.genes_3prime_end}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="+" {{print $1, $3, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' - > {output.genes_last_temp}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="-" {{print $1, $2, $2, $4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.genes_last_temp}) &>> {log};
        (bedtools sort -i {output.genes_last_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.genes_3prime_last}) &>> {log};
        """

###################################
#        Classify ncRNAs          #
###################################

rule uaRNA:
    input:
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz",
    resources:
        mem_mb = 58000,
    output:
        tmp = temp("NonCodingRNA/annotation/tmp/uaRNA.bed"),
        uaRNA = "NonCodingRNA/annotation/tmp/uaRNA.bed.gz",
    log:
        "logs/NonCodingRNA/ncRNA_annotation/uaRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.genes_tss} -wo > {output.tmp}) &> {log}
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.ncRNA_tss} -wo >> {output.tmp}) &>> {log}
        (bedtools sort -i {output.tmp} | gzip - > {output.uaRNA}) &>> {log}
        """

use rule uaRNA as uaRNA_annotated with:
    input:
        genes_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        ncRNA_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz",
    output:
        tmp = temp("NonCodingRNA/annotation/tmp/uaRNA.annotated.bed"),
        uaRNA = "NonCodingRNA/annotation/tmp/uaRNA.annotated.bed.gz",

rule incRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
    resources:
        mem_mb = 24000,
    output:
        protein_coding = "NonCodingRNA/annotation/tmp/protein_coding.bed.gz",
        incRNA = "NonCodingRNA/annotation/tmp/incRNA.bed.gz",
    log:
        "logs/NonCodingRNA/ncRNA_annotation/incRNAs.log"
    shell:
        """
        (zcat {input.genes_bed} | grep protein_coding - > {output.protein_coding}) &> {log}
        (bedtools intersect -v -a {input.ncRNA_bed} -b {output.protein_coding} -wa | gzip - > {output.incRNA}) &>> {log}
        """
        
rule incRNA_annotated: 
    input:
        lnc = "NonCodingRNA/annotation/tmp/lncRNA.bed",
        protein_coding = "NonCodingRNA/annotation/tmp/protein_coding.bed.gz"
    output:
        incRNA = "NonCodingRNA/annotation/tmp/incRNA.annotated.bed.gz",
    resources:
        mem_mb = 24000,
    log:
        "logs/NonCodingRNA/ncRNA_annotation/incRNAs.annotated.log"
    shell:
        """
        (bedtools intersect -v -a {input.lnc} -b {input.protein_coding} -wa | gzip - > {output.incRNA}) &> {log}
        """

rule GetLncRNAs:
    input:
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
        ncRNA = "NonCodingRNA/annotation/ncRNA.bed.gz"
    output:
        lnc = "NonCodingRNA/annotation/tmp/lncRNA.bed",
        lnc_bed = "NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    log:
        "logs/NonCodingRNA/ncRNA_annotation/getlncRNA.log"
    shell:
        """
        #(zcat {input.genes_bed} | grep 'pseudogene\\|lncRNA\\|snoRNA\\|snRNA' - | gzip - > {output.lnc}) &> {log};
        (zcat {input.genes_bed} | grep lncRNA - | awk '{{print $1, $2, $3, lncRNA_$4, $5, $6}}' FS='\\t' OFS='\\t' - > {output.lnc}) &> {log};
        (zcat {input.genes_bed} | grep pseudogene - | awk '{{print $1, $2, $3, pseudogene_$4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.lnc}) &> {log};
        (zcat {input.genes_bed} | grep snoRNA - | awk '{{print $1, $2, $3, snoRNA_$4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.lnc}) &> {log};
        (zcat {input.genes_bed} | grep snRNA - | awk '{{print $1, $2, $3, snRNA_$4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.lnc}) &> {log};
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -a - -b {output.lnc} -F 0.9 -wo | gzip - > {output.lnc_bed}) &>> {log}
        """
        
rule srtRNA:
    input:
        ncRNA = "NonCodingRNA/annotation/ncRNA.bed.gz",
        lncRNA_bed = "NonCodingRNA/annotation/tmp/lncRNA.bed",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 24000,
    output:
        tmp = temp("NonCodingRNA/annotation/tmp/srtRNA.bed"),
        srtRNA = "NonCodingRNA/annotation/tmp/srtRNA.bed.gz",
        tmp_ = temp("NonCodingRNA/annotation/tmp/rtRNA.bed"),
        rtRNA = "NonCodingRNA/annotation/tmp/rtRNA.bed.gz",
        tmp_a = temp("NonCodingRNA/annotation/tmp/srtRNA.annotated.bed"),
        srtRNA_a = "NonCodingRNA/annotation/tmp/srtRNA.annotated.bed.gz",
        tmp_a_ = temp("NonCodingRNA/annotation/tmp/rtRNA.annotated.bed"),
        rtRNA_a = "NonCodingRNA/annotation/tmp/rtRNA.annotated.bed.gz",
    log:
        "logs/NonCodingRNA/ncRNA_annotation/rtRNAs.log"
    shell:
        """
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a - -b {input.genes_bed} -wo -f 1  > {output.tmp}) &> {log}
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a - -b - -wo -f 1 >> {output.tmp}) &>> {log}
        (bedtools sort -i {output.tmp} | gzip - > {output.srtRNA}) &>> {log};
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a - -b {input.genes_bed} -wo > {output.tmp_}) &>> {log}
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a - -b - -wo >> {output.tmp_}) &>> {log}
        (bedtools sort -i {output.tmp_} | gzip - > {output.rtRNA}) &>> {log}
        (bedtools intersect -S -a {input.lncRNA_bed} -b {input.genes_bed} -wo -f 1  > {output.tmp_a}) &>> {log}
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a {input.lncRNA_bed} -b - -wo -f 1 >> {output.tmp_a}) &>> {log}
        (bedtools sort -i {output.tmp_a} | gzip - > {output.srtRNA_a}) &>> {log};
        (bedtools intersect -S -a {input.lncRNA_bed} -b {input.genes_bed} -wo > {output.tmp_a_}) &>> {log}
        (zcat {input.ncRNA} | awk '{{print $1, $2, $3, $4, $4, $6}}' OFS='\\t' - | bedtools intersect -S -a {input.lncRNA_bed} -b - -wo >> {output.tmp_a_}) &>> {log}
        (bedtools sort -i {output.tmp_a_} | gzip - > {output.rtRNA_a}) &>> {log}
        """




########       ########
########       ########
#######################
#######################
########       ########
########       ########

def GetTSSFile(wildcards):
    if wildcards.Annotation == 'ncRNA':
        return "NonCodingRNA/annotation/ncRNA.TSS.bed.gz"
    elif wildcards.Annotation == 'allGenes':
        return "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
        
def GetAnnotationFile(wildcards):
    if wildcards.Annotation == 'ncRNA':
        return "NonCodingRNA/annotation/ncRNA.bed.gz"
    elif wildcards.Annotation == 'allGenes':
        return "NonCodingRNA/annotation/allGenes.bed.gz"

def GetHistonePeaks(wildcards):
    if wildcards.Histone == 'H3K4ME1':
        return 'PeakCalling/H3K4ME1_peaks.narrowPeak'
    elif wildcards.Histone == 'H3K4ME3':
        return 'PeakCalling/H3K4ME3_peaks.broadPeak'
    elif wildcards.Histone == 'H3K27AC':
        return 'PeakCalling/H3K27AC_peaks.narrowPeak'
    

rule OvelapHistonePeaksWithTSS:
    input:
        TSS_file = GetTSSFile,
        histone_peaks = GetHistonePeaks
    output:
        'NonCodingRNA/annotation/histone_marks/{Annotation}.{Histone}.TSS_overlaps.bed.gz'
    log:
        'logs/NonCodingRNA/{Annotation}.{Histone}.TSS_overlaps.log'
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools intersect -a {input.TSS_file} -b {input.histone_peaks} -wa | bedtools sort -i - | gzip - > {output}) &> {log}
        """

use rule OvelapHistonePeaksWithTSS as OvelapHistonePeaksWithGene with:
    input:
        TSS_file = GetAnnotationFile,
        histone_peaks = GetHistonePeaks
    output:
        'NonCodingRNA/annotation/histone_marks/{Annotation}.{Histone}.overlaps.bed.gz'
    log:
        'logs/NonCodingRNA/{Annotation}.{Histone}.overlaps.log'


rule ClassifyNcRNAs:
    input:
        "NonCodingRNA/annotation/ncRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
        "NonCodingRNA/annotation/tmp/uaRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/rtRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/srtRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/incRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz",
        'NonCodingRNA/annotation/allGenes.TSS.bed.gz',
        'NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz',
        expand("NonCodingRNA/annotation/histone_marks/{Annotation}.{Histone}.{overlap_type}overlaps.bed.gz",
            Annotation = ["ncRNA", "allGenes"], Histone = ["H3K4ME1", "H3K4ME3", "H3K27AC"], 
            overlap_type=["", "TSS_"]
        )
    output:
        'NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz',
        'NonCodingRNA/annotation/allGenes.annotation.tab.gz',
        'NonCodingRNA/annotation/NonCodingRNA.bed.gz',
        'NonCodingRNA/annotation/tmp/tss.saf'
    log:
        "logs/NonCodingRNA/classify_ncRNA.log"
    resources:
        mem_mb = 58000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/NonCodingRNA/make_ncRNA_annotation.py &> {log}
        """



rule MakeSAFAnnotationForNonCodingRNAPhenotype:
    input:
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
    output:
        "NonCodingRNA/annotation/NonCodingRNA.saf",
    log:
        "logs/NonCodingRNA/getSAFforNonCodingRNA.log"
    shell:
        """
        echo -e 'GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output};
        zcat {input} | awk '{{print sep='\\t' $4, $1, $2, $3, $6 }}' FS='\\t' OFS='\\t' >> {output}
        """

rule featureCountsNonCodingRNA:
    input:
        bam = GetBamForPhenotype,
        saf = "NonCodingRNA/annotation/NonCodingRNA.saf",
    output:
        "featureCounts/{Phenotype}_ncRNA/Counts.txt",
    params:
        extraParams = PairedEndParams,
        strandParams = FeatureCountsNonCodingStrandParams
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/NonCodingRNA_annotation/featureCounts/{Phenotype}_ncRNA.log"
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o featureCounts/{wildcards.Phenotype}_ncRNA/Counts.txt {input.bam} &> {log};
        """


rule Prepare_polyA_ExpressionPhenotypes:
    input:
        "featureCounts/polyA.Expression/Counts.txt",
        "featureCounts/polyA.Expression_ncRNA/Counts.txt",
        "featureCounts/polyA.Expression_annotated_ncRNA/Counts.txt",
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
        "NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.CPM.bed.gz",
        "RPKM_tables/polyA.RPKM.bed.gz"
    log:
        "logs/NonCodingRNA/Prepare_polyA_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem = 48000,
    shell:
        """
        mkdir -p RPKM_tables/;
        Rscript scripts/NonCodingRNA/prepare_polyA_ncRNA_Phenotypes.R &> {log}
        """
        
rule Prepare_polyA_Subset_YRI_ExpressionPhenotypes:
    input:
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.CPM.bed.gz",
        "QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz",
    output:
        "QTLs/QTLTools/polyA.Expression_ncRNA.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_polyA_Expression_Phenotypes.Subset_YRI.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem = 48000,
    shell:
        """
        mkdir -p RPKM_tables/;
        Rscript scripts/NonCodingRNA/prepare_polyA_SubsetYRI_ncRNA_Phenotypes.R &> {log}
        """
        
rule Prepare_ml30_ExpressionPhenotypes:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "featureCounts/MetabolicLabelled.30min/Counts.txt",
        "featureCounts/MetabolicLabelled.30min_ncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.30min_annotated_ncRNA/Counts.txt",
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
        "NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/MetabolicLabelled.30min_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "RPKM_tables/MetabolicLabelled.30min.RPKM.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_ml30_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem = 48000,
    shell:
        """
        mkdir -p RPKM_tables/;
        Rscript scripts/NonCodingRNA/prepare_ml30_ncRNA_Phenotypes.R &> {log}
        """
        
        
rule Prepare_ml60_ExpressionPhenotypes:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "featureCounts/MetabolicLabelled.60min/Counts.txt",
        "featureCounts/MetabolicLabelled.60min_ncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.60min_annotated_ncRNA/Counts.txt",
        "NonCodingRNA/annotation/NonCodingRNA.bed.gz",
        "NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/MetabolicLabelled.60min_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "RPKM_tables/MetabolicLabelled.60min.RPKM.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_ml60_Expression_Phenotypes.log"
    resources:
        mem = 48000,
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        mkdir -p RPKM_tables/;
        Rscript scripts/NonCodingRNA/prepare_ml60_ncRNA_Phenotypes.R &> {log}
        """
        
#####################################################
# Get ProCap counts at TSS for validation of diQTLs #
# ###################################################

rule MakeUaRNAForAnnotatedLncRNAs:
    input:
        'NonCodingRNA/annotation/tmp/uaRNA.annotated.bed.gz',
        'RPKM_tables/chRNA.RPKM.bed.gz',
        'NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz',
        'QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz',
        'NonCodingRNA/annotation/tmp/tss.saf',
    output:
        'NonCodingRNA/annotation/Gencode.uaRNA.annotation.tab.gz',
        'NonCodingRNA/annotation/tmp/tss_all.saf',
    log:
        'logs/NonCodingRNA/get_Gencode_uaRNA_annotation.log'
    resources:
        mem_mb = 12000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/NonCodingRNA/make_annotated_uaRNA_TSS.py &> {log}
        """
        


rule MakeSAFAnnotationOfTSSForProCap:
    input:
        SAF='NonCodingRNA/annotation/tmp/tss_all.saf',
    log:
         'logs/NonCodingRNA/get_tss_saf.log'
    output:
        TSS = 'NonCodingRNA/annotation/uaRNA.TSS.saf'
    shell:
        """
        #(cat {input.SAF} > {output.TSS}) &>> {log};
        (head -n 1 {input.SAF} > {output.TSS}) &> {log};
        (tail -n+2 {input.SAF} | awk '{{print $1, $2, $3-800, $4+800, $5}}' OFS='\\t' FS='\\t' - >> {output.TSS}) &>> {log};
        (tail -n+2 {input.SAF} | awk '{{if ($5=="+") a="-"; else a="+"}} {{print $1"_reverse", $2, $3-800, $4+800, a}}' OFS='\\t' FS='\\t' - >> {output.TSS}) &>> {log};
        """


rule MakeProCapTSSQQNorm:
    input:
       'featureCounts/uaRNA_TSS/ProCap/Counts.txt' 
    output:    
        "QTLs/QTLTools/ProCap_uaRNA/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/NonCodingRNA/procap_tss_qqnorm.log"
    resources:
        mem_mb = 12000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/NonCodingRNA/PreparePhenotypeTable_ProCap_uaRNA.R {input} {output} &> {log}
        """











        
