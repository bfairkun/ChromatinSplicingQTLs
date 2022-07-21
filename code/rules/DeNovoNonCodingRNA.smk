rule MakeWindows:
    input:
        "../data/Chrome.sizes"
    output:
        plus = "NonCodingRNA/bed/{chrom}.windows.plus.bed.gz",
        minus = "NonCodingRNA/bed/{chrom}.windows.minus.bed.gz",
    params:
        win_size = 200,#50,
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
        max_junc_length = 20000
    shell:
        """
        (samtools view -h -F256 {input} | perl -lane 'print if (( $F[0] =~ /^@/ ) || ( abs($F[8]) <= {params.max_junc_length}  ))' | samtools view -bh > {output}) &> {log}
        """
        
#python scripts/filter_long_junc.py --input {input} --output {output} --max_junc_length {params.max_junc_length} &> {log}
        
        
rule IndexBam:
    input:
        "NonCodingRNA/bam/{IndID}.filtered.bam"
    output:
        temp("NonCodingRNA/bam/{IndID}.filtered.bam.bai")
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


rule MakeInputForHMM:
    input:
        expand("NonCodingRNA/counts/chRNA/{{chrom}}.{IndID}.{{strand}}.bed.gz",
               IndID = [x for x in chRNASeqSamples if x != "NA18855"]),
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
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
        #"ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
    output:
        "NonCodingRNA/annotation/allGenes.bed.gz"
    shell:
        """
        awk -F'\\t' '$3=="gene" {{print $1"\\t"$4"\\t"$5"\\t"$7"\\t"$9}}' {input} | awk -F'"' '{{print $1"\\t"$6"\\t"$4}}' | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4}}' | bedtools sort -i - | gzip - > {output}
        """
        #awk -F'\\t' '{{print "chr"$1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$4}}' {input} | bedtools sort -i - | gzip -  > {output}
    
def GetMergedDistance(wildcards):
    if wildcards.nstates == '3':
        return 1000
    elif wildcards.nstates == '2':
        return 2000

rule ProcessAnnotationsHMM:
    input:
        pred = "NonCodingRNA/tables/{chrom}.{strand}.predicted_{nstates}states.tab.gz",
        #genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz"
    output:
        TU_bed = "NonCodingRNA/tables/{chrom}.{strand}.TU_{nstates}states.tab.gz",
        #merged_bed = "NonCodingRNA/tables/{chrom}.{strand}.merged_{nstates}states.tab.gz",
        #ncRNA_bed = "NonCodingRNA/tables/{chrom}.{strand}.{nstates}states.ncRNA.bed.gz",
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        nstates = '2|3'
    params:
        strand = getStrandString,
        #merge_distance = 1000,#GetMergedDistance,
        #max_overlap = 0.1
    log:
        "logs/NonCodingRNA/process_hmm.{chrom}.{strand}.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/processHMM.py --input {input.pred} --output {output.TU_bed} --nstates {wildcards.nstates} &> {log};
        #(bedtools merge -i {output.TU_bed} -d {params.merge_distance} | awk -F'\\t' '{{print $0"\\t"$1"_"$2"_"$3"_{params.strand}\\t"$1"_"$2"_"$3"_{params.strand}\\t{params.strand}"}}' - | sort -u | gzip - > {output.merged_bed}) &>> {log};
        #(bedtools subtract -A -s -a {output.merged_bed} -b {input.genes_bed} -f {params.max_overlap} | sort -u | gzip - > {output.ncRNA_bed}) &>> {log}
        """
   
rule MergeNonCodingRNA:
    input:
        "NonCodingRNA/tables/{chrom}.{strand}.TU_{nstates}states.tab.gz",
    output:
        "NonCodingRNA/tables/{chrom}.{strand}.merged_{nstates}states.bed.gz",
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        nstates = '2|3'
    params:
        strand = getStrandString,
        merge_distance = 0,
        #merge_distance = 1000,
    log:
        "logs/NonCodingRNA/merge_hmm.{chrom}.{strand}.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools merge -i {input} -d {params.merge_distance} | awk -F'\\t' '{{print $0"\\t"$1"_"$2"_"$3"_{params.strand}\\t"$1"_"$2"_"$3"_{params.strand}\\t{params.strand}"}}' - | sort -u | gzip - > {output}) &>> {log};
        """
        
        
rule MergeChromosomesNonCodingRNA:
    input:
        expand("NonCodingRNA/tables/{chrom}.{strand}.merged_{{nstates}}states.bed.gz",
                chrom = chrom_list, strand = ['plus', 'minus']),
    output:
        tmp_hmm_bed = temp("NonCodingRNA/annotation/allHMM.{nstates}states.bed.gz"),
        hmm_bed = "NonCodingRNA/annotation/allHMM.{nstates}states.sorted.bed.gz",
    wildcard_constraints:
        nstates = '2|3'
    log:
        "logs/NonCodingRNA/hmm_annotation.{nstates}states.log"
    resources:
        mem_mb = 58000
    shell:
        """
        python scripts/merge_HMM_annotation.py --nstates {wildcards.nstates} &> {log};
        (bedtools sort -i {output.tmp_hmm_bed} | gzip - > {output.hmm_bed}) &>> {log};
        """
        
rule GetNonCodingRNAFromHMM:
    input:
        genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz",
        allHMM_bed = "NonCodingRNA/annotation/allHMM.2states.sorted.bed.gz",
    output:
        tmp_ncRNA_bed = temp("NonCodingRNA/annotation/ncRNA.2states.bed.gz"),
        #tmp_ncRNA_merged = temp("NonCodingRNA/annotation/ncRNA_merged.2states.bed.gz"),
        tmp_filtered_bed = temp("NonCodingRNA/annotation/ncRNA_filtered.2states.bed.gz"),
        #unfiltered_ncRNA_bed = temp("NonCodingRNA/annotation/ncRNA_unfiltered.2states.sorted.bed.gz"),
        filtered_bed_name = temp("NonCodingRNA/annotation/ncRNA_filtered.2states.sorted_renamed.bed.gz"),
        filtered_ncRNA_bed = "NonCodingRNA/annotation/ncRNA.2states.sorted.bed.gz",
        annotated_temp = temp("NonCodingRNA/annotation/AnnotatedTemp.bed"),
        annotated_HMM = "NonCodingRNA/annotation/AnnotatedHMM.bed.gz",
        temp_merged = temp("NonCodingRNA/annotation/allHMM.merged_temp.bed"),
        merged_allHMM = "NonCodingRNA/annotation/allHMM.merged.bed.gz",
    params:
        #length = 1000,
        length = 1,
        max_overlap = 0.1
    #wildcard_constraints:
    #    nstates = '2|3'
    log:
        "logs/NonCodingRNA/ncRNA_annotation.2states.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (bedtools intersect -f 0.1 -s -a {input.allHMM_bed} -b {input.genes_bed} -wa | bedtools sort -i - | bedtools merge -s -i - -d 1000 -c 4,5,6 -o first  > {output.annotated_temp}) &>> {log};
        (zcat {input.genes_bed} >> {output.annotated_temp}) &>> {log};
        (bedtools sort -i {output.annotated_temp} | bedtools merge -s -i - -c 4,5,6 -o first | gzip - > {output.annotated_HMM}) &>> {log}
        (bedtools subtract -A -s -a {input.allHMM_bed} -b {output.annotated_HMM} -f {params.max_overlap} | bedtools sort -i - | gzip - > {output.tmp_ncRNA_bed}) &>> {log};
        (bedtools merge -s -i {output.tmp_ncRNA_bed} -d 1000 -c 4,5,6 -o first | bedtools subtract -A -s -a - -b {output.annotated_HMM} -f {params.max_overlap} | bedtools sort -i - | gzip - > {output.tmp_filtered_bed}) &>> {log};
        (bedtools sort -i {output.tmp_filtered_bed} | gzip - > {output.filtered_bed_name}) &>> {log};
        (zcat {output.filtered_bed_name} | awk '{{print sep='\\t' $1, $2, $3, "ncRNA_" NR, $5, $6 }}' FS='\\t' OFS='\\t' | gzip - > {output.filtered_ncRNA_bed}) &>> {log};
        (zcat {output.filtered_ncRNA_bed} > {output.temp_merged}) &>> {log};
        (zcat {output.annotated_HMM} >> {output.temp_merged}) &>> {log};
        (bedtools sort -i {output.temp_merged} | gzip - > {output.merged_allHMM}) &>> {log}
        """
        #(bedtools merge -s -i {output.tmp_ncRNA_bed} -d 1000 -c 4,5,6 -o first | sort -u | gzip - > {output.tmp_ncRNA_merged}) &>> {log};
        #(bedtools sort -i {output.tmp_ncRNA_merged} | gzip - > {output.unfiltered_ncRNA_bed}) &>> {log};
        #(zcat {output.unfiltered_ncRNA_bed} | awk '$3-$2>={params.length}' - > {output.tmp_filtered_bed}) &>> {log} 
        
#rule GetCombinedBed:
#    input:
#        pc = "QTLs/QTLTools/chRNA.Expression.Splicing/OnlyFirstReps.RPKM.bed.gz",
#        nc = "QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.RPKM.bed.gz",
#    output:
#        unsorted = temp("NonCodingRNA/annotation/allTranscriptsUnsorted.bed"),
#        allTranscripts = "NonCodingRNA/annotation/allTranscripts.bed.gz"
#    log:
#        "logs/combineAnnotations.log"
#    shell:
#        """
#        (zcat {input.pc} | tail -n +2 - | awk '{{ print $1, $2, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' > {output.unsorted}) &>> {log} ;
#        (zcat {input.nc} | tail -n +2 - | awk '{{ print $1, $2, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' >> {output.unsorted}) &>> {log} ;
#        (bedtools sort -i {output.unsorted} | gzip - > {output.allTranscripts}) &>> {log} 
#        """
        

rule GetRPKMForAllHMM:
    input:
        "NonCodingRNA/Expression_HMM/chRNA.Expression_featureCounts/Counts.txt",
    output:
        "NonCodingRNA/Expression_HMM/OnlyFirstReps.RPKM.bed.gz",
        "NonCodingRNA/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz",
    log:
        "logs/PrepareAllHMMRPKM.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/AllHMM_RPKM.R &> {log}
        """
        
        
rule FilterNonCodingRNA:
    input:
        "NonCodingRNA/annotation/allHMM.merged.bed.gz",
        "NonCodingRNA/Expression_HMM/OnlyFirstReps.RPKM.bed.gz",
        "NonCodingRNA/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz",
    output:
        "NonCodingRNA/annotation/ncRNA.bed.gz",
    log:
        "logs/filter_ncRNAs.log"
    conda:
        "../envs/py_tools.yml"
    params:
        min_RPKM = 0.1,
        distance = 10000,
        min_correlation = 0.5, 
        RPKM_ratio = 2
    shell:
        """
        python scripts/filter_ncRNAs.py --min_RPKM {params.min_RPKM} --distance {params.distance} --min_correlation {params.min_correlation} --RPKM_ratio {params.RPKM_ratio} &> {log}
        """
        
        
#rule QQnormNonCodingRNAs:
#    input:
#        "QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.RPKM.filtered.bed.gz"
#    output:
#        "QTLs/QTLTools/chRNA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz"
#    log:
#        "logs/qqnorm_ncRNAs.log"
#    conda:
#        "../envs/r_essentials.yml"
#    shell:
#        """
#        Rscript scripts/qqnorm_ncRNAs.R &> {log}
#        """
    


#rule ProteinCodingOverlap:
#    input:
#        gene_list = "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
#        hmm_bed = "NonCodingRNA/annotation/allHMM.2states.sorted.bed.gz"
#    output:
#        annotated_hmm = "NonCodingRNA/annotation/AnnotatedHMM.bed.gz",
#        protein_coding = "NonCodingRNA/annotation/AnnotatedProteinCoding_overlap.bed.gz",
#        de_novo = "NonCodingRNA/annotation/DeNovoProteinCoding_overlap.bed.gz"
#    log:
#        "logs/NonCodingRNA/hmm_protein_coding.log"
#    resources:
#        mem_mb = 12000
#    shell:
#        """
#        (bedtools sort -i {input.hmm_bed} | bedtools merge -s -i - -d 1000 -c 4,5,6 -o first | gzip - > {output.annotated_hmm}) &> {log};
#        (bedtools intersect -s -a {input.gene_list} -b {output.annotated_hmm} -f 0.5 -wa | sort -u | gzip - > {output.protein_coding}) &>>{log};
#        (bedtools intersect -s -b {input.gene_list} -a {output.annotated_hmm} -f 0.5 -wa | sort -u | gzip - > {output.de_novo}) &>> {log}
#        """
    
rule SplitAnnotationByStrand:
    input:
        "NonCodingRNA/annotation/ncRNA.bed.gz"
    output:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
    log:
        "logs/NonCodingRNA/split_annotation.log"
    shell:
        """
        (zcat NonCodingRNA/annotation/ncRNA.bed.gz | awk -F'\\t' '$6=="+"' -  > {output.plus_bed}) &> {log};
        (zcat NonCodingRNA/annotation/ncRNA.bed.gz | awk -F'\\t' '$6=="-"' -  > {output.minus_bed}) &>> {log};
        """
        
        
rule GetTSSAnnotations:
    input:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz",
        chrom_sizes = "../data/Chrome.sizes"
    output:
        ncRNA_temp = temp("NonCodingRNA/annotation/ncRNA.TSS_temp.bed"),
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        genes_temp = temp("NonCodingRNA/annotation/allGenes.FirstTSS_temp.bed"),
        genes_tss = "NonCodingRNA/annotation/allGenes.FirstTSS.bed.gz"
    log:
        "logs/TSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 12000
    shell:
        """
        (awk -F'\\t' '{{print $1"\\t"$2"\\t"$2"\\t"$4"\\t"$5"\\t"$6}}' {input.plus_bed} > {output.ncRNA_temp}) &> {log};
        (awk -F'\\t' '{{print $1"\\t"$3"\\t"$3"\\t"$4"\\t"$5"\\t"$6}}' {input.minus_bed} >> {output.ncRNA_temp}) &>> {log};
        (bedtools sort -i {output.ncRNA_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.ncRNA_tss}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="+" {{print $1, $2, $2, $4, $5, $6}}' FS='\\t' OFS='\\t' - > {output.genes_temp}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="-" {{print $1, $3, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.genes_temp}) &>> {log};
        (bedtools sort -i {output.genes_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.genes_tss}) &>> {log};
        """
        
       
rule MakeSAFAnnotationForNonCodingRNA:
    input:
        "NonCodingRNA/annotation/ncRNA.bed.gz",
    output:
        "NonCodingRNA/annotation/ncRNA.saf",
    shell:
        """
        echo -e 'GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output};
        zcat {input} | awk '{{print sep='\\t' $4, $1, $2, $3, $6 }}' FS='\\t' OFS='\\t' >> {output}
        """
        
rule featureCountsNonCodingRNA:
    input:
        bam = GetBamForPhenotype,
        saf = "NonCodingRNA/annotation/ncRNA.saf",
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
        "logs/featureCounts/{Phenotype}_ncRNA.log"
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o featureCounts/{wildcards.Phenotype}_ncRNA/Counts.txt {input.bam} &> {log};
        """


rule Get3PrimeEndAnnotations:
    input:
        plus_bed = "NonCodingRNA/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA/deeptools/bed/ncRNA.minus.bed",
        utr_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.utr.bed",
        genes_bed = "NonCodingRNA/annotation/allGenes.bed.gz",
        chrom_sizes = "../data/Chrome.sizes"
    output:
        ncRNA_temp = temp("NonCodingRNA/annotation/ncRNA.3PrimeEnd_temp.bed"),
        genes_temp = temp("NonCodingRNA/annotation/allGenes.3PrimeEnd_temp.bed"),
        genes_last_temp = temp("NonCodingRNA/annotation/allGenes.Last3PrimeEnd_temp.bed"),
        ncRNA_3prime_end = "NonCodingRNA/annotation/ncRNA.3PrimeEnd.bed.gz",
        genes_3prime_end = "NonCodingRNA/annotation/allGenes.3PrimeEnd.bed.gz",
        genes_3prime_last = "NonCodingRNA/annotation/allGenes.Last3PrimeEnd.bed.gz",
    log:
        "logs/EndAnnotations.log",
    params:
        window = 200
    resources:
        mem_mb = 24000
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

#rule GetUpstreamAntisenseRNAs:
#    input:
#        genes_tss = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed"
#        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.2states.TSS.bed.gz"
#    output:
#    log:
#    params:
#    shell:
#        """
#        awk '{print "chr"$1, $2-999, $3+999, $6, $7, $4}' FS='\t' OFS='\t' #ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed | gzip - > notebooks_scratch/allGenesTSS.bed.gz
#        bedtools intersect -S -a notebooks_scratch/allGenesTSS.bed.gz -b notebooks_scratch/allGenesTSS.bed.gz -wb | uniq | gzip - > notebooks_scratch/allGenes.TSS.Antisense.bed.gz
#        bedtools intersect -S -a NonCodingRNA/annotation/ncRNA.2states.TSS.bed.gz -b notebooks_scratch/allGenesTSS.bed.gz -wb | uniq | gzip - > notebooks_scratch/allGenes.ncRNA.TSS.Antisense.bed.gz
#        bedtools intersect -S -a NonCodingRNA/annotation/ncRNA.2states.TSS.bed.gz -b NonCodingRNA/annotation/ncRNA.2states.TSS.bed.gz -wb | uniq | gzip - > notebooks_scratch/ncRNA.TSS.Antisense.bed.gz
#"""
 
use rule MakeSAFAnnotationForNonCodingRNA as MakeSAFAnnotationForAllHMM with:
    input:
        "NonCodingRNA/annotation/allHMM.merged.bed.gz",
    output:
        "NonCodingRNA/annotation/allHMM.merged.saf",
    log:
        "logs/getSAFforAllHMM.log"


        
rule featureCountsAllHMM:
    input:
        bam = GetBamForPhenotype,
        saf = "NonCodingRNA/annotation/allHMM.merged.saf",
    output:
        "NonCodingRNA/Expression_HMM/{Phenotype}_featureCounts/Counts.txt",
    log:
        "logs/featureCounts_for_allHMM{Phenotype}.log"
    params:
        extraParams = PairedEndParams,
        strandParams = FeatureCountsNonCodingStrandParams
    threads:
        8
    wildcard_constraints:
        Phenotype = "chRNA.Expression",
    resources:
        mem = 12000,
        cpus_per_node = 9,
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o NonCodingRNA/Expression_HMM/{wildcards.Phenotype}_featureCounts/Counts.txt {input.bam} &> {log};
        """
