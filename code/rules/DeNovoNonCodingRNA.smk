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
    
rule GetAllGenesGencode:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
        #"ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed"
    output:
        "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    shell:
        """
        awk -F'\\t' '$3=="gene" {{print $1"\\t"$4"\\t"$5"\\t"$7"\\t"$9}}' {input} | awk -F'"' '{{print $1"\\t"$2"\\t"$4}}' | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4}}' | bedtools sort -i - | gzip - > {output}
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
        tmp_ncRNA_bed = "NonCodingRNA/annotation/tmp/ncRNA.2states.bed.gz",
        #tmp_ncRNA_merged = "NonCodingRNA/annotation/tmp/ncRNA_merged.2states.bed.gz",
        tmp_filtered_bed = "NonCodingRNA/annotation/tmp/ncRNA_filtered.2states.bed.gz",
        #unfiltered_ncRNA_bed = "NonCodingRNA/annotation/ncRNA_unfiltered.2states.sorted.bed.gz",
        filtered_bed_name = "NonCodingRNA/annotation/tmp/ncRNA_filtered.2states.sorted_renamed.bed.gz",
        filtered_ncRNA_bed = "NonCodingRNA/annotation/ncRNA.2states.sorted.bed.gz",
        annotated_temp = "NonCodingRNA/annotation/tmp/AnnotatedTemp.bed",
        annotated_HMM = "NonCodingRNA/annotation/AnnotatedHMM.bed.gz",
        temp_merged = "NonCodingRNA/annotation/tmp/allHMM.merged_temp.bed",
        merged_allHMM = "NonCodingRNA/annotation/allHMM.merged.bed.gz",
        protein_coding_HMM = "NonCodingRNA/annotation/protein_coding_HMM.bed.gz",
    params:
        #length = 1000,
        length = 1,
        max_overlap = 0.15
    #wildcard_constraints:
    #    nstates = '2|3'
    log:
        "logs/NonCodingRNA/ncRNA_annotation.2states.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (zcat {input.genes_bed} | grep protein_coding - | bedtools intersect -f 0.1 -s -a {input.allHMM_bed} -b - -wa | bedtools sort -i - | bedtools merge -s -i - -d 1000 -c 4,5,6 -o first | gzip - > {output.protein_coding_HMM}) &> {log};
        (zcat {output.protein_coding_HMM} > {output.annotated_temp}) &>> {log};
        (zcat {input.genes_bed} >> {output.annotated_temp}) &>> {log};
        (bedtools sort -i {output.annotated_temp} | bedtools merge -s -i - -c 4,5,6 -o first | gzip - > {output.annotated_HMM}) &>> {log}
        (bedtools subtract -A -s -a {input.allHMM_bed} -b {output.protein_coding_HMM} -f {params.max_overlap} | bedtools sort -i - | gzip - > {output.tmp_ncRNA_bed}) &>> {log};
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
    resources:
        mem_mb = 32000
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
    resources:
        mem_mb = 32000
    params:
        min_RPKM = -1,
        distance = 10000,
        #min_correlation = 0.5, 
        #RPKM_ratio = 2
        min_correlation = 0.4, 
        RPKM_ratio = 5
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
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
        chrom_sizes = "../data/Chrome.sizes",
        tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        ncRNA_temp = temp("NonCodingRNA/annotation/ncRNA.TSS_temp.bed"),
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        genes_temp = temp("NonCodingRNA/annotation/allGenes.FirstTSS_temp.bed"),
        genes_tss = "NonCodingRNA/annotation/allGenes.FirstTSS.bed.gz",
        #all_genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    log:
        "logs/TSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 58000
    shell:
        """
        (awk -F'\\t' '{{print $1"\\t"$2"\\t"$2"\\t"$4"\\t"$5"\\t"$6}}' {input.plus_bed} > {output.ncRNA_temp}) &> {log};
        (awk -F'\\t' '{{print $1"\\t"$3"\\t"$3"\\t"$4"\\t"$5"\\t"$6}}' {input.minus_bed} >> {output.ncRNA_temp}) &>> {log};
        (bedtools sort -i {output.ncRNA_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.ncRNA_tss}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="+" {{print $1, $2, $2, $4, $5, $6}}' FS='\\t' OFS='\\t' - > {output.genes_temp}) &>> {log};
        (zcat {input.genes_bed} | awk '$6=="-" {{print $1, $3, $3, $4, $5, $6}}' FS='\\t' OFS='\\t' - >> {output.genes_temp}) &>> {log};
        (bedtools sort -i {output.genes_temp} | bedtools slop -b {params.window} -g {input.chrom_sizes} -i - | gzip - > {output.genes_tss}) &>> {log};
        """
#(awk '{{print "chr"$1, $2, $3, $6, $7, $4}}' FS='\\t' OFS='\\t' {input.tss_bed} | bedtools slop -b 999 -i - -g {input.chrom_sizes} | gzip - #> {output.all_genes_tss}) & >> {log};
 
rule GetAllTSSAnnotations:
    input:
        chrom_sizes = "../data/Chrome.sizes",
        tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        all_genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    log:
        "logs/AllTSSAnnotations.log",
    params:
        window = 1000
    resources:
        mem_mb = 58000
    output:
        all_genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz"
    shell:
        """
        (awk '{{print "chr"$1, $2, $3, $7, $6, $4}}' FS='\\t' OFS='\\t' {input.tss_bed} | bedtools slop -b 999 -i - -g {input.chrom_sizes} | gzip - > {output.all_genes_tss}) & >> {log};
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
        "logs/EndAnnotations.log",
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
        mem_mb = 12000,
        cpus_per_node = 9,
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o NonCodingRNA/Expression_HMM/{wildcards.Phenotype}_featureCounts/Counts.txt {input.bam} &> {log};
        """
        
#################################################

rule uaRNA:
    input:
        ncRNA_tss = "NonCodingRNA/annotation/ncRNA.TSS.bed.gz",
        genes_tss = "NonCodingRNA/annotation/allGenes.TSS.bed.gz",
    resources:
        mem_mb = 58000,
    output:
        uaRNA = "NonCodingRNA/annotation/tmp/uaRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/uaRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.genes_tss} -wo > NonCodingRNA/annotation/tmp/uaRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.ncRNA_tss} -wo >> NonCodingRNA/annotation/tmp/uaRNA.bed) &>> {log}
        (gzip NonCodingRNA/annotation/tmp/uaRNA.bed) &>> {log}
        """

rule rtRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        rtRNA = "NonCodingRNA/annotation/tmp/rtRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/rtRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_bed} -wo > NonCodingRNA/annotation/tmp/rtRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_bed} -wo >> NonCodingRNA/annotation/tmp/rtRNA.bed) &>> {log}
        (gzip NonCodingRNA/annotation/tmp/rtRNA.bed) &>> {log}
        """

rule srtRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        rtRNA = "NonCodingRNA/annotation/tmp/srtRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/srtRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_bed} -wo -f 1  > NonCodingRNA/annotation/tmp/srtRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_bed} -wo -f 1 >> NonCodingRNA/annotation/tmp/srtRNA.bed) &>> {log}
        (gzip NonCodingRNA/annotation/tmp/srtRNA.bed) &>> {log}
        """

rule coRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        ncRNA_3prime_end = "NonCodingRNA/annotation/ncRNA.3PrimeEnd.bed.gz",
        genes_3prime_last = "NonCodingRNA/annotation/allGenes.Last3PrimeEnd.bed.gz",
    resources:
        mem_mb = 58000,
    output:
        coRNA = "NonCodingRNA/annotation/tmp/coRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/coRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_3prime_last} -wo > NonCodingRNA/annotation/tmp/coRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_3prime_end} -wo >> NonCodingRNA/annotation/tmp/coRNA.bed) &>> {log}
        (gzip NonCodingRNA/annotation/tmp/coRNA.bed) &>> {log}
        """

rule incRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        incRNA = "NonCodingRNA/annotation/tmp/incRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/incRNAs.log"
    shell:
        """
        (bedtools intersect -v -a {input.ncRNA_bed} -b {input.genes_bed} -wa | gzip - > {output.incRNA}) &> {log}
        """
        
rule ctRNA:
    input:
        ncRNA_bed = "NonCodingRNA/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        ctRNA = "NonCodingRNA/annotation/tmp/ctRNA.bed.gz",
    log:
        "logs/ncRNA_annotation/ctRNAs.log"
    shell:
        """
        (bedtools intersect -s -a {input.ncRNA_bed} -b {input.genes_bed} -wo | gzip - > {output.ctRNA}) &> {log}
        """
        
rule GetLncRNAs:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA/annotation/tmp/lncRNA.bed.gz",
        out = "NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    log:
        "logs/ncRNA_annotation/getlncRNA.log"
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' | awk '$7=="lncRNA" {{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
         (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule GetSnoRNAs:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA/annotation/tmp/snoRNA.bed.gz",
        out = "NonCodingRNA/annotation/tmp/snoRNA.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    log:
        "logs/ncRNA_annotation/getsnoRNA.log"
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' | awk '$7=="snoRNA" {{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
         (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule GetPsuedogenes:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA/annotation/tmp/pseudogenes.bed.gz",
        out = "NonCodingRNA/annotation/tmp/pseudogenes.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    log:
        "logs/ncRNA_annotation/getPseudogenes.log"
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' - | grep pseudogene - | awk '{{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
        (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule ClassifyNcRNAs:
    input:
        "NonCodingRNA/annotation/ncRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz",
        "NonCodingRNA/annotation/tmp/uaRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/rtRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/srtRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/snoRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/coRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/incRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/ctRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/lncRNA.ncRNA.bed.gz",
        "NonCodingRNA/annotation/tmp/pseudogenes.ncRNA.bed.gz",
    output:
        "NonCodingRNA/annotation/ncRNA.annotation.tab.gz"
    log:
        "logs/classify_ncRNA.log"
    resources:
        mem_mb = 58000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/classify_ncRNA.py &> {log}
        """
        
        
# coRNA: colliding RNA. ncRNA overlaps 3' end of the gene on the reverse strand, but not the entire gene.
# ctRNA: cotranscribed RNA. Less than 15% of a ncRNA overlaps an annotated gene in the sense strand.
# lncRNA: long non-coding RNA. > 90% of a Gencode annotated lncRNA is contained within the ncRNA, but this
#         corresponds to < 15% of the ncRNA. More than one lncRNA can be associated to a ncRNA.
# pseudogene: same as lncRNA, but with annotated pseudogenes.
# rtRNA: reverse transcribed RNA. A ncRNA is mostly contained within the reverse strand of a gene, but it's not
#        already classified as another class of ncRNA. rna_type is srtRNA (strict rtRNA) if the entire ncRNA
#        is contained within the gene.
# uaRNA: upstream antisense RNA. ncRNA TSS is within 1000 bp of another gene or ncRNA, on the reverse strand.
#        A ncRNA can be close to more than one annotated gene. If an annotated gene has two or more ncRNAs,
#        only the closest ncRNA is annotated as uaRNA for that gene.

#rule AnnotateNonCodingRNA:









