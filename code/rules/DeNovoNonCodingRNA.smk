#######################################
#    Prepare input for HMM and run    #
#######################################

def GetWinSizePerHMM_pass(wildcards):
    if wildcards.HMM_pass == '':
        return 200
    elif wildcards.HMM_pass == '_merged':
        return 50

rule MakeWindows:
    input:
        "../data/Chrome.sizes"
    output:
        plus = "NonCodingRNA{HMM_pass}/bed/{chrom}.windows.plus.bed.gz",
        minus = "NonCodingRNA{HMM_pass}/bed/{chrom}.windows.minus.bed.gz",
    params:
        win_size = GetWinSizePerHMM_pass,
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        HMM_pass = '_annotation|_merged'
    shell:
        """
        bedtools makewindows -g {input} -w {params} | awk '$1=="{wildcards.chrom}" {{print $1"\\t"$2"\\t"$3"\\t.\\t.\\t+"}}' > {output.plus};
        bedtools makewindows -g {input} -w {params} | awk '$1=="{wildcards.chrom}" {{print $1"\\t"$2"\\t"$3"\\t.\\t.\\t-"}}' > {output.minus}
        """
        
rule FilterLongJunc:
    input:
        "Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam",
    output:
        temp("NonCodingRNA{HMM_pass}/bam/{IndID}.filtered.bam")
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        HMM_pass = '_annotation|_merged'
    resources:
        mem_mb = 12000
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/filter_bam.{IndID}.log"
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
        "NonCodingRNA{HMM_pass}/bam/{IndID}.filtered.bam"
    output:
        temp("NonCodingRNA{HMM_pass}/bam/{IndID}.filtered.bam.bai")
    log:
        "logs/NonCodingRNA{HMM_pass}/index.{IndID}_filtered.bam.log"
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
        bam = "NonCodingRNA{HMM_pass}/bam/{IndID}.filtered.bam",
        bai = "NonCodingRNA{HMM_pass}/bam/{IndID}.filtered.bam.bai"
    output:
        temp("NonCodingRNA{HMM_pass}/bam/{IndID}.{chrom}.bam")
    wildcard_constraints:
        IndID = "|".join(chRNASeqSamples),
        chrom = "|".join(chrom_list),
        HMM_pass = '_annotation|_merged'
    resources:
        mem_mb = 58000
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/split_bam.{IndID}.{chrom}.log"
    shell:
        """
        (samtools view -bh -F 256 -f 64 {input.bam} {wildcards.chrom} > {output}) &> {log}
        """

rule WindowCountsChRNA:
    input:
        faidx = "../data/Chrome.sizes",
        bam = "NonCodingRNA{HMM_pass}/bam/{IndID}.{chrom}.bam",
        chrom_bed = "NonCodingRNA{HMM_pass}/bed/{chrom}.windows.{strand}.bed.gz"
    output:
        "NonCodingRNA{HMM_pass}/counts/chRNA/{chrom}.{IndID}.{strand}.bed.gz"
    wildcard_constraints:
        IndID = "|".join([x for x in chRNASeqSamples if x != "NA18855"]),
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        HMM_pass = '_annotation|_merged'
    resources:
        mem_mb = 58000
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/counts_chRNA.{IndID}.{chrom}.{strand}.log"
    shell:
        """
        (samtools view -bh -F 256 -f 64 {input.bam} | bedtools intersect -sorted -S -g {input.faidx} -a {input.chrom_bed} -b - -c | gzip - > {output}) &> {log}
        """

def GetMergeFlag(wildcards):
    if wildcards.HMM_pass == '_merged':
        return '--merge'
    else:
       return ''

rule MakeInputForHMM:
    input:
        expand("NonCodingRNA{{HMM_pass}}/counts/chRNA/{{chrom}}.{IndID}.{{strand}}.bed.gz",
               IndID = [x for x in chRNASeqSamples if x != "NA18855"]),
    output:
        "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.counts.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        HMM_pass = '_annotation|_merged'
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/make_input_for_hmm.{chrom}.{strand}.log"
    resources:
        mem_mb = 58000
    params:
        extra_params = GetMergeFlag
    shell:
        """
        python scripts/MakeInputForHMM.py --chrom {wildcards.chrom} --counts_dir NonCodingRNA{wildcards.HMM_pass}/counts/ --strand {wildcards.strand} --output {output} {params.extra_params} &> {log}
        """

rule RunHMM:
    input:
        "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.counts.tab.gz"
    output:
        "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.predicted_{states}states.tab.gz"
    wildcard_constraints:
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        states = '2|3',
        HMM_pass = '_annotation|_merged'
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/run_hmm.{chrom}.{strand}.{states}.log"
    resources:
        mem_mb = 58000
    shell:
        """
        Rscript scripts/runHMM.R {wildcards.chrom} {wildcards.strand} {wildcards.states} NonCodingRNA{wildcards.HMM_pass}/tables/ &> {log};
        gzip NonCodingRNA{wildcards.HMM_pass}/tables/{wildcards.chrom}.{wildcards.strand}.predicted_{wildcards.states}states.tab
        """



########################################
#          Process HMM output          #
########################################


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
        "NonCodingRNA{HMM_pass}/annotation/allGenes.bed.gz"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        awk -F'\\t' '$3=="gene" {{print $1"\\t"$4"\\t"$5"\\t"$7"\\t"$9}}' {input} | awk -F'"' '{{print $1"\\t"$6"\\t"$4}}' | awk -F'\\t' '{{print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4}}' | bedtools sort -i - | gzip - > {output}
        """
        #awk -F'\\t' '{{print "chr"$1"\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t"$4}}' {input} | bedtools sort -i - | gzip -  > {output}
    
rule GetAllGenesGencode:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    output:
        "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
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
        pred = "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.predicted_{nstates}states.tab.gz",
    output:
        TU_bed = "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.TU_{nstates}states.tab.gz",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        nstates = '2|3'
    params:
        strand = getStrandString,
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/process_hmm.{chrom}.{strand}.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/processHMM.py --input {input.pred} --output {output.TU_bed} --nstates {wildcards.nstates} &> {log};
        """
   
rule MergeNonCodingRNA:
    input:
        "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.TU_{nstates}states.tab.gz",
    output:
        "NonCodingRNA{HMM_pass}/tables/{chrom}.{strand}.merged_{nstates}states.bed.gz",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
        chrom = "|".join(chrom_list),
        strand = 'minus|plus',
        nstates = '2|3'
    params:
        strand = getStrandString,
        merge_distance = 0,
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/merge_hmm.{chrom}.{strand}.{nstates}states.log"
    resources:
        mem_mb = 12000
    shell:
        """
        (bedtools merge -i {input} -d {params.merge_distance} | awk -F'\\t' '{{print $0"\\t"$1"_"$2"_"$3"_{params.strand}\\t"$1"_"$2"_"$3"_{params.strand}\\t{params.strand}"}}' - | sort -u | gzip - > {output}) &>> {log};
        """
        
        
rule MergeChromosomesNonCodingRNA:
    input:
        expand("NonCodingRNA{{HMM_pass}}/tables/{chrom}.{strand}.merged_{{nstates}}states.bed.gz",
                chrom = chrom_list, strand = ['plus', 'minus']),
    output:
        tmp_hmm_bed = temp("NonCodingRNA{HMM_pass}/annotation/allHMM.{nstates}states.bed.gz"),
        hmm_bed = "NonCodingRNA{HMM_pass}/annotation/allHMM.{nstates}states.sorted.bed.gz",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
        nstates = '2|3'
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/hmm_annotation.{nstates}states.log"
    resources:
        mem_mb = 58000
    params:
        GetMergeFlag
    shell:
        """
        python scripts/merge_HMM_annotation.py --nstates {wildcards.nstates} {params} &> {log};
        (bedtools sort -i {output.tmp_hmm_bed} | gzip - > {output.hmm_bed}) &>> {log};
        """
        
rule GetNonCodingRNAFromHMM:
    input:
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/allGenes.bed.gz",
        allHMM_bed = "NonCodingRNA{HMM_pass}/annotation/allHMM.2states.sorted.bed.gz",
    output:
        tmp_ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/ncRNA.2states.bed.gz",
        #tmp_ncRNA_merged = "NonCodingRNA{HMM_pass}/annotation/tmp/ncRNA_merged.2states.bed.gz",
        tmp_filtered_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/ncRNA_filtered.2states.bed.gz",
        #unfiltered_ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA_unfiltered.2states.sorted.bed.gz",
        filtered_bed_name = "NonCodingRNA{HMM_pass}/annotation/tmp/ncRNA_filtered.2states.sorted_renamed.bed.gz",
        filtered_ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.2states.sorted.bed.gz",
        annotated_temp = "NonCodingRNA{HMM_pass}/annotation/tmp/AnnotatedTemp.bed",
        annotated_HMM = "NonCodingRNA{HMM_pass}/annotation/AnnotatedHMM.bed.gz",
        temp_merged = "NonCodingRNA{HMM_pass}/annotation/tmp/allHMM.merged_temp.bed",
        merged_allHMM = "NonCodingRNA{HMM_pass}/annotation/allHMM.merged.bed.gz",
        protein_coding_HMM = "NonCodingRNA{HMM_pass}/annotation/protein_coding_HMM.bed.gz",
    params:
        #length = 1000,
        length = 1,
        max_overlap = 0.15,
        overlap_flag = '' # "-f 0.1"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/ncRNA_annotation.2states.log"
    resources:
        mem_mb = 58000
    shell:
        """
        (zcat {input.genes_bed} | grep protein_coding - | bedtools intersect {params.overlap_flag} -s -a {input.allHMM_bed} -b - -wa | bedtools sort -i - | bedtools merge -s -i - -d 1000 -c 4,5,6 -o first | gzip - > {output.protein_coding_HMM}) &> {log};
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


##################################################
#    Get expression to filter trailing ncRNAs    #
##################################################

rule MakeSAFAnnotationForNonCodingRNA:
    input:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
    output:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.saf",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        echo -e 'GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output};
        zcat {input} | awk '{{print sep='\\t' $4, $1, $2, $3, $6 }}' FS='\\t' OFS='\\t' >> {output}
        """
        

use rule MakeSAFAnnotationForNonCodingRNA as MakeSAFAnnotationForAllHMM with:
    input:
        "NonCodingRNA{HMM_pass}/annotation/allHMM.merged.bed.gz",
    output:
        "NonCodingRNA{HMM_pass}/annotation/allHMM.merged.saf",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/getSAFforAllHMM.log"


rule featureCountsAllHMM:
    input:
        bam = GetBamForPhenotype,
        saf = "NonCodingRNA{HMM_pass}/annotation/allHMM.merged.saf",
    output:
        "NonCodingRNA{HMM_pass}/Expression_HMM/{Phenotype}_featureCounts/Counts.txt",
    log:
        "logs/NonCodingRNA{HMM_pass}/featureCounts_for_allHMM{Phenotype}.log"
    params:
        extraParams = PairedEndParams,
        strandParams = FeatureCountsNonCodingStrandParams
    threads:
        8
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
        Phenotype = "chRNA.Expression",
    resources:
        mem_mb = 12000,
        cpus_per_node = 9,
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o NonCodingRNA{wildcards.HMM_pass}/Expression_HMM/{wildcards.Phenotype}_featureCounts/Counts.txt {input.bam} &> {log};
        """


rule GetRPKMForAllHMM:
    input:
        "NonCodingRNA{HMM_pass}/Expression_HMM/chRNA.Expression_featureCounts/Counts.txt",
    output:
        "NonCodingRNA{HMM_pass}/Expression_HMM/OnlyFirstReps.RPKM.bed.gz",
        "NonCodingRNA{HMM_pass}/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/PrepareAllHMMRPKM.log"
    conda:
        "../envs/r_essentials.yml"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    resources:
        mem_mb = 32000
    shell:
        """
        Rscript scripts/AllHMM_RPKM.R {wildcards.HMM_pass} &> {log}
        """
        
        
rule FilterNonCodingRNA:
    input:
        "NonCodingRNA{HMM_pass}/annotation/allHMM.merged.bed.gz",
        "NonCodingRNA{HMM_pass}/Expression_HMM/OnlyFirstReps.RPKM.bed.gz",
        "NonCodingRNA{HMM_pass}/Expression_HMM/OnlyFirstReps.qqnorm.bed.gz",
    output:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/filter_ncRNAs.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 32000
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    params:
        merge = GetMergeFlag,
        min_RPKM = -1,
        distance = 20000,
        min_correlation = 0.3, 
        RPKM_ratio = 2
    shell:
        """
        python scripts/filter_ncRNAs.py --min_RPKM {params.min_RPKM} --distance {params.distance} --min_correlation {params.min_correlation} --RPKM_ratio {params.RPKM_ratio} {params.merge} &> {log}
        """
        

### For now, only standard HMM_pass for phenotype.
        
rule featureCountsNonCodingRNA:
    input:
        bam = GetBamForPhenotype,
        saf = "NonCodingRNA_annotation/annotation/ncRNA.saf",
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


########################################
#          Get TSS and 3' end          #
########################################

    
rule SplitAnnotationByStrand:
    input:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz"
    output:
        plus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.minus.bed",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/NonCodingRNA{HMM_pass}/split_annotation.log"
    shell:
        """
        (zcat NonCodingRNA{wildcards.HMM_pass}/annotation/ncRNA.bed.gz | awk -F'\\t' '$6=="+"' -  > {output.plus_bed}) &> {log};
        (zcat NonCodingRNA{wildcards.HMM_pass}/annotation/ncRNA.bed.gz | awk -F'\\t' '$6=="-"' -  > {output.minus_bed}) &>> {log};
        """
        
        
rule GetTSSAnnotations:
    input:
        plus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.minus.bed",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz",
        chrom_sizes = "../data/Chrome.sizes",
        tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        ncRNA_temp = temp("NonCodingRNA{HMM_pass}/annotation/ncRNA.TSS_temp.bed"),
        ncRNA_tss = "NonCodingRNA{HMM_pass}/annotation/ncRNA.TSS.bed.gz",
        genes_temp = temp("NonCodingRNA{HMM_pass}/annotation/allGenes.FirstTSS_temp.bed"),
        genes_tss = "NonCodingRNA{HMM_pass}/annotation/allGenes.FirstTSS.bed.gz",
        #all_genes_tss = "NonCodingRNA{HMM_pass}/annotation/allGenes.TSS.bed.gz"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/TSSAnnotations.log",
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


rule GetAllTSSAnnotations:
    input:
        chrom_sizes = "../data/Chrome.sizes",
        tss_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.tss.bed",
    output:
        all_genes_tss = "NonCodingRNA{HMM_pass}/annotation/allGenes.TSS.bed.gz"
    log:
        "logs/NonCodingRNA{HMM_pass}/AllTSSAnnotations.log",
    params:
        window = 1000
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    resources:
        mem_mb = 58000
    output:
        all_genes_tss = "NonCodingRNA{HMM_pass}/annotation/allGenes.TSS.bed.gz"
    shell:
        """
        (awk '{{print "chr"$1, $2, $3, $7, $6, $4}}' FS='\\t' OFS='\\t' {input.tss_bed} | bedtools slop -b 999 -i - -g {input.chrom_sizes} | gzip - > {output.all_genes_tss}) & >> {log};
        """
       

rule Get3PrimeEndAnnotations:
    input:
        plus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.plus.bed",
        minus_bed = "NonCodingRNA{HMM_pass}/deeptools/bed/ncRNA.minus.bed",
        utr_bed = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.utr.bed",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz",
        chrom_sizes = "../data/Chrome.sizes"
    output:
        ncRNA_temp = temp("NonCodingRNA{HMM_pass}/annotation/ncRNA.3PrimeEnd_temp.bed"),
        genes_temp = temp("NonCodingRNA{HMM_pass}/annotation/allGenes.3PrimeEnd_temp.bed"),
        genes_last_temp = temp("NonCodingRNA{HMM_pass}/annotation/allGenes.Last3PrimeEnd_temp.bed"),
        ncRNA_3prime_end = "NonCodingRNA{HMM_pass}/annotation/ncRNA.3PrimeEnd.bed.gz",
        genes_3prime_end = "NonCodingRNA{HMM_pass}/annotation/allGenes.3PrimeEnd.bed.gz",
        genes_3prime_last = "NonCodingRNA{HMM_pass}/annotation/allGenes.Last3PrimeEnd.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/EndAnnotations.log",
    params:
        window = 200
    resources:
        mem_mb = 58000
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
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
        ncRNA_tss = "NonCodingRNA{HMM_pass}/annotation/ncRNA.TSS.bed.gz",
        genes_tss = "NonCodingRNA{HMM_pass}/annotation/allGenes.TSS.bed.gz",
    resources:
        mem_mb = 58000,
    output:
        uaRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/uaRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/uaRNAs.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.genes_tss} -wo > NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/uaRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_tss} -b {input.ncRNA_tss} -wo >> NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/uaRNA.bed) &>> {log}
        (gzip NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/uaRNA.bed) &>> {log}
        """

rule rtRNA:
    input:
        ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        rtRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/rtRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/rtRNAs.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_bed} -wo > NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/rtRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_bed} -wo >> NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/rtRNA.bed) &>> {log}
        (gzip NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/rtRNA.bed) &>> {log}
        """

rule srtRNA:
    input:
        ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        rtRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/srtRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/srtRNAs.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_bed} -wo -f 1  > NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/srtRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_bed} -wo -f 1 >> NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/srtRNA.bed) &>> {log}
        (gzip NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/srtRNA.bed) &>> {log}
        """

rule coRNA:
    input:
        ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        ncRNA_3prime_end = "NonCodingRNA{HMM_pass}/annotation/ncRNA.3PrimeEnd.bed.gz",
        genes_3prime_last = "NonCodingRNA{HMM_pass}/annotation/allGenes.Last3PrimeEnd.bed.gz",
    resources:
        mem_mb = 58000,
    output:
        coRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/coRNA.bed.gz",
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/coRNAs.log"
    shell:
        """
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.genes_3prime_last} -wo > NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/coRNA.bed) &> {log}
        (bedtools intersect -S -a {input.ncRNA_bed} -b {input.ncRNA_3prime_end} -wo >> NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/coRNA.bed) &>> {log}
        (gzip NonCodingRNA{wildcards.HMM_pass}/annotation/tmp/coRNA.bed) &>> {log}
        """

rule incRNA:
    input:
        ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        incRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/incRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/incRNAs.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (bedtools intersect -v -a {input.ncRNA_bed} -b {input.genes_bed} -wa | gzip - > {output.incRNA}) &> {log}
        """
        
rule ctRNA:
    input:
        ncRNA_bed = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        genes_bed = "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz"
    resources:
        mem_mb = 58000,
    output:
        ctRNA = "NonCodingRNA{HMM_pass}/annotation/tmp/ctRNA.bed.gz",
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/ctRNAs.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (bedtools intersect -s -a {input.ncRNA_bed} -b {input.genes_bed} -wo | gzip - > {output.ctRNA}) &> {log}
        """
        
rule GetLncRNAs:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA{HMM_pass}/annotation/tmp/lncRNA.bed.gz",
        out = "NonCodingRNA{HMM_pass}/annotation/tmp/lncRNA.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/getlncRNA.log"
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' | awk '$7=="lncRNA" {{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
         (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule GetSnoRNAs:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA{HMM_pass}/annotation/tmp/snoRNA.bed.gz",
        out = "NonCodingRNA{HMM_pass}/annotation/tmp/snoRNA.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/getsnoRNA.log"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' | awk '$7=="snoRNA" {{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
         (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule GetPsuedogenes:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
        ncRNA = "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz"
    output:
        bed = "NonCodingRNA{HMM_pass}/annotation/tmp/pseudogenes.bed.gz",
        out = "NonCodingRNA{HMM_pass}/annotation/tmp/pseudogenes.ncRNA.bed.gz"
    resources:
        mem_mb = 24000,
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    log:
        "logs/NonCodingRNA{HMM_pass}/ncRNA_annotation/getPseudogenes.log"
    shell:
        """
        (awk '$3=="gene" {{print $1, $4, $5, $7, $9}}' FS='\\t' OFS='\\t' {input.gtf} | awk '{{print $1, $2, $4, $6}}' FS='"' OFS='\\t' - | grep pseudogene - | awk '{{print $1, $2, $3, $6, $8, $4}}' FS='\\t' OFS='\\t' - | gzip - > {output.bed}) &> {log};
        (bedtools intersect -a {input.ncRNA} -b {output.bed} -F 0.9 -wo | gzip - > {output.out}) &>> {log}
        """
        
rule ClassifyNcRNAs:
    input:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/allGenes.Gencode.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/uaRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/rtRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/srtRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/snoRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/coRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/incRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/ctRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/lncRNA.ncRNA.bed.gz",
        "NonCodingRNA{HMM_pass}/annotation/tmp/pseudogenes.ncRNA.bed.gz",
    output:
        "NonCodingRNA{HMM_pass}/annotation/ncRNA.annotation.tab.gz"
    wildcard_constraints:
        HMM_pass = '_annotation|_merged',
    params:
        GetMergeFlag
    log:
        "logs/NonCodingRNA{HMM_pass}/classify_ncRNA.log"
    resources:
        mem_mb = 58000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/classify_ncRNA.py {params} &> {log}
        """


###################################
#   Overlap with histone marks    #
###################################

def GetTSSFile(wildcards):
    if wildcards.Annotation == 'ncRNA':
        return "NonCodingRNA_annotation/annotation/ncRNA.TSS.bed.gz"
    elif wildcards.Annotation == 'allGenes':
        return "NonCodingRNA_annotation/annotation/allGenes.TSS.bed.gz"
        
def GetAnnotationFile(wildcards):
    if wildcards.Annotation == 'ncRNA':
        return "NonCodingRNA_annotation/annotation/ncRNA.bed.gz"
    elif wildcards.Annotation == 'allGenes':
        return "NonCodingRNA_annotation/annotation/allGenes.bed.gz"

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
        'NonCodingRNA_annotation/annotation/histone_marks/{Annotation}.{Histone}.TSS_overlaps.bed.gz'
    log:
        'logs/NonCodingRNA_annotation/{Annotation}.{Histone}.TSS_overlaps.log'
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
        'NonCodingRNA_annotation/annotation/histone_marks/{Annotation}.{Histone}.overlaps.bed.gz'
    log:
        'logs/NonCodingRNA_annotation/{Annotation}.{Histone}.overlaps.log'
        
rule HistoneAnnotation:
    input:
        'NonCodingRNA_annotation/annotation/ncRNA.annotation.tab.gz',
        'NonCodingRNA_annotation/annotation/allGenes.TSS.bed.gz',
        expand("NonCodingRNA_annotation/annotation/histone_marks/{Annotation}.{Histone}.{overlap_type}overlaps.bed.gz",
            Annotation = ["ncRNA", "allGenes"], Histone = ["H3K4ME1", "H3K4ME3", "H3K27AC"], overlap_type=["", "TSS_"]
        )
    output:
        'NonCodingRNA_annotation/annotation/ncRNA.histone.tab.gz',
        'NonCodingRNA_annotation/annotation/allGenes.histone.tab.gz'
    log:
        'logs/NonCodingRNA_annotation/histone_annotation.log'
    resources:
        mem_mb = 12000
    shell:
        """
        python scripts/make_histone_ncRNA_annotations.py &> {log}
        """
        
        
###################################
#   Prepare phenotypes for QC     #
###################################


rule Prepare_polyA_ExpressionPhenotypes:
    input:
        "featureCounts/polyA.Expression/Counts.txt",
        "featureCounts/polyA.Expression_ncRNA/Counts.txt",
        "featureCounts/polyA.Expression_lncRNA/Counts.txt",
        "featureCounts/polyA.Expression_snoRNA/Counts.txt",
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "NonCodingRNA_annotation/annotation/ncRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.RPKM.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_polyA_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/prepare_polyA_ncRNA_Phenotypes.R &> {log}
        """
        
rule Prepare_polyA_Subset_YRI_ExpressionPhenotypes:
    input:
        "QTLs/QTLTools/polyA.Expression_ncRNA/OnlyFirstReps.RPKM.bed.gz",
        "QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz",
    output:
        "QTLs/QTLTools/polyA.Expression_ncRNA.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_polyA_Expression_Phenotypes.Subset_YRI.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/prepare_polyA_SubsetYRI_ncRNA_Phenotypes.R &> {log}
        """
        
rule Prepare_ml30_ExpressionPhenotypes:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "featureCounts/MetabolicLabelled.30min/Counts.txt",
        "featureCounts/MetabolicLabelled.30min_ncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.30min_lncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.30min_snoRNA/Counts.txt",
        "NonCodingRNA_annotation/annotation/ncRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/MetabolicLabelled.30min_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/MetabolicLabelled.30min_ncRNA/OnlyFirstReps.RPKM.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_ml30_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/prepare_ml30_ncRNA_Phenotypes.R &> {log}
        """
        
        
rule Prepare_ml60_ExpressionPhenotypes:
    input:
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt",
        "featureCounts/MetabolicLabelled.60min/Counts.txt",
        "featureCounts/MetabolicLabelled.60min_ncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.60min_lncRNA/Counts.txt",
        "featureCounts/MetabolicLabelled.60min_snoRNA/Counts.txt",
        "NonCodingRNA_annotation/annotation/ncRNA.annotation.tab.gz"
    output:
        "QTLs/QTLTools/MetabolicLabelled.60min_ncRNA/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/MetabolicLabelled.60min_ncRNA/OnlyFirstReps.RPKM.bed.gz",
    log:
        "logs/NonCodingRNA/Prepare_ml60_Expression_Phenotypes.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/prepare_ml60_ncRNA_Phenotypes.R &> {log}
        """
        







