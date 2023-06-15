rule CollectLR_Bed12:
    input:
        expand("LongReads/TerminiBeds5/{sample}.bed.gz", sample=long_read_samples),
        expand("LongReads/ClosestAnnotatedTermini/{GTFType}_{TES_or_TSS}/{sample}.bed.gz", GTFType = ["GTFTools", "GTFTools_BasicAnnotations"], TES_or_TSS = ["TSS", "TES"], sample=long_read_samples)

rule Bam2Bed12:
    input:
        bam = "/project2/yangili1/cfbuenabadn/iso-seq/demux.5p--YG_{sample}_3p.bam.fastq.hg38.sort.bam",
        bai = "/project2/yangili1/cfbuenabadn/iso-seq/demux.5p--YG_{sample}_3p.bam.fastq.hg38.sort.bam.bai",
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        bed = "LongReads/bed12/{sample}.bed.gz",
        tbi = "LongReads/bed12/{sample}.bed.gz.tbi",
    log:
        "logs/LongReads/Bam2Bed12.{sample}.log"
    wildcard_constraints:
        sample = '|'.join(long_read_samples) 
    shell:
        """
        (bedtools bamtobed -i {input.bam} -bed12 -split | awk -F'\\t' -v OFS='\\t' '$11 !~ ",0,"' | bedtools getfasta -bedOut -s -split -bed - -name -fi {input.fa} | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        (tabix -p bed {output.bed}) &>> {log}
        """
        
use rule Bam2Bed12 as Bam2Bed12_ONT with:
    input:
        bam = "Alignments/minimap2_Alignment/{sample}.sort.bam",
        bai = "Alignments/minimap2_Alignment/{sample}.sort.bam.bai",
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        bed = "LongReads/bed12/{sample}.bed.gz",
        tbi = "LongReads/bed12/{sample}.bed.gz.tbi",
    wildcard_constraints:
        sample = '|'.join(long_read_samples_ONT) 

rule Bam2Bed12_5_and_3Termini:
    input:
        bed = "LongReads/bed12/{sample}.bed.gz",
    output:
        bed5 = "LongReads/TerminiBeds5/{sample}.bed.gz",
        bed3 = "LongReads/TerminiBeds3/{sample}.bed.gz"
    shell:
        """
        zcat {input} | awk -F'\\t' -v OFS='\\t' '$6 == "+" {{ print $1, $2, $2+1, $4, $5, $6 }} $6 == "-" {{ print $1, $3-1, $3, $4, $5, $6 }}' | bedtools sort -i - | gzip - > {output.bed5} 
        zcat {input} | awk -F'\\t' -v OFS='\\t' '$6 == "+" {{ print $1, $3-1, $3, $4, $5, $6 }} $6 == "-" {{ print $1, $2, $2+1, $4, $5, $6 }}' | bedtools sort -i - | gzip - > {output.bed3} 
        """

rule GetAnnotated_TES_TSS_Site:
    input:
        "ReferenceGenome/Annotations/{GTFType}/gencode.v34.chromasomal.utr.bed"
    output:
        bed5= "LongReads/ReferenceFeatures/{GTFType}/gencode.v34.chromasomal.SingleNt_TSS.bed",
        bed3="LongReads/ReferenceFeatures/{GTFType}/gencode.v34.chromasomal.SingleNt_TES.bed"
    wildcard_constraints:
        GTFType = "GTFTools_BasicAnnotations|GTFTools"
    shell:
        """
        cat {input} | awk -F'\\t' -v OFS='\\t' '$6=="5UTR" {{print "chr"$1, $2, $3, $5, ".", $4}}' | bedtools sort -i - | uniq | awk -F'\\t' -v OFS='\\t' '$6 == "+" {{ print $1, $2, $2+1, $4, $5, $6 }} $6 == "-" {{ print $1, $3-1, $3, $4, $5, $6 }}' | bedtools sort -i - | gzip - > {output.bed5}
        cat {input} | awk -F'\\t' -v OFS='\\t' '$6=="3UTR" {{print "chr"$1, $2, $3, $5, ".", $4}}' | bedtools sort -i - | uniq | awk -F'\\t' -v OFS='\\t' '$6 == "+" {{ print $1, $3-1, $3, $4, $5, $6 }} $6 == "-" {{ print $1, $2, $2+1, $4, $5, $6 }}' | bedtools sort -i - | gzip - > {output.bed3}
        """

def GetReadsTermini(wildcards):
    if wildcards.TES_or_TSS == "TSS":
        return "LongReads/TerminiBeds5/{sample}.bed.gz" 
    elif wildcards.TES_or_TSS == "TES":
        return"LongReads/TerminiBeds3/{sample}.bed.gz" 

rule ClosestTES_TSS:
    input:
        bed_annotated = "LongReads/ReferenceFeatures/{GTFType}/gencode.v34.chromasomal.SingleNt_{TES_or_TSS}.bed",
        bed_reads_termini = GetReadsTermini
    wildcard_constraints:
        GTFType = "GTFTools_BasicAnnotations|GTFTools",
        TES_or_TSS = "TSS|TES"
    output:
        "LongReads/ClosestAnnotatedTermini/{GTFType}_{TES_or_TSS}/{sample}.bed.gz"
    shell:
        """
        bedtools closest -s -D a -t first -a {input.bed_reads_termini} -b {input.bed_annotated} | gzip - > {output}
        """
