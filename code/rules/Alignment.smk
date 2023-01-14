
rule MakeHisat2Index:
    input:
        "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2"
    log:
        "logs/MakeIndex.log"
    conda:
        "../envs/hisat2.yml"
    threads: 8
    resources:
        mem_mb = 42000
    shell:
        """
        hisat2-build {input} {input} -p {threads} &> {log}
        """

rule gtf2leafcutter:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        multiext("ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf", "_all_exons.txt.gz", "_all_introns.bed.gz", "_fiveprime.bed.gz", "_threeprime.bed.gz")
    log:
        "logs/gtf2leafcutter.log"
    shell:
        """
        scripts/leafcutter/leafviz/gtf2leafcutter.pl -o {input.gtf} {input.gtf} &> {log}
        """

def GetGTFToolsAnnotation(wildcards):
    if wildcards.GTFTools == "GTFTools":
        return "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    elif wildcards.GTFTools == "GTFTools_BasicAnnotations":
        return "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf"

rule GTFTools:
    input:
        GetGTFToolsAnnotation
    output:
        tss = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.tss.bed",
        exons = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.exons.bed",
        introns  = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.introns.bed",
        genes = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.genes.bed",
        splicesite  = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.splicesite.bed",
        utr  = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.utr.bed",
        maskedint  = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.maskedint.bed",
        independentintron  = "ReferenceGenome/Annotations/{GTFTools}/gencode.v34.chromasomal.independentintron.bed",
    wildcard_constraints:
        GTFTools = "GTFTools|GTFTools_BasicAnnotations"
    log:
        "logs/{GTFTools}.log"
    shell:
        """
        python2.7 scripts/GTFtools_0.8.0/gtftools.py -e {output.exons} -i {output.introns} -d {output.independentintron} -k {output.maskedint} --splice_site {output.splicesite} -u {output.utr} -g {output.genes} -t {output.tss} -w 1 {input} 2> {log}
        """

rule STAR_make_index:
    """
    did not work on bigmem2. Never figured out why (the log file didn't
    indicate anything). Ran on login node with success.
    """
    input:
        fasta = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
    log:
        "logs/STAR_make_index.log"
    params:
        genomeDir = "ReferenceGenome/STARIndex/"
    threads: 4
    resources:
        mem = "42G",
        partition = "bigmem2",
        ntasks = 5
    shell:
        """
        module load STAR/2.7.7a
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN {threads} --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule MakeDummyVcfForSTAR_WASP_MODE:
    """
    https://github.com/alexdobin/STAR/issues/772
    The vcf for STAR WASP mode should contain genotypes for a single
    (personalized) sample, or contain genotypes for a dummy sample with
    heterozygous genotypes at all variant positions.
    """
    input:
        "Genotypes/1KG_GRCh38/{chrom}.vcf.gz"
    params:
        IndIDs = lambda wildcards: ",".join(Fastq_samples.loc[(Fastq_samples['Include']==True) & (Fastq_samples['Phenotype']==wildcards.Phenotype), ['IndID']].drop_duplicates()['IndID'])
    log:
        "logs/MakeDummyVcfForSTAR_WASP_MODE/{Phenotype}/{chrom}.log"
    output:
        vcf = temp("ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/{chrom}.vcf.gz"),
        tbi = temp("ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/{chrom}.vcf.gz.tbi")
    shell:
        """
        (bcftools view -s {params.IndIDs} --force-samples -c 1:minor {input} | awk -F'\\t' -v OFS='\\t' '$1~"^##" {{print}} $1=="#CHROM" {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,"DummyID"}} $1!~"^#" {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0|1"}}' | bgzip -c /dev/stdin > {output}) &> {log}
        tabix -p vcf {output.vcf}
        """

rule MergeVcfForSTAR_WASP_MODE:
    input:
        vcf = expand("ReferenceGenome/STAR_WASP_Vcfs/{{Phenotype}}/{chrom}.vcf.gz", chrom=autosomes),
        tbi = expand("ReferenceGenome/STAR_WASP_Vcfs/{{Phenotype}}/{chrom}.vcf.gz.tbi", chrom=autosomes)
    output:
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf",
    shell:
        """
        bcftools concat -O v {input.vcf} > {output.vcf}
        """

rule STAR_Align_WASP:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz",
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf"
    output:
        bam = temp("Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam"),
        align_log = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align_WASP/{Phenotype}/{IndID}.{Rep}.log"
    params:
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        WASP_params = "--waspOutputMode SAMtag --outSAMattributes NH HI AS nM XS vW --varVCFfile",
        JunctionScore = GetSTARJunctionScoreParams
    resources:
        cpus_per_node = 9,
        mem = 58000,
    wildcard_constraints:
        Phenotype = "Expression.Splicing|chRNA.Expression.Splicing"
    shell:
        """
        module load STAR/2.7.7a
        STAR --readMapNumber {params.readMapNumber} {params.JunctionScore} --outFileNamePrefix Alignments/STAR_Align/{wildcards.Phenotype}/{wildcards.IndID}/{wildcards.Rep}/ --genomeDir ReferenceGenome/STARIndex/ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 {params.WASP_params} {input.vcf} --limitBAMsortRAM 16000000000 {params.ENCODE_params} --outSAMstrandField intronMotif  &> {log}
        """

use rule STAR_Align_WASP as STAR_Align_WASP_SE with:
    input:
        index = "ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "FastqFastpSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz",
        R2 = [],
        vcf = "ReferenceGenome/STAR_WASP_Vcfs/{Phenotype}/WholeGenome.vcf"
    wildcard_constraints:
        Phenotype = "MetabolicLabelled.30min|MetabolicLabelled.60min|ProCap"

rule FilterBAM_WaspTags:
    input:
        bam = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Aligned.sortedByCoord.out.bam"
    output:
        bam = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam",
        bai = "Alignments/STAR_Align/{Phenotype}/{IndID}/{Rep}/Filtered.bam.bai",
    log:
        "logs/FilterBAM_WaspTags/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 8000
    shell:
        """
        (cat <(samtools view -H {input.bam}) <(samtools view {input.bam} | grep -v "vW:i:[2-7]") | samtools view -bh > {output.bam}) &> {log}
        samtools index {output.bam} &>> {log}
        """

rule Hisat2_Align:
    input:
        Ref = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2",
        R1 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R1.fastq.gz",
        R2 = "FastqFastp/{Phenotype}/{IndID}/{Rep}.R2.fastq.gz",
    output:
        bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam"),
        bai = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai")
    threads: 8
    conda:
        "../envs/hisat2.yml"
    params:
        MaxPE_InsertLen = 1000,
    wildcard_constraints:
        Phenotype = "H3K4ME3|H3K4ME1|H3K27AC"
    log:
        "logs/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 58000,
        cpus_per_node = 9
    shell:
        """
        (hisat2 -1 {input.R1} -2 {input.R2} -x ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa -p {threads} --no-spliced-alignment --no-discordant --maxins {params.MaxPE_InsertLen} | samtools view -bh -F256 | samtools sort > {output.bam} ) &> {log}
        samtools index {output.bam}
        """
        
def GetDNaseIReplicates(wildcards):
    replicates = [str(x) for x in list(
    DNaseSamples_df.loc[DNaseSamples_df.IndID == wildcards.IndID].RepNumber
    )]
    return expand("FastqFastpSE/DNaseISensitivity/{{IndID}}/{Rep}.SE.fastq.gz", Rep=replicates)
    
    
rule MergeRemapFastqForDNaseISensitivity:
    input:
        GetDNaseIReplicates
    output:
        temp("FastqFastpSE/DNaseISensitivity/{IndID}/merged.SE_.fastq.gz")
    log:
        "logs/MergeHornetFastq/{IndID}/log"
    shell:
        """
        (zcat FastqFastpSE/DNaseISensitivity/{wildcards.IndID}/*.SE.fastq.gz | gzip - > {output}) &> {log}
        """

rule Hisat2_Align_SE:
    input:
        Ref = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2",
        R1 = "FastqFastpSE/{Phenotype}/{IndID}/{Rep}.SE_.fastq.gz",
    output:
        bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam"),
        bai = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai")
    wildcard_constraints:
        Phenotype = "DNaseISensitivity",
        Rep='merged'
    threads: 8
    conda:
        "../envs/hisat2.yml"
    params:
        MaxPE_InsertLen = 1000,
    log:
        "logs/Hisat2_Align_SE/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 58000,
        cpus_per_node = 9
    shell:
        """
        (hisat2 -U {input.R1} -x ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa -p {threads} --no-spliced-alignment --no-discordant --maxins {params.MaxPE_InsertLen} | samtools view -bh -F256 | samtools sort > {output.bam} ) &> {log}
        samtools index {output.bam}
        """


# rule FilterHisat2PE_dups:
#     input:
#         bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam",
#         bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai"
#     output:
#         bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.dedup.bam",
#         bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.dedup.bam.bai"
#     log:
#         "logs/FilterHisat2PE_dups/{Phenotype}/{IndID}.{Rep}.log"
#     shell:
#         """
#         samtools view -bh -F 3072 {input.bam} > {output.bam} 2> {log}
#         samtools index {output.bam}
#         """

rule MakeHornetSnpList:
    input:
        "Genotypes/1KG_GRCh38/{chrom}.vcf.gz"
    params:
        IndIDs = lambda wildcards: ",".join(Fastq_samples.loc[(Fastq_samples['Include']==True) & (Fastq_samples['Phenotype']==wildcards.Phenotype), ['IndID']].drop_duplicates()['IndID'])
    log:
        "logs/MakeHornetSnpList/{Phenotype}/{chrom}.log"
    output:
        snps = temp("ReferenceGenome/HornetSnpLists/{Phenotype}/chr{chrom}.snps.txt.gz"),
    shell:
        """
        (bcftools view -H -s {params.IndIDs} -m2 -M2 -v snps --force-samples -c 1:minor {input} | awk -F'\\t' '{{ print $2, $4, $5 }}' | gzip - > {output} ) &> {log}
        """

rule Hornet_find_intersecting_snps:
    """
    keep: reads to keep that do not intersect snps
    remap: reads to remap. In fastq, the snp allele is switched
    """
    input:
        bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam",
        bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai",
        snp_lists = expand("ReferenceGenome/HornetSnpLists/{{Phenotype}}/chr{chrom}.snps.txt.gz", chrom=autosomes)
    output:
        R1 = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq1.gz"),
        R2 = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq2.gz"),
        keep = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.keep.bam"),
        remap = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.to.remap.bam"),
    log:
        "logs/Hornet_find_intersecting_snps/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 12000
    wildcard_constraints:
        Phenotype = "H3K4ME3|H3K4ME1|H3K27AC"
    shell:
        """
        python scripts/Hornet/mapping/find_intersecting_snps.py -p {input.bam} ReferenceGenome/HornetSnpLists/{wildcards.Phenotype} &> {log}
        """
        
rule Hornet_find_intersecting_snps_SE:
    input:
        bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam",
        bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.bam.bai",
        snp_lists = expand("ReferenceGenome/HornetSnpLists/{{Phenotype}}/chr{chrom}.snps.txt.gz", chrom=autosomes)
    output:
        R1 = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq.gz"),
        keep = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.keep.bam"),
        remap = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.to.remap.bam"),
    wildcard_constraints:
        Phenotype = "DNaseISensitivity",
        Rep='merged'
    log:
        "logs/Hornet_find_intersecting_snps_SE/{Phenotype}/{IndID}.{Rep}.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem = 12000
    shell:
        """
        python scripts/Hornet/mapping/find_intersecting_snps.py {input.bam} ReferenceGenome/HornetSnpLists/{wildcards.Phenotype} &> {log}
        """

use rule Hisat2_Align as Hornet_Hisat2_Align with:
    input:
        Ref = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2",
        R1 = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq1.gz",
        R2 = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq2.gz",
    output:
        bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.bam"),
        bai = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.bam.bai")
    wildcard_constraints:
        Phenotype = "H3K4ME3|H3K4ME1|H3K27AC"
    log:
        "logs/Hornet_Hisat2_Align/{Phenotype}/{IndID}.{Rep}.log"
        
        

        
use rule Hisat2_Align_SE as Hornet_Hisat2_Align_SE with:
    input:
        Ref = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        ht1 = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.1.ht2",
        R1 = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.remap.fq.gz",
    wildcard_constraints:
        Phenotype = "DNaseISensitivity",
        Rep="merged"
    log:
        "logs/Hornet_Hisat2_Align_SE/{Phenotype}/{IndID}.{Rep}.log"
    output:
        bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.bam"),
        bai = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.bam.bai")

rule Hornet_filter_remapped:
    """
    check that remapped reads with allele siwtched map to same place. If they
    do, write to kept. Merge kept with keep from Hornet_find_intersecting_snps
    rule
    """
    input:
        remap_bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.bam",
        keep_bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.keep.bam",
        to_remap_bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.to.remap.bam",
    output:
        remap_bam_name_sorted = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.namesorted.bam"),
        kept_bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.hornet.remapped.kept.bam"),
        keep_bam_sorted = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.keep.sorted.bam"),
        kept_bam_sorted = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.kept.sorted.bam"),
        merged_bam = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.sorted.bam"),
        merged_bai = temp("Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.sorted.bam.bai"),
    log:
        "logs/Hornet_filter_remapped/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 16000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        (samtools sort -n {input.remap_bam} > {output.remap_bam_name_sorted}) &> {log}
        python scripts/Hornet/mapping/filter_remapped_reads.py {input.to_remap_bam} {output.remap_bam_name_sorted} {output.kept_bam} &>> {log}
        (samtools sort {output.kept_bam} > {output.kept_bam_sorted}) &>> {log}
        (samtools sort {input.keep_bam} > {output.keep_bam_sorted}) &>> {log}
        samtools merge {output.merged_bam} {output.kept_bam_sorted} {output.keep_bam_sorted} &>> {log}
        samtools index {output.merged_bam} &>> {log}
        """

use rule Hornet_filter_remapped as Hornet_filter_remapped_DNase with:
    wildcard_constraints:
        Rep="merged",
        Phenotype = "DNaseISensitivity"

rule MarkDups:
    input:
        bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.sorted.bam",
        bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.sorted.bam.bai",
    output:
        bam = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam",
        bai = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam.bai",
    log:
        "logs/MarkDups/{Phenotype}/{IndID}.{Rep}.log"
    shell:
        """
        (samtools sort -n {input.bam} | samtools fixmate -m - - |  samtools sort | samtools markdup - -  > {output.bam} ) &> {log}
        samtools index {output.bam}
        """
        
use rule MarkDups as MarkDupsDNase with:
    wildcard_constraints:
        Rep="merged",
        Phenotype = "DNaseISensitivity"


