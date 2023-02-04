rule Bam2JuncPerRead:
    input:
        bam = "/project2/yangili1/cfbuenabadn/iso-seq/demux.5p--YG_{sample}_3p.bam.fastq.hg38.sort.bam",
        bai = "/project2/yangili1/cfbuenabadn/iso-seq/demux.5p--YG_{sample}_3p.bam.fastq.hg38.sort.bam.bai",
    output:
        junc = "LongReads/Junctions/{sample}.junc.gz",
    log:
        "logs/LongReads/Junctions/Bam2JuncPerRead.{sample}.log"
    conda:
        "../envs/py_tools.yml"
    wildcard_constraints:
        sample = '|'.join(long_read_samples) 
    shell:
        """
        python scripts/Bam2JuncPerRead.py --input {input.bam} --output {output.junc} &> {log};
        """
        
use rule Bam2JuncPerRead as Bam2JuncPerRead_ONT with:
    input:
        bam = "Alignments/minimap2_Alignment/{sample}.sort.bam",
        bai = "Alignments/minimap2_Alignment/{sample}.sort.bam.bai",
    wildcard_constraints:
        sample = '|'.join(long_read_samples_ONT) 
        
rule GetSpliceJunctionAnnotation:
    input:
        "/project2/yangili1/yangili/chRNA/annotation_leaf_JAN28.txt.gz",
    output:
        "LongReads/Annotations/Annotation.bed.gz"
    log:
        "logs/LongReads/Annotations/Annotation.log"
    shell:
        """
        (cat {input} | awk '{{print $1, $2, $3, $6, $7, $4, $8}}' FS='\\t' OFS='\\t' - | sort -u - | bedtools sort -i - | gzip - > {output}) &> {log}
        """
        
rule AnnotateJuncFiles:
    input:
        junc = "LongReads/Junctions/{sample}.junc.gz",
        annot = "LongReads/Annotations/Annotation.bed.gz"
    output:
        "LongReads/Junctions/{sample}.annotated.junc.gz"
    wildcard_constraints:
        sample = '|'.join(all_long_reads) 
    log:
        "logs/LongReads/Junctions/Annotate.{sample}.log"
    shell:
        """
        (bedtools intersect -r -f 1 -a {input.junc} -b {input.annot} -wb | awk '{{print $1, $2, $3, $1":"$2"-"$3":"$6, $5, $6, $10, $11, $13}}' FS='\\t' OFS='\\t' - | gzip - > {output}) &> {log}
        """
        
def GetDownloadLinkFuncsONT(LinkType):
    def F(wildcards):
        df_subset = long_read_samples_df.loc[
                (long_read_samples_df['IndID'] == wildcards.IndID) &
                (long_read_samples_df['Phenotype'] == wildcards.Phenotype) 
                ]
        return df_subset['R1_' + LinkType].fillna('').tolist()
    return F
        
rule DownloadONTFastqFromLink:
    """
    Download w aspera if possible, then try ftp link, then try fasterq-dump
    """
    input:
        aspera_key = config['aspera_key']
    output:
        fastq = "Fastq/ONT/{Phenotype}/{IndID}/R1.fastq.gz",
    log:
        "logs/DownloadONTFastqFromAsperaLink/{Phenotype}.{IndID}.R1.log"
    shadow: "shallow"
    params:
        ftp_link = GetDownloadLinkFuncsONT('ftp'),
        aspera_link  = GetDownloadLinkFuncsONT('aspera'),
    wildcard_constraints:
        Phenotype = '|'.join(long_read_samples_df.Phenotype.unique()),
        IndID = '|'.join(long_read_samples_df.IndID.unique())
    shell:
        """
        #Use aspera if aspera key file and aspera links parameters are defined (not empty strings)
        if [[ ! -z "{input.aspera_key}" && ! -z "{params.aspera_link}" ]]; then
            for link in {params.aspera_link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                ascp -v -QT -l 300m -P33001 -i {input.aspera_key} era-fasp@${{link}} $tmpfile &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        else
            for link in {params.ftp_link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                wget -O $tmpfile ${{link}} &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        fi
        """
        
rule CreateMinimapIndex:
    input:
        "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa"
    output:
        "ReferenceGenome/minimap2_index/minimap2.mmi"
    log:
        "logs/minimap2/minimap2.index.log"
    resources:
        mem_mb = 36000
    shell:
        """
        (minimap2 -d {output} {input}) &> {log}
        """

rule RunMinimap2:
    input:
        #index = "ReferenceGenome/minimap2_index/minimap2.mmi",
        ref = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        fastq = "Fastq/ONT/{Phenotype}/{IndID}/R1.fastq.gz"
    output:
        sam = temp("Alignments/minimap2_Alignment/{Phenotype}.{IndID}.sam"),
        bam = temp("Alignments/minimap2_Alignment/{Phenotype}.{IndID}.bam"),
    log:
        "logs/minimap2/minimap2.alignment.{Phenotype}.{IndID}.log"
    wildcard_constraints:
        Phenotype = '|'.join(list(long_read_samples_df.Phenotype)),
        IndID = '|'.join(list(long_read_samples_df.IndID))
    resources:
        mem_mb = 36000
    shell:
        """
        (minimap2 -ax splice {input.ref} {input.fastq} > {output.sam}) &> {log};
        (samtools view -S -b {output.sam} > {output.bam}) &>> {log}
        """

rule minimap2SortAndIndexBam:
    input:
        "Alignments/minimap2_Alignment/{sample}.bam",
    output:
        bam = "Alignments/minimap2_Alignment/{sample}.sort.bam",
        bai = "Alignments/minimap2_Alignment/{sample}.sort.bam.bai"
    log:
        "logs/minimap2/{sample}.index.log"
    wildcard_constraints:
        sample = "|".join(long_read_samples_ONT)
    shell:
        """
        samtools sort -o {output.bam} {input};
        samtools index {output.bam} &> {log}
        """
        
rule GetNMDPerJunctionsAndSampledAvgs:
    input:
        expand("LongReads/Junctions/{sample}.annotated.junc.gz", sample = long_read_samples)
    output:
        'LongReads/Analysis/IsoSeq.nmd.tab.gz',
        'LongReads/Analysis/IsoSeq.stable.tab.gz',
        'LongReads/Analysis/IsoSeq.nmd_avg.tab.gz',
        'LongReads/Analysis/IsoSeq.stable_avg.tab.gz',
        'LongReads/Analysis/NMD_KD.nmd.tab.gz',
        'LongReads/Analysis/NMD_KD.stable.tab.gz',
        'LongReads/Analysis/NMD_KD.nmd_avg.tab.gz',
        'LongReads/Analysis/NMD_KD.stable_avg.tab.gz',
        'LongReads/Analysis/Churchman.nmd.tab.gz',
        'LongReads/Analysis/Churchman.stable.tab.gz',
        'LongReads/Analysis/Churchman.nmd_avg.tab.gz',
        'LongReads/Analysis/Churchman.stable_avg.tab.gz',
    log:
        'logs/LongReads/subsample.log'
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 36000
    shell:
        """
        (python scripts/sample_long_read_NMD_junctions.py) &> {log}
        """
    
rule GetNMDPerJunctionByQuartile:
    input:
        expand("LongReads/Junctions/{sample}.annotated.junc.gz", sample = long_read_samples),
        'QTLs/QTLTools/Expression.Splicing.Subset_YRI/OnlyFirstRepsUnstandardized.qqnorm.bed.gz',
    output:
        'LongReads/Analysis/IsoSeq.ByQuartile.tab.gz',
        'LongReads/Analysis/NMD_KD.ByQuartile.tab.gz',
        'LongReads/Analysis/Churchman.ByQuartile.tab.gz'
    log:
        'logs/LongReads/analysis_by_quartile.log'
    resources:
        mem_mb = 36000
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        (python scripts/sample_long_read_by_expression_quartiles.py) &> {log}
        """
        