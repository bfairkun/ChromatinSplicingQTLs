rule Bam2JuncPerRead:
    input:
        bam = "/project2/yangili1/cfbuenabadn/iso-seq/demux.5p--YG_{sample}_3p.bam.fastq.hg38.sort.bam",
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
        
rule GetSpliceJunctionAnnotation:
    input:
        "/project2/yangili1/yangili/chRNA/annotation_leaf_JAN19.txt.gz",
    output:
        "LongReads/Annotations/Annotation.bed.gz"
    log:
        "logs/LongReads/Annotations/Annotation.log"
    shell:
        """
        (cat {input} | awk '{{print $1, $2, $3, $6, $6, $4}}' FS='\\t' OFS='\\t' - | sort -u - | bedtools sort -i - | gzip - > {output}) &> {log}
        """
        
rule AnnotateJuncFiles:
    input:
        junc = "LongReads/Junctions/{sample}.junc.gz",
        annot = "LongReads/Annotations/Annotation.bed.gz"
    output:
        "LongReads/Junctions/{sample}.annotated.junc.gz"
    log:
        "logs/LongReads/Junctions/Annotate.{sample}.log"
    shell:
        """
        (bedtools intersect -r -f 1 -a {input.junc} -b {input.annot} -wb | awk '{{print $1, $2, $3, $1":"$2"-"$3":"$6, $5, $6, $10}}' FS='\\t' OFS='\\t' - | gzip - > {output}) &> {log}
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
        
#cat /project2/yangili1/yangili/chRNA/annotation_leaf_JAN19.txt.gz | awk '{print $1, $2, $3, $6, $6, $4}' #FS='\t' OFS='\t' - | sort -u | bedtools sort -i - | gzip - > YangAnnotation.bed.gz        
        
#                annot = "SplicingAnalysis/IntronTypeAnnotations.JoinedWithLeafcutterJuncList.txt.gz"
#                junc_filtered = "LongReads/Junctions/{sample}.junc.filtered.bed.gz",
#        """
#        (bedtools intersect -s -r -f 1 -a {output.junc} -b {input.annot} -wb | gzip - > {output.junc_filtered}) &>> {log}
#        """
        