localrules: DownloadFasterqDump_SE

rule DownloadFasterqDump_SE:
    """
    FasterqDump for single end read datasets. can't be parallelized on midway
    compute nodes. is painfully slow relative to aspera. Using aspera links are
    better...
    This rule is probably not needed because aspera from ENA (which is mostly snychronized with SRA) is faster and more reliable. There is a different rule for using aspera to download fastq files
    """
    output:
        "FastqSE/{Phenotype}/{IndID}/{Rep}.SE.fastq.gz"
    log:
        "logs/DownloadFasterqDump_SE/{Phenotype}/{IndID}.{Rep}.log"
    params:
        SRA_Run = GetDownloadLinkFuncs('SRA_Run'),
    shadow: "shallow"
    threads: 1
    resources:
        tasks_per_node = 7
    shell:
        """
        fasterq-dump -S -v -e {threads} {params.SRA_Run} &> {log}
        for accession in {params.SRA_Run};
        do
            gzip -c $accession.fastq >> {output}
            rm $accession.fastq
        done
        """
