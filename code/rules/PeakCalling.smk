rule Macs2PeakCalling_narrow:
    input:
        bams = GetBamForPhenotype
    output:
        peaks = "PeakCalling/{Phenotype}_peaks.narrowPeak",
        saf = "PeakCalling/{Phenotype}_peaks.narrowPeak.saf"
    log:
        "logs/Macs2PeakCalling_narrow/{Phenotype}.log"
    resources:
        mem = 58000
    params:
        "--tempdir /scratch/midway2/cnajar/ --outdir PeakCalling/ --name {Phenotype}"
    shell:
        """
        macs2 callpeak {params} -f BAMPE --name {wildcards.Phenotype} -t {input.bams} &> {log}
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }} {{ print $4, $1, $2, $3, $6 }}' {output.peaks} > {output.saf}
        """

use rule Macs2PeakCalling_narrow as Macs2PeakCalling_broad with:
    output:
        peaks = "PeakCalling/{Phenotype}_peaks.broadPeak",
        saf = "PeakCalling/{Phenotype}_peaks.broadPeak.saf"
    params:
        "--tempdir /scratch/midway2/cnajar/ --outdir PeakCalling/ --name {Phenotype} --broad"
    log:
        "logs/Macs2PeakCalling_broad/{Phenotype}.log"

use rule Macs2PeakCalling_narrow as Macs2PeakCalling_narrow_perind with:
    input:
        bams = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam"
    output:
        peaks = "PeakCallingPerInd/{Phenotype}/{IndID}.{Rep}_peaks.narrowPeak",
        saf = "PeakCallingPerInd/{Phenotype}/{IndID}.{Rep}_peaks.narrowPeak.saf"
    params:
        "--tempdir /scratch/midway2/cnajar/ --outdir PeakCallingPerInd/{Phenotype}/ --name {IndID}.{Rep}"
    log:
        "logs/Macs2PeakCalling_narrow_perind/{Phenotype}/{IndID}.{Rep}.log"

use rule Macs2PeakCalling_narrow as Macs2PeakCalling_broad_perind with:
    input:
        bams = "Alignments/Hisat2_Align/{Phenotype}/{IndID}.{Rep}.wasp_filterd.markdup.sorted.bam"
    output:
        peaks = "PeakCallingPerInd/{Phenotype}/{IndID}.{Rep}_peaks.broadPeak",
        saf = "PeakCallingPerInd/{Phenotype}/{IndID}.{Rep}_peaks.broadPeak.saf"
    params:
        "--tempdir /scratch/midway2/cnajar/ --outdir PeakCallingPerInd/{Phenotype}/ --name {IndID}.{Rep} --broad"
    log:
        "logs/Macs2PeakCalling_broad_perind/{Phenotype}/{IndID}.{Rep}.log"

