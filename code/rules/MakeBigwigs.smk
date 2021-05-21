# rule CalculateNormFactorsForBigwig:
#     input:
#         "featureCounts/{Phenotype}/Counts.txt"
#     output:
#         "featureCounts/{Phenotype}.tsv"
#     log:
#         "logs/CalculateNormFactorsForBigwig.log"
#     conda:
#         "../envs/r_essentials.yml"
#     shell:
#         """
#         shell
#         """

rule MakeBigwigs:
    """
    Scale bigwig to base coverage per billion chromosomal reads
    """
    input:
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        bam = GetBamForBigwig,
        bai = GetBaiForBigwig
    params:
        GetBigwigParams
    output:
        bw = "bigwigs/{Phenotype}/{IndID}.{Rep}.bw"
    log:
        "logs/MakeBigwigs/{Phenotype}/{IndID}.{Rep}.log"
    resources:
        mem = 48000
    shell:
        """
        scripts/GenometracksByGenotype/BamToBigwig.sh {input.fai} {input.bam} {output.bw} {params} -scale $(bc <<< "scale=3;1000000000/$(samtools idxstats {input.bam} | awk '$1 ~ "^chr" {{sum+=$2}} END{{printf sum}}')") &> {log}
        """
