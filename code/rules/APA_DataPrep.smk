rule MakeHg19APA_Tables:
    input:
        "../data/APA_dat/APApeak_BrianasPermutationTestResults.Nuclear.tsv.gz",
        "../data/APA_dat/APApeak_BrianasPermutationTestResults.Total.tsv.gz",
        "../data/APA_dat/APApeak_Phenotype_GeneLocAnno.Nuclear.tsv.gz",
        "../data/APA_dat/APApeak_Phenotype_GeneLocAnno.Total.tsv.gz",
    output:
        "APA_Processing/PhenotypeTables/Nuclear.hg19.bed",
        "APA_Processing/PhenotypeTables/Total.hg19.bed"
    log:
        "logs/MakeHg19APA_Tables.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/Prep_hg19_APA_dat.R &> {log}
        """

rule liftoverAPA_Tables:
    input:
        bed = "APA_Processing/PhenotypeTables/{Fraction}.hg19.bed",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    wildcard_constraints:
        Fraction = "Nuclear|Total"
    conda:
        "../envs/liftover.yml"
    output:
        "QTLs/QTLTools/APA_{Fraction}/OnlyFirstReps.qqnorm.bed.gz"
    shell:
        """
        cat <(cat {input.bed} | head -1) <(liftOver -bedPlus=6 {input.bed} {input.chain} /dev/stdout /dev/null) | gzip - > {output}
        """

rule GatherLiftedOverAPA_Tables:
    input:
        expand("QTLs/QTLTools/APA_{Fraction}/OnlyFirstReps.sorted.qqnorm.bed.gz", Fraction=["Nuclear", "Total"])

# rule GetSeqContainingPutativePASMotifs:
#     input:
#         bed = expand("QTLs/QTLTools/APA_{Fraction}/OnlyFirstReps.sorted.qqnorm.bed.gz", Fraction=["Nuclear", "Total"]),
#         fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
#         fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
#     output:
#     conda:
#         ""
#     shell:
#         """
#         zcat {input.bed} | awk -v OFS='\\t' '$1 !~ "^#" {{print $1,$2,$3,$4,$5,$6}}' | sort | uniq | grep -P '^chr[0-9]+\\t' |  bedtools sort -i - | bedtools flank -s -l 0  -r 50 -g {input.fai} | bedtools shift -i - -m 50 -p -50 -g {input.fai} | shuf -n 1000 | bedtools sort -i - | bedtools getfasta -s -fi {input.fa} -bed - | meme /dev/stdin -oc scratch/memetest -w 7 -nmotifs 1
#         """



