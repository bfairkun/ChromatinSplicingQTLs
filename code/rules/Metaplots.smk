# rule GetInterestingRegionsForMetaplot:
#     input:
#         # ColocOutput
#     output:
#         beds = expand("Metaplots/Regions/InterestingRegions/{Region}.bed", Region="chromatin_eQTLs", "splice_not_chromatinQTLs")
#     log:
#         "logs/GetInterestingRegionsForMetaplot.log"
#     shell:
#         """
#         Rscript
#         """

rule CreateBigwigLists:
    input:
        colors = "../data/ColorsForPhenotypes.xlsx",
        bw = rules.GatherAllBigwigs.output,
    output:
        bw_list = "Metaplots/bwList.tsv",
        bw_list_allunstranded = "Metaplots/bwList.allunstranded.tsv",
        bw_groups = "Metaplots/bwGroups.tsv",
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/MakeBigwigList.R
        """

rule GtfSubset:
    input:
        gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf",
        genes = "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        bed12 = "ExpressionAnalysis/polyA/ExpressedGeneList.bed12.bed.gz",
        tbi = "ExpressionAnalysis/polyA/ExpressedGeneList.bed12.bed.gz.tbi",
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        grep -f <(awk -F'\\t' '{{print $4}}' {input.genes} ) {input.gtf} | bedparse gtf2bed --extraFields gene_id | sort -u -k13,13 | awk -F'\\t' -v OFS='\\t' '{{$4=$13; print $0}}' | cut -d$'\\t' -f1-12 | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed12}
        tabix -p bed {output.bed12}
        """


rule CreateMetaQTL_Bigwigs:
    input:
        bw = rules.GatherAllBigwigs.output,
        bw_list = "Metaplots/bwList.allunstranded.tsv",
        bw_groups = "../data/Metaplots/GroupsFiles/{GroupsFile}.tsv",
        vcf = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz.tbi",
        bed = "../data/Metaplots/Regions/{Region}.bed"
    output:
        ## I should define the output
        directory("Metaplots/IntermediateBigwigs/{Region}.{GroupsFile}")
    log:
        "logs/CreateMetaQTL_Bigwigs/{Region}.{GroupsFile}.log"
    params:
        flanking_bp = 10000
    conda:
        "../scripts/GenometracksByGenotype/GenometracksByGenotype.conda_env.yml"
    resources:
        # mem = much_more_mem_after_first_attempt
        mem = 48000
    shell:
        """
        mkdir -p {output}
        python scripts/GenometracksByGenotype/AggregateForQTLMetaplot.py --QTLsBed {input.bed} --FlankingRegionLengthBp {params.flanking_bp} --VCF {input.vcf} --OutputPrefix Metaplots/IntermediateBigwigs/{wildcards.Region}.{wildcards.GroupsFile}/ --BigwigList {input.bw_list} --GroupSettingsFile {input.bw_groups} &> {log}
        """

rule GatherMetaQTL_Bigwigs:
    input:
        expand("Metaplots/IntermediateBigwigs/{Region}.{GroupsFile}", Region=["sQTL.eQTL.Annotated_basic", "sQTL.eQTL.Annotated_NMD", "sQTL.eQTL.Annotated_Not_basic", "sQTL.eQTL.Unannotated"], GroupsFile="H3K27AC_AndAllRNA")
