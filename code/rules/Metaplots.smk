rule GetInterestingRegionsForMetaplot:
    input:
        # ColocOutput
    output:
        beds = expand("Metaplots/Regions/InterestingRegions/{Region}.bed", Region="chromatin_eQTLs", "splice_not_chromatinQTLs")
    log:
        "logs/GetInterestingRegionsForMetaplot.log"
    shell:
        """
        Rscript
        """

rule CreateBigwigLists:
    input:
        bw = rules.GatherAllBigwigs.output,
        colors = "../data/ColorsForPhenotypes.xlsx"
    output:
        bw_list = "Metaplots/bwList.tsv",
        bw_list_allunstranded = "Metaplots/bwList.allunstranded.tsv"
        bw_groups = "Metaplots/bwGroups.tsv",
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript
        """


rule CreateMetaQTL_Bigwigs:
    input:
        bw = rules.GatherAllBigwigs.output,
        bw_list = "Metaplots/bwList.tsv",
        bw_groups = "Metaplots/bwGroups.tsv"
        vcf = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz",
        tbi = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz.tbi",
        bed = "Metaplots/Regions/InterestingRegions/{Region}.bed"
    output:
        ## I should define the output
        # "Metaplots/Regions/InterestingRegions_bigwigs/{Region}/"
    log:
        "logs/CreateMetaQTL_Bigwigs/{Region}.log"
    params:
        flanking_bp = 10000
    conda:
        "../scripts/GenometracksByGenotype/GenometracksByGenotype.conda_env.yml"
    shell:
        """
        python scripts/GenometracksByGenotype/AggregateForQTLMetaplot.py --QTLsBed {input.bed} --FlankingRegionLengthBp {params.flanking_bp} --VCF {input.vcf} --OutputPrefix Metaplots/Regions/InterestingRegions_bigwigs/{wildcards.Region}/ --BigwigList {input.bw_list} --GroupSettingsFile {input.bw_groups} &> {log}
        """
