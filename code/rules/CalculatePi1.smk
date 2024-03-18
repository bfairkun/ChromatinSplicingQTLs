
# use tabix genewise qtltools output to ascertain p value in trait2 for top SNP determined from trait1 discoveries.
# There are so many potential trait pairs, so i am going to break this down into chunks to parralelize
rule CalculatePi1_GetTraitPairs_AllTraits:
    input:
        Nominal = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.tabix.gz", Pass="NominalPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=PhenotypesToColoc),
        Permutation = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz", Pass="PermutationPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=PhenotypesToColoc)
    params:
        NumChunks = NumPvalsForPi1Chunks
    output:
        expand("pi1/PairwiseTraitsToCompare/{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))
    conda:
        "../envs/r_essentials.yml"
    resources:
        mem_mb = 32000
    log:
        "logs/CalculatePi1_GetTraitPairs_AllTraits.log"
    shell:
        """
        Rscript scripts/CalculatePi1_GetTraitPairs_AllTraits.R {params.NumChunks} pi1/PairwiseTraitsToCompare/ {input.Permutation} &> {log}
        """

rule GatherChunks:
    input:
        expand("pi1/PairwiseTraitsToCompare/{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))

rule GetPvalsForPi1AllTraitPairs:
    input:
        TraitsToCompare = "pi1/PairwiseTraitsToCompare/{chunk}.txt.gz",
        tabix_QTLsOut = expand("QTLs/QTLTools/{Phenotype}/NominalPassForColoc.txt.tabix.gz", Phenotype=PhenotypesToColoc)
    output:
        "pi1/PairwiseTraitsToCompare/P.{chunk}.txt.gz"
    log:
        "logs/GetPvalsForPi1AllTraitPairs/{chunk}.log"
    shell:
        """
        python scripts/CalculatePi1_GetAscertainmentP_AllPairs.py {input.TraitsToCompare} {output} {input.tabix_QTLsOut} &> {log}
        """

# use rule GetPvalsForPi1AllTraitPairs as GetPvalsForPi1AllTraitPairs_Unstandardized with:
#     input:
#         TraitsToCompare = "pi1/PairwiseTraitsToCompare/{chunk}.txt.gz",
#         tabix_QTLsOut = ["QTLs/QTLTools/polyA.Splicing/NominalPassForColocUnstandardized.txt.tabix.gz", "QTLs/QTLTools/polyA.Splicing/NominalPassForColocUnstandardized.txt.tabix.gz"


rule GatherPvalsForPi1AllTraitPairs:
    input:
        # expand("scratch/PairwisePi1Traits.P.{chunk}.txt.gz", chunk=range(1, NumPvalsForPi1Chunks+1))
        expand("pi1/PairwiseTraitsToCompare/P.{chunk}.txt.gz", chunk=range(1, NumPvalsForPi1Chunks+1))
    output:
        "pi1/PairwisePi1Traits.P.all.txt.gz"
    shell:
        """
        cat <(zcat {input[0]} | head -1) <(zcat {input} | grep -v -P '^PC1\\t') | gzip - > {output}
        """

rule PlotPi1Heatmaps:
    input:
        QTLs = "pi1/PairwisePi1Traits.P.all.txt.gz",
        Peaks = expand("Misc/PeaksClosestToTSS/{ChIP_Phenotype}_assigned.tsv.gz", ChIP_Phenotype = ["H3K27AC", "H3K4ME3", "H3K4ME1"])
    output:
        dat = "pi1/DatForHeatmapPlot.tsv.gz",
        P = "pi1/DatForHeatmapPlot.pdf"
    log:
        "logs/PlotPi1Heatmaps.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/CalculatePi1_PlotHeatmap.R {output.dat} {output.P} &> {log}
        """

rule Plot_eQTL_QQ:
    input:
        QTLs = "pi1/PairwisePi1Traits.P.all.txt.gz",
        TestSNPs = "QTLs/QTLTools/Expression.Splicing/NominalPassForColoc.RandomSamplePvals.txt.gz",
        sQTLs = "SplicingAnalysis/sQTLs_p_and_u.tsv.gz"
    output:
        dat = "Misc/eQTL_qq/dat.tsv.gz",
        P = "Misc/eQTL_qq/Plot.pdf"
    log:
        "logs/Plot_eQTL_QQ.log"
    conda:
        "../envs/r_scattermore.yml"
    shell:
        """
        Rscript scripts/Plot_eQTL_QQ.R {output.dat} {output.P} &> {log}
        """


# rule Plot

# rule CalculatePi1:
#     input:
