
# use tabix genewise qtltools output to ascertain p value in trait2 for top SNP determined from trait1 discoveries.
# There are so many potential trait pairs, so i am going to break this down into chunks to parralelize
NumPvalsForPi1Chunks = 10
rule CalculatePi1_GetTraitPairs_AllTraits:
    input:
        Nominal = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.tabix.gz", Pass="NominalPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=MyPhenotypes),
        Permutation = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz", Pass="PermutationPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=MyPhenotypes)
    params:
        NumChunks = NumPvalsForPi1Chunks
    output:
        expand("pi1/PairwiseTraitsToCompare/{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/CalculatePi1_GetTraitPairs_AllTraits.R {params.NumChunks} pi/PairwiseTraitsToCompare/
        """

rule GatherChunks:
    input:
        expand("pi1/PairwiseTraitsToCompare/{n}.txt.gz", n=range(1, 1+NumPvalsForPi1Chunks))

rule GetPvalsForPi1AllTraitPairs:
    input:
        "scratch/PairwisePi1Traits.{chunk}.txt.gz"
    output:
        "scratch/PairwisePi1Traits.P.{chunk}.txt.gz"
    log:
        "logs/GetPvalsForPi1AllTraitPairs/{chunk}.log"
    shell:
        """
        python scripts/CalculatePi1_GetAscertainmentP_AllPairs.py {input} {output} &> {log}
        """

rule GatherPvalsForPi1AllTraitPairs:
    input:
        expand("scratch/PairwisePi1Traits.P.{chunk}.txt.gz", chunk=range(1, NumPvalsForPi1Chunks+1))

rule GetPvalsForPiStats:
    """
    
    """
    input:
        Nominal = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.tabix.gz", Pass="NominalPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=MyPhenotypes),
        Permutation = expand("QTLs/QTLTools/{Phenotype}/{Pass}{QTLsGenotypeSet}{FeatureCoordinatesRedefinedFor}.txt.gz", Pass="PermutationPass", QTLsGenotypeSet="", FeatureCoordinatesRedefinedFor="ForColoc", Phenotype=MyPhenotypes)
    output:
        "QC/PvalsForPi1_{FDR}.txt.gz"
    params:
        phenotypes = ' '.join(MyPhenotypes)
    log:
        "logs/GetPvalsForPiStats/{FDR}.log"
    shell:
        """
        python scripts/GatherTopSNPPvalsCrossTraits.py {output} {wildcards.FDR} {params.phenotypes} &> {log}
        """
