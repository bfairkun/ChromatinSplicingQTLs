
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
