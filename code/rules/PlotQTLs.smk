rule MakeBigwigListTsv:
    """
    track for all test features, samples tsv file, and track for all coloc features colored by cluster/gene, and filter junc file with colocalized introns
    """
    input:
        # beds = expand("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz", Phenotype = PhenotypesToColoc),
        samplelist = "config/samples.tsv",
        # bigwigs = GatherAllBigwigs,
        # hyprccoloc_results = "../output/hyprcoloc_results/ForColoc/hyprcoloc.results.txt.gz",
        hyprccoloc_results = "scratch/hyprcoloc.results.txt.gz",
        leafcutter_numers = "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    output:
        "PlotQTLs/BigwigList.tsv"
    # conda:
    #     "../envs/r_slopes.yml"
    params:
        PhenotypesToColoc = " ".join(PhenotypesToColoc)
    log:
        "logs/MakeBigwigListTsv.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/MakeQTLPlotTables.R PlotQTLs/ {input.hyprccoloc_results} {input.leafcutter_numers} '{params.PhenotypesToColoc}' &> {log}
        """

# PrefixOut <- args[1]
# hyprcoloc_results_fn <- args[2]
# leafcutter_count_numers <- args[3]
# Phenotypes <- unlist(strsplit(args[4], ' '))

# rule MakeColocTest
# 
