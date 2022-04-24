rule MakeBigwigListTsv:
    """
    track for all test features, samples tsv file, and track for all coloc features colored by cluster/gene, and filter junc file with colocalized introns
    """
    input:
        beds = expand("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz", Phenotype = PhenotypesToColoc),
        samplelist = "config/samples.tsv",
        colors = "../data/ColorsForPhenotypes.xlsx",
        # bigwigs = GatherAllBigwigs,
        hyprccoloc_results = "../output/hyprcoloc_results/ForColoc/MolColocStandard/hyprcoloc.results.txt.gz",
    output:
        "PlotQTLs/ColocFeatures.bed",
        "PlotQTLs/ColocTestFeatures.bed",
        "PlotQTLs/bwList.tsv",
        "PlotQTLs/bwList.Groups.tsv",
    conda:
        "../envs/r_slopes.yml"
    params:
        PhenotypesToColoc = " ".join(PhenotypesToColoc)
    log:
        "logs/MakeBigwigListTsv.log"
    shell:
        """
        Rscript scripts/MakeQTLPlotTables.R PlotQTLs/ {input.hyprccoloc_results}  '{params.PhenotypesToColoc}' &> {log}
        """

rule indexFullBeds:
    input:
        coloc_feats = "PlotQTLs/ColocFeatures.bed",
        test_feats = "PlotQTLs/ColocTestFeatures.bed",
    output:
        coloc_feats = "PlotQTLs/ColocFeatures.bed.gz",
        test_feats = "PlotQTLs/ColocTestFeatures.bed.gz",
        coloc_feats_tbi = "PlotQTLs/ColocFeatures.bed.gz.tbi",
        test_feats_tbi = "PlotQTLs/ColocTestFeatures.bed.gz.tbi",
    shell:
        """
        bgzip {input.coloc_feats}
        bgzip {input.test_feats}
        tabix -p bed {output.coloc_feats}
        tabix -p bed {output.test_feats}
        """


rule CollectNormalizedPsiTables:
    input:
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz.tbi", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"])

rule PlotColocalizedQTLs:
    input:
        hyprccoloc_results = "../output/hyprcoloc_results/ForColoc/MolColocStandard/hyprcoloc.results.txt.gz",
        bwList = "PlotQTLs/bwList.tsv",
        groups = "PlotQTLs/bwList.Groups.tsv",
        coloc_feats = "PlotQTLs/ColocFeatures.bed.gz.tbi",
        test_feats = "PlotQTLs/ColocTestFeatures.bed.gz.tbi",
    log:
        "logs/PlotColocalizedQTLs/{n}.log"
    output:
        touch("PlotQTLs/ColocQTLs/touchfile.{n}.txt")
    conda:
        "../scripts/GenometracksByGenotype/GenometracksByGenotype.conda_env.yml"
    shell:
        """
        mkdir -p PlotQTLs/Temp/
        python scripts/PlotQTLs.py PlotQTLs/Temp/ PlotQTLs/ColocQTLs/ {wildcards.n} 50 &> {log}
        """

rule GatherPlotColocalizedQTLs:
    input:
        expand("PlotQTLs/ColocQTLs/touchfile.{n}.txt", n=[i+1 for i in range(50)])

