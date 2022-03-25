rule MakeBigwigListTsv:
    """
    track for all test features, samples tsv file, and track for all coloc features colored by cluster/gene, and filter junc file with colocalized introns
    """
    input:
        "config/samples.tsv"
    output:
        "PlotQTLs/BigwigList.tsv"
    log:
        "logs/MakeBigwigListTsv.log"
    shell:
        """
        shell
        """

# rule MakeColocTest
# 
