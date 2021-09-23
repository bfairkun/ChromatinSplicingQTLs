rule Count_sQTLs_with_SpliceSiteSNPs:
    """
    An example rule using a jupyter notebook.
    See the official documentation for more:
    https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration
    """
    input:
        # As per usual, the notebook can rely on input files
        "QTLs/QTLTools/polyA.Splicing/PermutationPass.FDR_Added.txt.gz",
        "QTLs/QTLTools/polyA.Splicing/OnlyFirstReps.sorted.qqnorm.bed.gz",
        "QTLs/QTLTools/Expression.Splicing/Genotypes/WholeGenome.vcf.gz"
    output:
        # And the notebook can write out output files
        "QTLs/QTLTools/polyA.Splicing/PermutationPass.FDR_Added.SS_SNPs.Annotated.txt.gz",
    log:
        # Save a copy of the ipynb to ../docs/ to render on the website for easy sharing
        notebook = "../docs/20210921_CountSpliceSiteSNPsInSQTLs.py.ipynb"
    conda:
        # Use a conda env defined in yml to run notebook
        # Unlike input, output, and log directives, the path to yml needs to be
        # relative to the .smk file containing this rule (eg in
        # /code/rules/*.smk)
        "../envs/jupyter.yml"
    notebook:
        # Path to notebook needs to be relative to the .smk file (eg in /code/rules/*.smk)
        # Unlike input, output, and log directives, the path to notebook  needs to be
        # relative to the .smk file containing this rule (eg in
        # /code/rules/*.smk).
        # (But reference relative paths in notebook relative to the code
        # directory (eg the working directory that snakemake gets executed in)
        "../notebooks/20210921_CountSpliceSiteSNPsInSQTLs.py.ipynb"

