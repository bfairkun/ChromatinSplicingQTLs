
rule DownloadGruberQTLs:
    output:
        "PlotGruberQTLs/Data/localQTLs_{pheno}.FDR0.1.bedpe"
    shell:
        """
        wget -O-  http://mitra.stanford.edu/kundaje/portal/chromovar3d/QTLs/localQTL/localQTLs_{wildcards.pheno}.FDR0.1.bedpe.gz | zcat | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,$4,$5,$6,$7,".", ".", "."}}' > {output}
        """

rule LiftOverGruberQTLs:
    input:
        bedpe = "PlotGruberQTLs/Data/localQTLs_{pheno}.FDR0.1.bedpe",
        chain = "../data/hg19ToHg38.over.chain.gz"
    output:
        "PlotGruberQTLs/Data/localQTLs_{pheno}.FDR0.1.hg38.bedpe"
    log:
        "logs/PlotGruberExamples/LiftOverGruberQTLs.{pheno}.log"
    shell:
        """
        python2.7 scripts/liftOverBedpe/liftOverBedpe.py --i {input.bedpe} --o {output} --v T --h F --lift liftOver --chain {input.chain} &> {log}
        """

rule GetGruberPeaks:
    input:
        "PlotGruberQTLs/Data/localQTLs_{pheno}.FDR0.1.hg38.bedpe"
    output:
        "PlotGruberQTLs/Data/{pheno}.peaks.bed"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '{{print $4,$5,$6}}' {input} | sort | uniq | bedtools sort -i - > {output}
        """

rule gtf_to_bed:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    output:
        "PlotGruberQTLs/Data/genes.gtf.bed12"
    shell:
        """
        python scripts/fromgtfTobed12.py {input} --useGene --mergeTranscripts --output {output}
        """

# rule MakeGeneAndPeakTracks:
#     input:
#         gene = "PlotGruberQTLs/Data/genes.gtf.bed12",
#         peaks = "PlotGruberQTLs/Data/H3K4ME3.peaks.bed"
#     conda:
#         "../envs/GenometracksByGenotype.yml"
#     output:
#         "PlotGruberQTLs/Data/Test/H3K4me3.tracks.ini"
#     shell:
#         "make_tracks_file -f {input.peaks} -o {output}"


rule MakePeaksTracks:
    input:
        "PlotGruberQTLs/Data/{pheno}.peaks.bed"
    output:
        "PlotGruberQTLs/Data/Tracks/{pheno}.tracks.ini"
    run:
        template = """
[{pheno} peaks]
file = {peaks_bed}
title = {pheno} peaks
height = 1
color = darkblue
labels = false
fontsize = 10
file_type = bed
gene_rows = 3
        """ 
        with open(output[0],'w') as fh:
            fh.write(template.format(pheno=wildcards.pheno, peaks_bed=input[0]))

rule MakePlotShellscript:
    input:
        vcf = "Genotypes/GEUVADIS_Lappalainnen.vcf.gz",
        bw = expand( "Bigwigs/GEUVADIS_RNAseq/{RNASeqSample}.bw", RNASeqSample=GEUVADIS_line_fastq_dict.keys()),
        gene = "../data/pygenometracks/gtf_tracks.ini",
        peaks ="PlotGruberQTLs/Data/Tracks/{pheno}.tracks.ini" ,
        H3K4me3QTL = "PlotGruberQTLs/Data/localQTLs_{pheno}.FDR0.1.hg38.bedpe",
        eQTL = "../data/QTLBase.GEUVADIS.eQTLs.hg38.txt.gz"
    output:
        TopQTLsScript ="PlotGruberQTLs/Plots/Plot{pheno}TopQTL.sh",
        SharedQTLsScript = "PlotGruberQTLs/Plots/Plot{pheno}QTL_eQTLs.sh"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/VisualizeGrubertHistoneQTLs.R {wildcards.pheno} 
        """

rule PlotFromShellScript:
    input:
        TopQTLsScript ="PlotGruberQTLs/Plots/Plot{pheno}TopQTL.sh",
        SharedQTLsScript = "PlotGruberQTLs/Plots/Plot{pheno}QTL_eQTLs.sh"
    output:
        "logs/PlotGruberQTLs/PlotFromShellScript.{pheno}.log"
    conda:
        "../envs/GenometracksByGenotype.yml"
    shell:
        """
        bash {input.TopQTLsScript} &> {output}
        bash {input.SharedQTLsScript} >> {output} 2>&1
        """
 



