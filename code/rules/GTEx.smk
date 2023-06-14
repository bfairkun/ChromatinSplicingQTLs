rule DownloadGTExCounts:
    output:
        "GTEx/data/gene_reads_2017-06-05_v8_{tissue}.gct.gz"
    log:
        "logs/GTEx/download_gene_counts.{tissue}.log"
    wildcard_constraints:
        tissue = '|'.join(gtex_tissues)
    shell:
        """
        (wget -O {output} https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_{wildcards.tissue}.gct.gz) &> {log}
        """
    
rule DownloadFromGTEx_VCF:
    input:
        manifest = "GTEx/manifest/file-manifest.json",
    output:
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz", 
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz.tbi",
    log:
        'logs/GTEx/download_vcf.log' 
    resources:
        mem_mb = 4000
    shell:
        """
        ({config[client]} download-multiple --profile=AnVIL --manifest={input.manifest} --download-path=/project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/GTEx/data/ --protocol=s3) &> {log}
        """

rule PrepareGTExPhenotypesForQTLTools:
    input:
        "GTEx/data/gene_reads_2017-06-05_v8_{tissue}.gct.gz",
        "ExpressionAnalysis/polyA/ExpressedGeneList.txt"
    output:
        "GTEx/QTLs/{tissue}/cpm.bed.gz",
        "GTEx/QTLs/{tissue}/qqnorm.bed.gz"
    log:
        "logs/GTEx/QTLs/prepare_input.{tissue}.log"
    resources:
        mem_mb = 8000
    shell:
        """
        Rscript scripts/PrepareGTExPhenotypes.R {input} {output} &> {log}
        """
        
use rule SortQTLtoolsPhenotypeTable as SortQTLtoolsPhenotypeTable_GTEx with:
    input:
        "GTEx/QTLs/{tissue}/{norm}.bed.gz"
    output:
        bed = "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz",
        tbi = "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz.tbi"
    log:
        "logs/GTEx/QTLs/sort.{tissue}.{norm}.log"
    wildcard_constraints:
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'
        
use rule PhenotypePCs as PhenotypePCs_GTEx with:
    input:
        "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz",
    output:
        "GTEx/QTLs/{tissue}/{norm}.sorted.bed.pca",
    log:
        "logs/GTEx/QTLs/PCs.{tissue}.{norm}.log"
    wildcard_constraints:
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'

use rule PlotPhenotypePCs as PlotPhenotypePCs_GTEx with:
    input:
        "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz",
    output:
        "GTEx/QTLs/{tissue}/{norm}.sorted.bed.pca.pdf"
    log:
        "logs/GTEx/QTLs/PCs_plot.{tissue}.{norm}.log"
    wildcard_constraints:
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'
        


use rule QTLtools_generalized as QTLtools_GTEx with:
    input:
        vcf = "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz",
        tbi = "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        bed = "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz",
        bed_tbi = "GTEx/QTLs/{tissue}/{norm}.sorted.bed.gz.tbi",
        cov = "GTEx/QTLs/{tissue}/{norm}.sorted.bed.pca"
    output:
        temp("GTEx/QTLs/{tissue}/{Pass}Chunks/{QTLTools_chunk_n}{norm}.txt")
    log:
        "logs/GTEx/QTLs/{tissue}/{Pass}Chunks/{QTLTools_chunk_n}{norm}.log"
    params:
        WindowFlag = "--window 100000",
        OtherFlags = "",
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = "",
    wildcard_constraints:
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "PermutationPass|NominalPass",
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'
 
use rule Gather_QTLtools_cis_pass as Gather_QTLtools_GTEx with:
    input:
        expand("GTEx/QTLs/{{tissue}}/{{Pass}}Chunks/{QTLTools_chunk_n}{{norm}}.txt", QTLTools_chunk_n=ChunkNumbers)
    output:
        "GTEx/QTLs/{tissue}/{Pass}.{norm}.txt.gz"
    wildcard_constraints:
        Pass = "PermutationPass|NominalPass",
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'
    log:
        "logs/GTEx/QTLs/{tissue}/{Pass}.{norm}.collect.log"


use rule AddQValueToPermutationPass as AddQValueToPermutationPass_GTEx with:
    input:
        "GTEx/QTLs/{tissue}/{Pass}.{norm}.txt.gz"
    output:
        table = "GTEx/QTLs/{tissue}/{Pass}.{norm}.FDR_Added.txt.gz"
    log:
        "logs/GTEx/QTLs/AddQValueToPermutationPass/{tissue}.{norm}.collect.{Pass}.log"
    wildcard_constraints:
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm',
        Pass = 'PermutationPass'



use rule tabixNominalPassQTLResults as tabixNominalPassQTLResults_GTEx with:
    input:
        "GTEx/QTLs/{tissue}/NominalPass.{norm}.txt.gz"
    wildcard_constraints:
        tissues = '|'.join(gtex_tissues),
        norm = 'cpm|qqnorm'
    output:
        txt = "GTEx/QTLs/{tissue}/NominalPass.{norm}.txt.tabix.gz",
        tbi = "GTEx/QTLs/{tissue}/NominalPass.{norm}.txt.tabix.gz.tbi"
    log:
        "logs/GTEx/QTLs/tabixNominalPassQTLResults/{tissue}.NominalPass.{norm}.log"



rule collect_gtex:
    input:
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz", 
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz.tbi",
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz",
        "GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz.tbi",
        expand("GTEx/data/gene_reads_2017-06-05_v8_{tissue}.gct.gz", tissue = gtex_tissues),
        expand("GTEx/QTLs/{tissue}/cpm.bed.gz", tissue = gtex_tissues),
        expand("GTEx/QTLs/{tissue}/NominalPass.{norm}.txt.tabix.gz", tissue=gtex_tissues, norm = ['cpm', 'qqnorm']),
        expand("GTEx/QTLs/{tissue}/PermutationPass.{norm}.FDR_Added.txt.gz", tissue=gtex_tissues, norm = ['cpm', 'qqnorm']),
        
        