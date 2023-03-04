rule MakeGeneEnsembl2Symbol:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf"
    output:
        "NonCodingRNA/annotation/genes.txt.gz"
    log:
        "logs/NonCodingRNA/genes.txt.log"
    shell:
        """
        (awk '$3=="gene"' {input} | awk -F'"' '{{print $2, $6}}' OFS='\\t' - | gzip - > {output}) &> {log}
        """
        
rule SummarizeUaRNA:
    input:
        'NonCodingRNA/annotation/NonCodingRNA.annotation.tab.gz',
        'NonCodingRNA/annotation/Gencode.uaRNA.annotation.tab.gz',
        'NonCodingRNA/annotation/ncRNA.bed.gz',
        expand('QTLs/QTLTools/{Phenotype}/{Pass}.txt.gz',
               Phenotype = ["chRNA.Expression_ncRNA", "chRNA.Expression.Splicing"],
               Pass = ["NominalPass", "PermutationPass.FDR_Added"]
            ),
        'NonCodingRNA/annotation/allGenes.TSS_bp.bed.gz',
        'NonCodingRNA/annotation/NonCodingRNA.bed.gz',
        'NonCodingRNA/annotation/tmp/allGenes.Gencode.bed.gz',
        'RPKM_tables/chRNA.RPKM.bed.gz',
    output:
        "NonCodingRNA/QTLs/summary.uaRNA.tab.gz"
    log:
        "logs/NonCodingRNA/Summarize.uaRNA.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 36000
    shell:
        """
        python scripts/NonCodingRNA/summarize_uaRNA_table.py &> {log}
        """
        
rule uaRNAMetaplotInput:
    input:
        'NonCodingRNA/QTLs/summary.uaRNA.tab.gz',
        'scripts/GenometracksByGenotype/PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz',
        'NonCodingRNA/annotation/genes.txt.gz'
    output:
        'NonCodingRNA/QTLs/apQTLs.ncQTL.uaGene.bed',
        'NonCodingRNA/QTLs/diQTLs.ncQTL.uaGene.bed',
        'NonCodingRNA/QTLs/apQTLs.eQTL.uaGene.bed',
        'NonCodingRNA/QTLs/diQTLs.eQTL.uaGene.bed',
        'NonCodingRNA/QTLs/apQTLs.ncQTL.uaRNA.bed',
        'NonCodingRNA/QTLs/diQTLs.ncQTL.uaRNA.bed',
        'NonCodingRNA/QTLs/apQTLs.eQTL.uaRNA.bed',
        'NonCodingRNA/QTLs/diQTLs.eQTL.uaRNA.bed',
        'NonCodingRNA/QTLs/apQTLs.ncQTL.uaGene.bed12',
        'NonCodingRNA/QTLs/diQTLs.ncQTL.uaGene.bed12',
        'NonCodingRNA/QTLs/apQTLs.eQTL.uaGene.bed12',
        'NonCodingRNA/QTLs/diQTLs.eQTL.uaGene.bed12'
    log:
        "logs/NonCodingRNA/uaRNA.metaplot.input.log"
    conda:
        "../envs/py_tools.yml"
    resources:
        mem_mb = 36000
    shell:
        """
        python scripts/NonCodingRNA/get_uaRNA_metaplot_input.py &> {log}
        """
    
rule RunAggregateForMetaplot:
    input:
        QTLsBed = 'NonCodingRNA/QTLs/{QTLs}.{var}.{phe}.bed',
        vcf = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/Genotypes/1KG_GRCh38/Autosomes.vcf.gz",
        bwList = '../data/ncQTL.bwList.tsv',
        Groups = '../data/ncQTL.Groups.tsv'
    output:
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA-+-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA-+-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA-+-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA---High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA---Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA---Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-polyA.RNA-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-polyA.RNA-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-polyA.RNA-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K27AC-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K27AC-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K27AC-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K4ME3-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K4ME3-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K4ME3-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K36ME3-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K36ME3-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-H3K36ME3-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-ProCap-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-ProCap-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-ProCap-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-CTCF-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-CTCF-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-CTCF-.-Low.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-DNaseISensitivity-.-High.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-DNaseISensitivity-.-Mid.bw',
        'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-DNaseISensitivity-.-Low.bw',
    wildcard_constraints:
        QTLs = 'apQTLs|diQTLs',
        var = 'ncQTL|eQTL',
        phe = 'uaGene|uaRNA'
    log:
        "logs/NonCodingRNA/metaplots/Aggregate.{QTLs}.{var}.{phe}.log"
    conda:
        "../envs/GenometracksByGenotype.yml"
    resources:
        mem_mb = 36000
    shell:
        """
        python scripts/GenometracksByGenotype/AggregateForQTLMetaplot.py --QTLsBed {input.QTLsBed} --FlankingRegionLengthBp 10000 --Workdir /project2/yangili1/cfbuenabadn/ChromatinSplicingQTLs/code/ --VCF {input.vcf} --OutputPrefix NonCodingRNA/QTLs/metaplots/{wildcards.QTLs}.{wildcards.var}.{wildcards.phe}/ --BigwigList {input.bwList} --GroupSettingsFile {input.Groups} -vv &> {log}
        """


def GetMetaplotBed(wildcards):
    if wildcards.phe == 'uaRNA':
        return 'NonCodingRNA/QTLs/{QTLs}.{var}.uaRNA.bed'.format(QTLs = wildcards.QTLs, var = wildcards.var)
    else:
        return 'NonCodingRNA/QTLs/{QTLs}.{var}.uaGene.bed12'.format(QTLs = wildcards.QTLs, var = wildcards.var)

def GetBwListForCommand(wildcards):
    bw_list = []
    for strand in ['+', '-']:
        for Size in ['High', 'Mid', 'Low']:
            file_bw = 'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-chRNA-{strand}-{Size}.bw'
            file_bw = file_bw.format(QTLs = wildcards.QTLs, var = wildcards.var, 
                                     phe=wildcards.phe, strand = strand, Size = Size)
            bw_list.append(file_bw)
            
    for Phenotype in ['polyA.RNA', 'H3K27AC', 'H3K4ME3', 'H3K36ME3', 'ProCap', 'CTCF', 'DNaseISensitivity']:
        for Size in ['High', 'Mid', 'Low']:
            file_bw = 'NonCodingRNA/QTLs/metaplots/{QTLs}.{var}.{phe}/Summarised-{Phenotype}-.-{Size}.bw'
            file_bw = file_bw.format(QTLs = wildcards.QTLs, var = wildcards.var, 
                                     phe=wildcards.phe, Phenotype = Phenotype, Size = Size)
            bw_list.append(file_bw)
            
    out_bw = ' '.join(bw_list)
    return out_bw
        

rule ComputeMatrixForMetaplot:
    input:
        expand('NonCodingRNA/QTLs/metaplots/{{QTLs}}.{{var}}.{{phe}}/Summarised-chRNA-{strand}-{Size}.bw',
                strand = ['+', '-'], Size = ['High', 'Mid', 'Low']),
        rbed = GetMetaplotBed,
        bw_list = expand(
            'NonCodingRNA/QTLs/metaplots/{{QTLs}}.{{var}}.{{phe}}/Summarised-{Phenotype}-.-{Size}.bw', 
            Phenotype = ['polyA.RNA', 'H3K27AC', 'H3K4ME3', 'H3K36ME3', 'ProCap', 'CTCF', 'DNaseISensitivity'], 
            Size = ['High', 'Mid', 'Low']
        ),
    output:
        "NonCodingRNA/QTLs/metaplots/mat/{QTLs}.{var}.{phe}.mat.gz"
    wildcard_constraints:
        QTLs = 'apQTLs|diQTLs',
        var = 'ncQTL|eQTL',
        phe = 'uaGene|uaRNA'
    params:
        bwList_input = GetBwListForCommand
    log:
        "logs/NonCodingRNA/metaplots/ComputeMatrix.{QTLs}.{var}.{phe}.log"
    resources:
        mem_mb = 36000
    conda:
        "../envs/deeptools.yml"
    shell:
        """ # before it was reference-point
        computeMatrix scale-regions --regionBodyLength 5000 -a 10000 -b 10000 -R {input.rbed} -S {params.bwList_input} -o {output} &> {log}
        """
        
rule ComputeMatrixForMetaplotMetagene:
    input:
        expand('NonCodingRNA/QTLs/metaplots/{{QTLs}}.{{var}}.{{phe}}/Summarised-chRNA-{strand}-{Size}.bw',
                strand = ['+', '-'], Size = ['High', 'Mid', 'Low']),
        rbed = GetMetaplotBed,
        bw_list = expand(
            'NonCodingRNA/QTLs/metaplots/{{QTLs}}.{{var}}.{{phe}}/Summarised-{Phenotype}-.-{Size}.bw', 
            Phenotype = ['polyA.RNA', 'H3K27AC', 'H3K4ME3', 'H3K36ME3', 'ProCap', 'CTCF', 'DNaseISensitivity'], 
            Size = ['High', 'Mid', 'Low']
        ),
    output:
        "NonCodingRNA/QTLs/metaplots/mat/{QTLs}.{var}.{phe}_metagene.mat.gz"
    wildcard_constraints:
        QTLs = 'apQTLs|diQTLs',
        var = 'ncQTL|eQTL',
        phe = 'uaGene|uaRNA'
    params:
        bwList_input = GetBwListForCommand
    log:
        "logs/NonCodingRNA/metaplots/ComputeMatrix.{QTLs}.{var}.{phe}_metagene.log"
    resources:
        mem_mb = 36000
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        computeMatrix reference-point --regionBodyLength 5000 -a 10000 -b 10000 -R {input.rbed} --metagene -S {params.bwList_input} -o {output} &> {log}
        """
    
    
    
    
    
    