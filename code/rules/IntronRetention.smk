def GetSkippedSamples(wildcards):
    if wildcards.Phenotype == 'chRNA.Splicing':
        return '--skip_samples NA18855'
    else:
        return '--skip_samples NOSAMPLE'
        
    
def GetSPLICEqModeParams(wildcards):
    if wildcards.mode == "IER":
        return "--IERatio"
    else:
        return "-f 3"

rule RunSPLICEq_chRNA:
    input:
        bam = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/chRNA.Expression.Splicing/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/chRNA.{mode}/{IndID}.spliceq.tab"
    params:
        GetSPLICEqModeParams
    conda:
        "../envs/spliceq.yml"
    wildcard_constraints:
        IndID = '|'.join(chRNAsamples),
        mode = 'IER' 
    log:
        "logs/spliceq/chRNA.{mode}/{IndID}.log"
    shell:
        """
        SPLICE-q.py -b {input.bam} -g {input.Comprehensive_gtf} -o {output} {params} --MinCoverage 0 &> {log}
        """

use rule RunSPLICEq_chRNA as rule RunSPLICEq_polyA with:
    input:
        bam = 'Alignments/STAR_Align/Expression.Splicing/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/Expression.Splicing/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    resources:
        mem_mb = 50000
    output:
        "SplicingAnalysis/SPLICEq/polyA.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/polyA.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(polyAsamples),
        mode = 'IER'
        
use rule RunSPLICEq_chRNA as rule RunSPLICEq_M30 with:
    input:
        bam = 'Alignments/STAR_Align/MetabolicLabelled.30min/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/MetabolicLabelled.30min/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/MetabolicLabelled.30min.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/MetabolicLabelled.30min.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(metabolic30samples),
        mode = 'IER'
        
use rule RunSPLICEq_chRNA as rule RunSPLICEq_M60 with:
    input:
        bam = 'Alignments/STAR_Align/MetabolicLabelled.60min/{IndID}/1/Filtered.bam',
        bai = 'Alignments/STAR_Align/MetabolicLabelled.60min/{IndID}/1/Filtered.bam.bai',
        Comprehensive_gtf = "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf",
    output:
        "SplicingAnalysis/SPLICEq/MetabolicLabelled.60min.{mode}/{IndID}.spliceq.tab"
    log:
        "logs/spliceq/MetabolicLabelled.60min.{mode}/{IndID}.log"
    wildcard_constraints:
        IndID = '|'.join(metabolic60samples),
        mode = 'IER'
        
rule CollectSPLICEq_chRNA:
    input: 
        expand("SplicingAnalysis/SPLICEq/chRNA.{{mode}}/{IndID}.spliceq.tab", IndID=chRNAsamples)
    output:
        out_qqnorm = "QTLs/QTLTools/{Phenotype}.{mode}/OnlyFirstReps.qqnorm.bed.gz",
        out = "QTLs/QTLTools/{Phenotype}.{mode}/OnlyFirstReps.ratios.bed.gz",
    log:
        "logs/spliceq/{Phenotype}.{mode}.log"
    wildcard_constraints:
        mode = 'IER',
        Phenotype = 'chRNA'
    resources:
        mem_mb = 32000
    params:
        "NA18855"
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/collect_spliceq.py --spliceq_dir SplicingAnalysis/SPLICEq/{wildcards.Phenotype}.{wildcards.mode}/  --output {output.out} --output_qqnorm {output.out_qqnorm} --skip_samples {params} &> {log}
        """

use rule CollectSPLICEq_chRNA as CollectSPLICEq_polyA with:
    input: 
        expand("SplicingAnalysis/SPLICEq/polyA.{{mode}}/{IndID}.spliceq.tab", IndID=polyAsamples)
    wildcard_constraints:
        mode = 'IER',
        Phenotype = 'polyA'
    params:
        "None"
        
use rule CollectSPLICEq_chRNA as CollectSPLICEq_M30 with:
    input: 
        expand("SplicingAnalysis/SPLICEq/MetabolicLabelled.30min.{{mode}}/{IndID}.spliceq.tab", IndID=metabolic30samples)
    wildcard_constraints:
        mode = 'IER',
        Phenotype = 'MetabolicLabelled.30min'
    params:
        "None"


use rule CollectSPLICEq_chRNA as CollectSPLICEq_M60 with:
    input: 
        expand("SplicingAnalysis/SPLICEq/MetabolicLabelled.60min.{{mode}}/{IndID}.spliceq.tab", IndID=metabolic60samples)
    wildcard_constraints:
        mode = 'IER',
        Phenotype = 'MetabolicLabelled.60min'
    params:
        "None"


