### Rules to get non-coding RNA


rule get_introns_saf:
    input:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed",
    output:
        "../data/introns.saf"
    log:
        "logs/introns_reference.log"
    shell:
        """
        python scripts/get_introns_saf.py &> {log}
        """

def FeatureCountsNonCodingStrandParams(wildcards):
    if wildcards.Phenotype == 'chRNA.Expression':
        return "-s 2"
    else:
        return ""
        
def SAFeRNAForPhenotype(wildcards):
    if wildcards.Phenotype == 'chRNA.Expression':
        return "../data/eRNA_both_strands.saf"
    else:
        return "../data/eRNA.saf"



rule featureCountsIntrons:
    input:
        bam = GetBamForPhenotype,
        introns = "../data/introns.saf",
    output:
        "featureCounts/{Phenotype}_introns/Counts.txt",
    params:
        extraParams = PairedEndParams,
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"])
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}_introns.log"
    shell:
        """
        featureCounts {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.introns} -o featureCounts/{wildcards.Phenotype}_introns/Counts.txt {input.bam} &> {log};
        """

 
rule GetAdditionalNonCodingRNAFromFeatureCounts:
    input:
        fCRNA = "featureCounts/{Phenotype}/Counts.txt",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf"
    output:
        "featureCounts/{Phenotype}_snoRNA/Counts.txt",
        "featureCounts/{Phenotype}_lncRNA/Counts.txt"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"]),
        #ncRNA = "|".join(["snoRNA", "lncRNA"])
    log:
        snoRNA_log = "logs/{Phenotype}.get_snoRNA.log",
        lncRNA_log = "logs/{Phenotype}.get_lncRNA.log",
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/GetNonCodingRNAFromFeatureCounts.py --phenotype {wildcards.Phenotype} --ncRNA snoRNA &> {log.snoRNA_log};
        python scripts/GetNonCodingRNAFromFeatureCounts.py --phenotype {wildcards.Phenotype} --ncRNA lncRNA &> {log.lncRNA_log}
        """
        ### Here pseudogenes are included with lncRNAs. There are probably more efficient ways to do this;
        ### Might change later

