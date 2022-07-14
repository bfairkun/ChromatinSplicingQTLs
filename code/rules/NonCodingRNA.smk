### Rules to get non-coding RNA

rule get_eRNA_saf:
    input:
        "ProCapAnalysis/CountTable.hg38.features.bed",
    output:
        "../data/eRNA.saf",
        "../data/eRNA_both_strands.saf",
    log:
        "logs/eRNA_reference.log"
    shell:
        """
        python scripts/get_eRNA_saf.py &> {log}
        """

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

rule get_cheRNA_saf:
    output:
        "../data/cheRNA_K562_GSE83531.saf",
    log:
        "logs/cheRNA_reference.log"
    shell:
        """
        python scripts/get_cheRNA_saf.py &> {log}
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

rule featureCountsNonCoding:
    input:
        bam = GetBamForPhenotype,
        eRNA = SAFeRNAForPhenotype,
        cheRNA = "../data/cheRNA_K562_GSE83531.saf",
    output:
        "featureCounts/{Phenotype}_eRNA/Counts.txt",
        "featureCounts/{Phenotype}_cheRNA/Counts.txt",
    params:
        extraParams = PairedEndParams,
        strandParams = FeatureCountsNonCodingStrandParams
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"])
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}"
    shell:
        """
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.eRNA} -o featureCounts/{wildcards.Phenotype}_eRNA/Counts.txt {input.bam} &> {log}.eRNA.log;
        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.cheRNA} -o featureCounts/{wildcards.Phenotype}_cheRNA/Counts.txt {input.bam} &> {log}.cheRNA.log
        """



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
    shell:
        """
        python scripts/GetNonCodingRNAFromFeatureCounts.py --phenotype {wildcards.Phenotype} --ncRNA snoRNA &> {log.snoRNA_log};
        python scripts/GetNonCodingRNAFromFeatureCounts.py --phenotype {wildcards.Phenotype} --ncRNA lncRNA &> {log.lncRNA_log}
        """
        ### Add transcribed_unprocessed_pseudogene, other annotated ncRNAs

#rule MakeSAFForNonCodingRNA:
#    input:
#        "NonCodingRNA_200/annotation/ncRNA_filtered.2states.sorted.bed.gz",
#    output:
#        "NonCodingRNA_200/annotation/ncRNA_filtered.2states.sorted.saf",
#    shell:
#        """
#        echo -e 'GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output};
#        zcat {input} | awk '{{print sep='\\t' "ncRNA_" NR, $1, $2, $3, $6 }}' FS='\\t' OFS='\\t' >> {output}
#        """
        
#rule featureCountsNonCodingRNA:
#    input:
#        bam = GetBamForPhenotype,
#        saf = "NonCodingRNA_200/annotation/ncRNA_filtered.2states.sorted.saf",
#    output:
#        "featureCounts/{Phenotype}_ncRNA/Counts.txt",
#    params:
#        extraParams = PairedEndParams,
#        strandParams = FeatureCountsNonCodingStrandParams
#    threads:
#        8
#    wildcard_constraints:
#        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
#    resources:
#        mem = 12000,
#        cpus_per_node = 9,
#    log:
#        "logs/featureCounts/{Phenotype}_ncRNA.log"
#    shell:
#        """
#        featureCounts {params.extraParams} {params.strandParams} -F SAF -T {threads} --ignoreDup --primary -a {input.saf} -o featureCounts/{wildcards.Phenotype}_ncRNA/Counts.txt {input.bam} &> {log};
#        """

        
        
        
        
        
        
        
        
        



