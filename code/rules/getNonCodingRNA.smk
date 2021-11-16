### Rules to get non-coding RNA

rule get_eRNA_saf:
    input:
        "QTLs/QTLTools/ProCap/OnlyFirstReps.qqnorm.bed.gz",
    output:
        "ReferenceGenome/Annotations/eRNA.saf",
    log:
        "logs/eRNA_reference.log"
    shell:
        """
        python scripts/get_eRNA_saf.py &> {log}
        """

rule featureCountsNonCoding:
    input:
        bam = GetBamForPhenotype,
        eRNA = "ReferenceGenome/Annotations/eRNA.saf",
        cheRNA = "ReferenceGenome/Annotations/cheRNA_K562_GSE83531.saf",
    output:
        "featureCounts/{Phenotype}_eRNA/Counts.txt",
        "featureCounts/{Phenotype}_cheRNA/Counts.txt",
    params:
        extraParams = PairedEndParams,
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}"
    shell:
        """
        featureCounts {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.eRNA} -o featureCounts/{wildcards.Phenotype}_eRNA/Counts.txt {input.bam} &> {log}.eRNA.log;
        featureCounts {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.cheRNA} -o featureCounts/{wildcards.Phenotype}_cheRNA/Counts.txt {input.bam} &> {log}.cheRNA.log
        """
      

