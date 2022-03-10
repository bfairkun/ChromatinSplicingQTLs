### Rules to get non-coding RNA

rule get_eRNA_saf:
    input:
        "ProCapAnalysis/CountTable.hg38.features.bed",
    output:
        "../data/eRNA.saf",
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


rule featureCountsNonCoding:
    input:
        bam = GetBamForPhenotype,
        eRNA = "../data/eRNA.saf",
        cheRNA = "../data/cheRNA_K562_GSE83531.saf",
    output:
        "featureCounts/{Phenotype}_eRNA/Counts.txt",
        "featureCounts/{Phenotype}_cheRNA/Counts.txt",
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
        "logs/featureCounts/{Phenotype}"
    shell:
        """
        featureCounts {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.eRNA} -o featureCounts/{wildcards.Phenotype}_eRNA/Counts.txt {input.bam} &> {log}.eRNA.log;
        featureCounts {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.cheRNA} -o featureCounts/{wildcards.Phenotype}_cheRNA/Counts.txt {input.bam} &> {log}.cheRNA.log
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
        
#rule SubsetYRIncRNA:
#    input:
#        fCRNA = "featureCounts/{Phenotype}_{ncRNA}/Counts.txt",
#        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf"
#    output:
#        "featureCounts/{Phenotype}_{ncRNA}.Subset_YRI/Counts.txt",
#    wildcard_constraints:
#        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"]),
#        ncRNA = "|".join(["cheRNA", "eRNA"])
#    log:
#        "logs/{Phenotype}_{ncRNA}.Subset_YRI.log",
#    shell:
#        """
#        python scripts/Subset_YRI.py --phenotype {wildcards.Phenotype} --ncRNA {wildcards.ncRNA} &> {log};
#        """
 
#rule MergeFeatureCountsAllRNA:
#    input:
#        mRNA = "featureCounts/{Phenotype}/Counts.txt",
#        eRNA = "featureCounts/{Phenotype}_eRNA/Counts.txt",
#        cheRNA = "featureCounts/{Phenotype}_cheRNA/Counts.txt",
#        samples = "config/samples.tsv",
#        igsr_samples = "../data/igsr_samples.tsv.gz"
#    output:
#        "QTLs/QTLTools/{Phenotype}.AllRNA.Subset_YRI/Counts.txt",
#    wildcard_constraints:
#        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "ProCap"])
#    log:
#        "logs/MergeFeatureCountsAllRNA.{Phenotype}.log"
#    shell:
#        """
#        python scripts/MergeFeatureCountsAllRNA.py --phenotype {wildcards.Phenotype} &> {log}
#        """


### Excluding ProCap to avoid confusion

#rule PrepareAllRNACountsForQTLTools:
#    input:
#        featureCounts = "featureCounts/{Phenotype}_{ncRNA}/Counts.txt",
#    output:
#        FirstReps = "QTLs/QTLTools/{Phenotype}_{ncRNA}/OnlyFirstReps.qqnorm.bed.gz"
#    log:
#        "logs/Prepare_{Phenotype}.{ncRNA}_PhenotypeTable.log"
#    conda:
#        "../envs/r_essentials.yml"
#    wildcard_constraints:
#        Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min"]),
#        ncRNA = "|".join(["eRNA", "cheRNA", "snoRNA", "lncRNA"])
#    conda:
#        "../envs/r_essentials.yml"
#    shell:
#        """
#        Rscript scripts/PreparePhenotypeTablesNonCodingRNA.R {input.featureCounts} {output.FirstReps} &> {log}
#        """

rule Prepare_ncRNA_RNA_seq_PhenotypeTable_AndMakeGeneList:
    input:
        featureCounts = "featureCounts/{Phenotype}_{ncRNA}/Counts.txt",
        gtf = "ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        genes = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.genes.bed",
    output:
        GeneList = "ExpressionAnalysis/{Phenotype}_{ncRNA}/ExpressedGeneList.txt",
        PhenotypesBed = "QTLs/QTLTools/{Phenotype}_{ncRNA}/AllReps.qqnorm.bed.gz",
        FirstReps = "QTLs/QTLTools/{Phenotype}_{ncRNA}/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
       Phenotype = "|".join(["polyA.Expression", "chRNA.Expression", "MetabolicLabelled.30min", "MetabolicLabelled.60min"]),
       ncRNA = "|".join(["eRNA", "cheRNA", "snoRNA", "lncRNA"])
    log:
        "logs/Prepare_{Phenotype}_{ncRNA}_RNA_seq_PhenotypeTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PreparePhenotypeTablesFromFeatureCounts.R {input.featureCounts} rename_STAR_alignment_samples {input.gtf} {input.genes} {output.GeneList} {output.PhenotypesBed} {output.FirstReps} &> {log}
        """
