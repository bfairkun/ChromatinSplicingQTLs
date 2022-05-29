        
rule GetSpliceSitePhenotypes:
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        "QTLs/QTLTools/{Phenotype}.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/{Phenotype}.5PrimeSS/OnlyFirstReps.PSI.bed.gz",
        "QTLs/QTLTools/{Phenotype}.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/{Phenotype}.3PrimeSS/OnlyFirstReps.PSI.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing"])
    log:
        "logs/leafcutter_PreparePhenotypes/{Phenotype}.SpliceSites.log"
    params:
        flags = GetSkippedSamples
    conda:
        "../envs/py_tools.yml"
    resources:
        mem = 32000,
    shell:
        """
        python scripts/GetSpliceSitesFromLeafcutter.py --counts {input} --output QTLs/QTLTools/{wildcards.Phenotype} {params.flags} &> {log}
        """
    
    
# Ugly rule, but merging quick processes to avoid clogging slurm
rule MakeSpliceSiteAnnotation:
    input:
        "QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz",
        "QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz"
    output:
        "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed.gz",
        "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed.gz"
    params:
        up = "494",
        dn = "497"
        #up = "94",
        #dn = "97"
    resources:
        mem = 12000,
    shell:
        """
        zcat QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' > ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ; 
        zcat QTLs/QTLTools/chRNA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ;
        zcat QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ; 
        zcat QTLs/QTLTools/polyA.Splicing.5PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed ;
        sort -u ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed > ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed;
        gzip ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed;
        rm ReferenceGenome/Annotations/SpliceSites.5PrimeSS_.bed;
        zcat QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' > ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ; 
        zcat QTLs/QTLTools/chRNA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ;
        zcat QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="-" {{print $1"\\t"$2-{params.up}"\\t"$3+{params.dn}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ; 
        zcat QTLs/QTLTools/polyA.Splicing.3PrimeSS/OnlyFirstReps.qqnorm.bed.gz | awk -F'\\t' '$6=="+" {{print $1"\\t"$2-{params.dn}"\\t"$3+{params.up}"\\t"$4"\\t"$5"\\t"$6"\\t"$2"\\t"$3}}' >> ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed ;
        sort -u ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed > ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed;
        gzip ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed;
        rm ReferenceGenome/Annotations/SpliceSites.3PrimeSS_.bed;
        """
        
rule MakeChromatinSpliceSitePeakPhenotypes:
    input:
        SpliceSites = "ReferenceGenome/Annotations/SpliceSites.{Prime}.bed.gz",
        Chromatin = "QTLs/QTLTools/{chromatinPhenotype}/OnlyFirstReps.sorted.qqnorm.bed.gz"
    output:
        "QTLs/QTLTools/{chromatinPhenotype}.{Prime}.Peaks/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'H3K4ME1|H3K4ME3|H3K27AC'
    params:
        start = "85",
        end = "86",
        name = "82",
        ncols = "78"
    resources:
        mem = 52000,
    shell:
        """
        zcat {input.Chromatin} | head -n 1 > QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed;
        bedtools intersect -wo -F 0.1 -a {input.Chromatin} -b {input.SpliceSites} | awk -F"\\t" '{{printf($1"\\t"${params.start}"\\t"${params.end}"\\t"${params.name}"\\t"); for (x=5; x<={params.ncols}; x++) printf("%s\\t", $x);printf("\\n"); }}' >> QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed;
        gzip QTLs/QTLTools/{wildcards.chromatinPhenotype}.{wildcards.Prime}.Peaks/OnlyFirstReps.qqnorm.bed
        """
    
use rule MakeChromatinSpliceSitePeakPhenotypes as MakeChromatinSpliceSitePeakH3K36ME3 with:
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'H3K36ME3'
    params:
        start = "107",
        end = "108",
        name = "104",
        ncols = "100"

    
use rule MakeChromatinSpliceSitePeakPhenotypes as MakeChromatinSpliceSiteCTCF with:
    wildcard_constraints:
        Prime = '5PrimeSS|3PrimeSS',
        chromatinPhenotype = 'CTCF'
    params:
        start = "64",
        end = "65",
        name = "61",
        ncols = "57"


rule MakeSpliceSiteSAF:
    input:
        in5 = "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.bed.gz",
        in3 = "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.bed.gz"
    output:
        out5 = "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.saf",
        out3 = "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.saf",
    shell:
        """
        zcat {input.in5} | awk -F'\\t' '{{print $4"\\t"$1"\\t"$2"\\t"$3"\\t."}}' > {output.out5};
        zcat {input.in3} | awk -F'\\t' '{{print $4"\\t"$1"\\t"$2"\\t"$3"\\t."}}' > {output.out3};
        """

def  GetAnnotationsForSpliceSite(wildcards):
    if wildcards.SpliceSite == '5PrimeSS':
        return "ReferenceGenome/Annotations/SpliceSites.5PrimeSS.saf"
    elif wildcards.SpliceSite == '3PrimeSS':
        return "ReferenceGenome/Annotations/SpliceSites.3PrimeSS.saf"

rule featureCountsForSpliceSite:
    input:
        bam = GetBamForPhenotype,
        annotations = GetAnnotationsForSpliceSite
    output:
        "featureCounts/{Phenotype}.{SpliceSite}/Counts.txt"
    params:
        extraParams = GetFeatureCountsParams,
        paired = PairedEndParams
    threads:
        8
    wildcard_constraints:
        Phenotype = 'H3K4ME1|H3K4ME3|H3K27AC|H3K36ME3',
        SpliceSite = '5PrimeSS|3PrimeSS'
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.{SpliceSite}.log"
    shell:
        """
        featureCounts {params.paired} {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """
