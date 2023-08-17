rule ExtractJuncs:
    input:
        bam = GetBamForBigwig,
        bai = GetBaiForBigwig
    output:
        junc = temp(expand("SplicingAnalysis/leafcutter/juncfiles/chr{chrom}/{{Phenotype}}_{{IndID}}_{{Rep}}.junc", chrom=autosomes)),
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",
    params:
        # strand = GetLibStrandForRegtools
        strand = 0
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{Phenotype}/{IndID}.{Rep}.log"
    shell:
        """
        for chrom in {autosomes}
        do
            (regtools junctions extract -m 20 -s {params.strand} -r chr${{chrom}} {input.bam} > SplicingAnalysis/leafcutter/juncfiles/chr${{chrom}}/{wildcards.Phenotype}_{wildcards.IndID}_{wildcards.Rep}.junc ) &> {log}
        done
        cat {output.junc} > {output.junc_autosomes}
        """
        
rule ExtractJuncs_NMD_KD:
    input:
        bam = "Alignments/STAR_Align/NMD_KD/{Phenotype}/{IndID}/1/Filtered.bam",
        bai = "Alignments/STAR_Align/NMD_KD/{Phenotype}/{IndID}/1/Filtered.bam.bai", 
    output:
        junc = temp(expand("SplicingAnalysis/NMD_KD/juncfiles/chr{chrom}/{{Phenotype}}_{{IndID}}_{{Rep}}.junc", chrom=autosomes)),
        junc_autosomes = "SplicingAnalysis/NMD_KD/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",
    wildcard_constraints:
        Phenotype = "HeLa.scr|HeLa.UPF1.KD|HeLa.SMG6.KD|HeLa.SMG7.KD|HeLa.dKD",
        IndID = '|'.join(NMD_KD_accessions),
        Rep = '1'
    params:
        strand = 0
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs_NMD_KD/{Phenotype}/{IndID}.{Rep}.log"
    shell:
        """
        for chrom in {autosomes}
        do
            (regtools junctions extract -m 20 -s {params.strand} -r chr${{chrom}} {input.bam} > SplicingAnalysis/NMD_KD/juncfiles/chr${{chrom}}/{wildcards.Phenotype}_{wildcards.IndID}_{wildcards.Rep}.junc ) &> {log}
        done
        cat {output.junc} > {output.junc_autosomes}
        """

rule make_leafcutter_juncfile:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
    output:
        "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    params:
        SamplesToRemove = ""
    run:
        import os
        if params.SamplesToRemove:
            SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        else:
            SamplesToRemove=[]
        with open(output[0], "w") as out:
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

rule MakeLongJuncTable:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
    output:
        "SplicingAnalysis/CombinedJuncTables/All.tsv.gz"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{print "chrom", "start", "stop", "strand","Dataset", "IndID", "RepNumber", "Count"}} {{n=split(FILENAME, f, "/"); split(f[n], b, "_"); split($11, a, ","); print $1, $2+a[1], $3-a[2], $6, b[1],b[2],b[3], $5}}' {input} | gzip > {output}
        """

use rule MakeLongJuncTable as MakeLongJuncTable_NMD_KD with:
    input:
        expand("SplicingAnalysis/NMD_KD/juncfiles/autosomes/HeLa.scr_{IndID}_1.junc", IndID=['SRR4081222', 'SRR4081223',
        'SRR4081224', 'SRR4081237', 'SRR4081238', 'SRR4081239']),
        expand("SplicingAnalysis/NMD_KD/juncfiles/autosomes/HeLa.SMG6.KD_{IndID}_1.junc", IndID=['SRR4081231', 'SRR4081232', 'SRR4081233']),
        expand("SplicingAnalysis/NMD_KD/juncfiles/autosomes/HeLa.SMG7.KD_{IndID}_1.junc", IndID=['SRR4081240', 'SRR4081241', 'SRR4081242']),
        expand("SplicingAnalysis/NMD_KD/juncfiles/autosomes/HeLa.UPF1.KD_{IndID}_1.junc", IndID=['SRR4081225', 'SRR4081226', 'SRR4081227']),
        expand("SplicingAnalysis/NMD_KD/juncfiles/autosomes/HeLa.dKD_{IndID}_1.junc", IndID=['SRR4081246', 'SRR4081247', 'SRR4081248']),
    output:
        "SplicingAnalysis/CombinedJuncTables/NMD_KD.tsv.gz"

use rule MakeLongJuncTable as MakeLongJuncTable_YRI with:
    input:
        expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df_YRI['Phenotype'], IndID=RNASeqSamplesNoProcap_df_YRI['IndID'], Rep=RNASeqSamplesNoProcap_df_YRI['RepNumber']),
    output:
        "SplicingAnalysis/CombinedJuncTables/YRI.tsv.gz"

rule leafcutter_cluster:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    params:
        rundir = "-r SplicingAnalysis/leafcutter/clustering/autosomes/"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_cluster/autosomes.log"
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input.juncfile_list} {params.rundir} &> {log}
        """

def GetAnnotationTypeGtf(wildcards):
    if wildcards.AnnotationType == "basic":
        return "ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf"
    elif wildcards.AnnotationType == "comprehensive":
        return "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"

rule annotate_juncfiles:
    input:
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
        junc_autosomes = "SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",
        gtf = GetAnnotationTypeGtf
    output:
        "SplicingAnalysis/leafcutter/regtools_annotate/{AnnotationType}/{Phenotype}_{IndID}_{Rep}.bed.gz",
    wildcard_constraints:
        AnnotationType = "basic|comprehensive"
    conda:
        "../envs/regtools.yml"
    log:
        "logs/annotate_juncfiles/{AnnotationType}/{Phenotype}_{IndID}_{Rep}.log"
    shell:
        """
        (regtools junctions annotate {input.junc_autosomes} {input.fa} {input.gtf} | gzip - > {output} ) &> {log}
        """

rule ConcatUniqJuncs:
    input:
        regtools_annotate = expand ("SplicingAnalysis/leafcutter/regtools_annotate/{{AnnotationType}}/{Phenotype}_{IndID}_{Rep}.bed.gz",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
    output:
        "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bed.gz"
    wildcard_constraints:
        AnnotationType = "basic|comprehensive"
    log:
        "logs/ConcatUniqJuncs/{AnnotationType}.log"
    resources:
        mem = much_more_mem_after_first_attempt
    shell:
        """
        (cat <(printf "chrom\\tstart\\tend\\tstrand\\tsplice_site\\tacceptors_skipped\\texons_skipped\\tdonors_skipped\\tanchor\\tknown_donor\\tknown_acceptor\\tknown_junction\\tgene_names\\tgene_id\\n") <(zcat {input.regtools_annotate} | awk -v OFS='\\t' '$1 != "chrom" {{print $1,$2,$3, $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' | sort | uniq ) | gzip - > {output}) &> {log}
        """

rule NMD_transcript_tag_introns:
    input:
        "ReferenceGenome/Annotations/gencode.v34.chromasomal.annotation.gtf"
    output:
        NMD = "SplicingAnalysis/Annotations/NMD/NMD_trancsript_introns.bed.gz",
        NonNMD = "SplicingAnalysis/Annotations/NMD/NonNMD_trancsript_introns.bed.gz"
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        grep 'nonsense_mediated_decay' {input} | bedparse gtf2bed /dev/stdin | awk '{{OFS="\\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {{A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";}} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}}' | bed12ToBed6 -i /dev/stdin | awk -F'\\t' -v OFS='\\t' '{{print $4=""; print $0}}' | sort | uniq | gzip - > {output.NMD}
        grep -v 'nonsense_mediated_decay' {input} | bedparse gtf2bed /dev/stdin | awk '{{OFS="\\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {{A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";}} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}}' | bed12ToBed6 -i /dev/stdin | awk -F'\\t' -v OFS='\\t' '{{print $4=""; print $0}}' | sort | uniq | gzip - > {output.NonNMD}

        """

rule GatherConcatUniqJuncs:
    input:
        expand("SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bed.gz", AnnotationType = ["basic", "comprehensive"])

rule IntronsToGenes:
    input:
        expand("SplicingAnalysis/regtools_annotate_combined/comprehensive.bed.gz", AnnotationType = ["basic", "comprehensive"])

rule GetObservedSpliceSites:
    input:
        bed = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bed.gz",
        fa = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa",
    output:
        donors = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.5ss.bed.gz",
        acceptors = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.3ss.bed.gz",
        bptregions = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bptregion.bed.gz",
        donors_tbi = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.5ss.bed.gz.tbi",
        acceptors_tbi = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.3ss.bed.gz.tbi",
        bptregions_tbi = "SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bptregion.bed.gz.tbi"
    shell:
        """
        zcat {input.bed} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,".", $10,$4}}' | awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1, $2-3, $2+7, ".", $5, $6}} $6=="-" {{print $1, $3-8, $3+2, ".", $5, $6}}' | bedtools getfasta -s -bedOut -bed - -fi {input.fa} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3, "SpliceDonor_"$5, $1"_"$2"_"$3"_"$6"_"$7, $6}}' | sort | uniq | bedtools sort -i - | bgzip /dev/stdin -c > {output.donors}
        zcat {input.bed} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,".", $11,$4}}' | awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1, $3-11, $3-1, ".", $5, $6}} $6=="-" {{print $1, $2, $2+10, ".", $5, $6}}' | bedtools getfasta -s -bedOut -bed - -fi {input.fa} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3, "SpliceAcceptor_"$5, $1"_"$2"_"$3"_"$6"_"$7, $6}}' | sort | uniq | bedtools sort -i - | bgzip /dev/stdin -c > {output.acceptors}
        zcat {input.bed} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,".", $11,$4}}' | awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1, $3-41, $3-11, ".", $5, $6}} $6=="-" {{print $1, $2+10, $2+40, ".", $5, $6}}' | bedtools getfasta -s -bedOut -bed - -fi {input.fa} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3, "SpliceBranchpointRegion_"$5, $1"_"$2"_"$3"_"$6"_"$7, $6}}' | sort | uniq | bedtools sort -i - | bgzip /dev/stdin -c > {output.bptregions}
        tabix -p bed {output.donors}
        tabix -p bed {output.acceptors}
        tabix -p bed {output.bptregions}
        """

rule GetIntronFeatures:
    """
    independent introns that don't overlap any annotated exons, and are in
    expressed genes to quantify intron retention phenotypes
    """
    input:
        AllIntron = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.introns.bed",
        IndependentIntrons = "ReferenceGenome/Annotations/GTFTools/gencode.v34.chromasomal.independentintron.bed",
        eQTL_genes_qqnorm = "QTLs/QTLTools/Expression.Splicing/OnlyFirstReps.qqnorm.bed.gz",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",

    log:
        "logs/GetIntronFeatures.log"
    output:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.saf"
    shell:
        """
        (bedtools sort -i {input.AllIntron} | bedtools intersect -a {input.IndependentIntrons} -b - -f 1 -r  -sorted -wo | awk -F'\\t' '$4==$10 {{print $1, $2, $3, $4, $5, $11}}' | sort | uniq | grep -F -f <(zcat {input.eQTL_genes_qqnorm} | awk -F'\\t' '{{print $4}}') - | awk -v OFS='\\t' '{{$1="chr"$1; print $0}}' |  bedtools slop -b -3 -i - -g {input.fai} | awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }} {{print $4"_IntID."NR, $1, $2, $3, $6}}' > {output} ) &> {log}
        """

rule GetIntronRetentionSpliceSitesAndIntronsBed:
    input:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.saf"
    output:
        "SplicingAnalysis/Annotations/IntronRetentionTargets.bed"
    shell:
        """
        awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $2, $3,$4,$1, ".", $5}}' {input} > {output}
        """

rule CountSpliceSitesOverIR_Intron_Features:
    input:
        IRFeatures = "SplicingAnalysis/Annotations/IntronRetentionTargets.bed",
        fai = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai",
        vcf = "QTLs/QTLTools/Expression.Splicing.Subset_YRI/Genotypes/WholeGenome.vcf.gz"
    output:
        ss_bed = "SplicingAnalysis/Annotations/IntronRetentionTargets_ss.bed",
        vcf = "SplicingAnalysis/Annotations/YRI.Snps.OverInts.vcf.gz"
    shell:
        """
        bedtools slop -b 1 -i {input.IRFeatures} -g {input.fai} | bedtools flank -b 2 -g {input.fai} -i - | bedtools sort -i - > {output.ss_bed}
        bcftools view -O z -R {output.ss_bed} {input.vcf} > {output.vcf}
        """

rule featureCounts_IR:
    input:
        bam = GetBamForPhenotype,
        annotations = "SplicingAnalysis/Annotations/IntronRetentionTargets.saf",
    output:
        "SplicingAnalysis/IR/{Phenotype}/Counts.txt"
    threads:
        8
    params:
        extraParams = GetFeatureCountsParams,
        pairedEndParams = PairedEndParams
    resources:
        mem = 12000,
        cpus_per_node = 9
    log:
        "logs/featureCounts_IR/{Phenotype}.log"
    shell:
        """
        featureCounts {params.pairedEndParams} {params.extraParams} -F SAF -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

def Get_intron_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "SplicingAnalysis/IR/chRNA.Expression.Splicing/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "SplicingAnalysis/IR/Expression.Splicing/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.30min.IR":
        return "SplicingAnalysis/IR/MetabolicLabelled.30min/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.IR":
        return "SplicingAnalysis/IR/MetabolicLabelled.60min/Counts.txt"

def Get_gene_feature_Counts_ForIR(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "featureCounts/chRNA.Expression/Counts.txt"
    elif wildcards.Phenotype == "polyA.IR":
        return "featureCounts/polyA.Expression/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.30min.IR":
        return "featureCounts/MetabolicLabelled.30min/Counts.txt"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.IR":
        return "featureCounts/MetabolicLabelled.60min/Counts.txt"

def SampleMinimumIR_Counts(wildcards):
    if wildcards.Phenotype == "chRNA.IR":
        return "100000"
    elif wildcards.Phenotype == "polyA.IR":
        return "1"

rule featureCounts_IR_to_bedgz:
    input:
        IR_counts = Get_intron_feature_Counts_ForIR,
        Gene_counts = Get_gene_feature_Counts_ForIR
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
        plot = "SplicingAnalysis/IR/{Phenotype}/IR.Ratio.plot.pdf",
        ir = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.IR.bed.gz",
    wildcard_constraints:
        Phenotype = "|".join(["chRNA.IR", "polyA.IR", "MetabolicLabelled.30min.IR", "MetabolicLabelled.60min.IR"])
    log:
        "logs/featureCounts_IR_to_bedgz/{Phenotype}.log"
    params:
        SampleMinimumIR_Counts = "1" #SampleMinimumIR_Counts
    resources:
        mem_mb = 32000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/ProcessIRFeatureCounts.R {input.IR_counts} {input.Gene_counts} {output.bed} {output.plot} {output.ir} {params.SampleMinimumIR_Counts} &> {log}
        """

rule SplitLeafcutter_countsTable:
    input:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz"
    output:
        expand("SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.{Phenotype}.gz", Phenotype=[p for p in RNASeqPhenotypes if p != 'ProCap'])
    log:
        "logs/SplitLeafcutter_countsTable.log"
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/SplitLeafcutterPerindCounts.R {input} &> {log}
        """

def GetSplitLeafcutterCountTablesForPhenotype(wildcards):
    if wildcards.Phenotype == "polyA.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.Expression.Splicing.gz"
    elif wildcards.Phenotype == "chRNA.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.chRNA.Expression.Splicing.gz"
    elif wildcards.Phenotype == "MetabolicLabelled.30min.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.MetabolicLabelled.30min.gz"
    elif wildcards.Phenotype == "MetabolicLabelled.60min.Splicing":
        return "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz.MetabolicLabelled.60min.gz"
        

rule Leafcutter_countsTable_toPSI:
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        temp("QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed")
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
    log:
        "logs/Leafcutter_countsTable_toPSI/{Phenotype}.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/LeafcutterToPSIBed.R {input} {output} &> {log}
        """

rule sortAndIndexPSITable:
    input:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed"
    output:
        bed = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed.gz",
        tbi = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.PSI.bed.gz.tbi"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
    shell:
        """
        bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule leafcutter_PreparePhenotypes:
    """
    scripts in main leafcutter repo is buggy. Use fork from https://github.com/mdshw5/leafcutter.git
    """
    input:
        GetSplitLeafcutterCountTablesForPhenotype
    output:
        "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "chRNA.Splicing", "MetabolicLabelled.30min.Splicing", "MetabolicLabelled.60min.Splicing"])
    conda:
        "../envs/py27.yaml"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_PreparePhenotypes/{Phenotype}.log"
    shadow: "shallow"
    shell:
        """
        python2.7 scripts/prepare_phenotype_table.py {input} -p 20 &> {log}
        (awk -F'\\t' -v OFS='\\t' 'NR==1 {{$4="pid\\tgid\\tstrand"; print $0}} FNR!=1 {{$1="chr"$1; split($4, a, ":"); split(a[4], b, "_"); $4=$4"\\t"$1"_"a[4]"\\t"b[3]; print $0}}' {input}.qqnorm_chr* | gzip - > {output} ) &>> {log}
        """

rule Subset_YRI_phenotype_table:
    input:
        input_file = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
        igsr = '../data/igsr_samples.tsv.gz'
    output:
        "QTLs/QTLTools/{Phenotype}.Subset_YRI/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "polyA.IR", "polyA.IER", "polyA.Splicing.5PrimeSS", "polyA.Splicing.3PrimeSS"])
    log:
        "logs/Subsample_YRI.{Phenotype}.log"
    shell:
        """
        python scripts/subsample.Splicing_YRI.py --input {input.input_file} --output {output} &> {log}
        """
        
rule Subset_EUR_phenotype_table:
    input:
        input_file = "QTLs/QTLTools/{Phenotype}/OnlyFirstReps.qqnorm.bed.gz",
        igsr = '../data/igsr_samples.tsv.gz'
    output:
        "QTLs/QTLTools/{Phenotype}.Subset_EUR/OnlyFirstReps.qqnorm.bed.gz"
    wildcard_constraints:
        Phenotype = "|".join(["polyA.Splicing", "polyA.IR", "polyA.IER", "polyA.Splicing.5PrimeSS", "polyA.Splicing.3PrimeSS"])
    log:
        "logs/Subsample_YRI.{Phenotype}.log"
    shell:
        """
        python scripts/subsample.Splicing_EUR.py --input {input.input_file} --output {output} &> {log}
        """

rule MakeNormalizedPsiTables:
    input:
        numers = "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    output:
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"]),
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.{Phenotype}.bed", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"])
    log:
        "logs/MakeNormalizedPsiTables.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/MakeNormalizedPSI.Tables.R SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI. {input.numers} &> {log}
        """

rule BgzipAndTabixPsiTables:
    input:
        "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed"
    output:
        bed = "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz",
        tbi = "SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz.tbi"
    log:
        "logs/BgzipAndTabixPsiTables/{Phenotype}.log"
    shell:
        """
        bgzip {input}
        tabix -p bed {output.bed}
        """

rule MakeMorePSI_Tables:
    input:
        juncs = "SplicingAnalysis/CombinedJuncTables/All.tsv.gz",
        Annotations = "../data/IntronAnnotationsFromYang.Updated.tsv.gz",
        clustered_juncs = "../code/SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    output:
        expand("SplicingAnalysis/NormalizedPsiTables/{PsiStat}_{Dataset}.bed", PsiStat = ["CountsOverClusterSum", "CountsOverClusterMax", "CountsOverMedianProductive", "CountsOverMedianProductiveCoercedRange"], Dataset = ["chRNA.Expression.Splicing", "Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])
    conda:
        "../envs/r_2.yaml"
    resources:
        mem_mb = 48000
    log:
        "logs/MakeMorePSI_Tables.log"
    shell:
        """
        (Rscript scripts/MakeNormalizedPSI.Tables.ManyWays.R SplicingAnalysis/NormalizedPsiTables/ {input.juncs} {input.Annotations} ) &> {log}
        """

use rule BgzipAndTabixPsiTables as BgzipAndTabixMorePsiTables with:
    input:
        "SplicingAnalysis/NormalizedPsiTables/{PsiStat}_{Dataset}.bed" 
    log:
        "logs/BgzipAndTabixMorePsiTables/{PsiStat}/{Dataset}.log"
    output:
        bed = "SplicingAnalysis/NormalizedPsiTables/{PsiStat}_{Dataset}.bed.gz",
        tbi = "SplicingAnalysis/NormalizedPsiTables/{PsiStat}_{Dataset}.bed.gz.tbi"

rule GatherBgzipedTabixPsiTables:
    input:
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.{Phenotype}.bed.gz", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"]),
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed.gz", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"]),
        expand("SplicingAnalysis/NormalizedPsiTables/{PsiStat}_{Dataset}.bed.gz.tbi", PsiStat = ["CountsOverClusterSum", "CountsOverClusterMax", "CountsOverMedianProductive", "CountsOverMedianProductiveCoercedRange"], Dataset = ["chRNA.Expression.Splicing", "Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min"])

rule SumJuncCountsAcrossYRISamples:
    input:
        junc_tables = expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.JunctionCounts.{Phenotype}.bed.gz", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"]),
        juncs_list = "SplicingAnalysis/regtools_annotate_combined/comprehensive.bed.gz"
    output:
        Summarised = "SplicingAnalysis/regtools_annotate_combined/Comprehensive.YRI.Sample.Sum.Counts.tsv.gz",
        LongTable = "SplicingAnalysis/regtools_annotate_combined/Comprehensive.YRI.Sample.LongTable.Counts.tsv.gz"
    log:
        "logs/SumJuncCountsAcrossYRISamples.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/SumJuncCountsAcrossYRISamples.R {output.Summarised} {output.LongTable} {input.juncs_list} {input.junc_tables} &> {log}
        """

rule GetPSI_of_GWAS_sQTLs:
    input:
        Gwas_sQTLs = "../output/sQTLsThatColocWithGWAS.tsv.gz",
        PSI = "SplicingAnalysis/NormalizedPsiTables/CountsOverClusterMax_{RNAseqDataset}.bed.gz",
        tbi = "SplicingAnalysis/NormalizedPsiTables/CountsOverClusterMax_{RNAseqDataset}.bed.gz.tbi",
        VCF = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz",
        VCF_tbi = "Genotypes/1KG_GRCh38/Autosomes.vcf.gz.tbi"
    output:
        "SplicingAnalysis/GWAS_sQTL_PSI/{RNAseqDataset}.tsv.gz"
    log:
        "logs/GetPSI_of_GWAS_sQTLs/{RNAseqDataset}.log"
    conda:
        "../scripts/GenometracksByGenotype/GenometracksByGenotype.conda_env.yml"
    shell:
        """
        python scripts/GetPSI_For_Gwas_sQTLs.py {input.Gwas_sQTLs} {input.PSI} {input.VCF} {output} &> {log}
        """

rule GatherGetPSI_of_GWAS_sQTLs:
    input:
        expand("SplicingAnalysis/GWAS_sQTL_PSI/{RNAseqDataset}.tsv.gz", RNAseqDataset = ["Expression.Splicing", "chRNA.Expression.Splicing"])





###########################################
#          Out of order splicing          #
###########################################

rule QQnormOutOfOrder:
    input:
        "/project2/yangili1/yangili/chRNA_order/splicing_order/splice_order/table.txt"
    output:
        "QTLs/QTLTools/chRNA.Splicing.Order/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/qqnorm_splicing_order.log"
    resources:
        mem_mb = 58000
    params:
        min_obs = 35
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/PrepareOutOfOrderQQnorm.py --input_file {input} --output {output} --min_obs {params.min_obs} &> {log}
        """
    

rule QQnormRNAEditing:
    input:
        "/project2/yangili1/cdai/aicher/code/chRNA/Results/hs38/GatherEditing/EL.txt",
    output:
        "QTLs/QTLTools/chRNA.RNA.Editing/OnlyFirstReps.qqnorm.bed.gz"
    log:
        "logs/qqnorm_rna_editing.log"
    resources:
        mem_mb = 58000
    params:
        min_obs = 35
    conda:
        "../envs/py_tools.yml"
    shell:
        """
        python scripts/PrepareRNAEditingQQnorm.py --input_file {input} --output {output} --min_obs {params.min_obs} &> {log}
        """



