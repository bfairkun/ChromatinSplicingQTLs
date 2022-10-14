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

rule leafcutter_cluster:
    input:
        juncs = expand ("SplicingAnalysis/leafcutter/juncfiles/autosomes/{Phenotype}_{IndID}_{Rep}.junc",  zip, Phenotype=RNASeqSamplesNoProcap_df['Phenotype'], IndID=RNASeqSamplesNoProcap_df['IndID'], Rep=RNASeqSamplesNoProcap_df['RepNumber']),
        juncfile_list = "SplicingAnalysis/leafcutter/juncfilelist.autosomes.txt"
    output:
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind.counts.gz",
        "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    resources:
        mem_mb = 16000
    log:
        "logs/leafcutter_cluster/autosomes.log"
    shell:
        """
        python scripts/davidaknowles_leafcutter/scripts/leafcutter_cluster_regtools_py3.py -j {input.juncfile_list} -r SplicingAnalysis/leafcutter/clustering/autosomes/ &> {log}
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
        grep 'nonsense_mediated_decay' {input} | bedparse gtf2bed /dev/stdin --extraFields gene_id,gene_name,transcript_type | gzip - > {output.NMD}
        grep -v 'nonsense_mediated_decay' {input} | bedparse gtf2bed /dev/stdin --extraFields gene_id,gene_name,transcript_type | gzip - > {output.NonNMD}

        """

rule GatherConcatUniqJuncs:
    input:
        expand("SplicingAnalysis/regtools_annotate_combined/{AnnotationType}.bed.gz", AnnotationType = ["basic", "comprehensive"])

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
        
        

rule MakeNormalizedPsiTables:
    input:
        numers = "SplicingAnalysis/leafcutter/clustering/autosomes/leafcutter_perind_numers.counts.gz"
    output:
        expand("SplicingAnalysis/leafcutter/NormalizedPsiTables/PSI.{Phenotype}.bed", Phenotype = ["Expression.Splicing", "MetabolicLabelled.30min", "MetabolicLabelled.60min", "chRNA.Expression.Splicing"])
    log:
        "logs/MakeNormalizedPsiTables.log"
    conda:
        "envs/r_2.yaml"
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
    



