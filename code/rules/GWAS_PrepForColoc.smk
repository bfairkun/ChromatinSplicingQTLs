# rule GWASFromPhoenix:
#     input:
# /project/yangili1/zpmu/GWAS_loci/liftOver/hg38_RA_GWASmeta_TransEthnic_v2.txt.gz
# /project/yangili1/zpmu/GWAS_loci/liftOver/hg38_sle_Vyse_chr1-22.txt.gz
# 

rule MS_GWAS_to_bed:
    input:
        "/project2/yangili1/bjf79/gwas_summary_stats/discovery_metav3.0.meta.gz"
    output:
        "gwas_summary_stats/hg19_summarystat_beds/IMSGC2019.bed.gz"
    shell:
        """
        zcat {input} | awk -v OFS='\\t' 'NR>1 && $8!="NA" {{print "chr"$1, $2, $2+1, $7, $3, $4, $5, $7, $8}}' | gzip - > {output}
        """

rule liftover_Gwas_stats:
    """
    Convert hg19 bed of summary stats to hg38_summarystats. For convenience with downstream rules, make sure the input bed has Pvalues in column4. Other summary stats can be in later columns
    """
    input:
        bed = "gwas_summary_stats/hg19_summarystat_beds/{accession}.bed.gz",
        chain = "ReferenceGenome/Chains/hg19ToHg38.over.chain.gz"
    output:
        temp("gwas_summary_stats/hg38lifted_summarystat_beds/{accession}.bed")
    conda:
        "../envs/crossmap.yml"
    shadow: "shallow"
    wildcard_constraints:
        accession = '|'.join(["IMSGC2019"])
    shell:
        """
        CrossMap.py bed {input.chain} {input.bed} gwas_summary_stats/hg38lifted_summarystat_beds/{wildcards.accession}.bed
        """

rule CreateBedForNonGwasCatalogStats:
    """
    ensure that columns1-3 are standard bed format, and column4 is P value.
    Rest of summary stats can be in later columns
    """
    input:
        "gwas_summary_stats/hg38lifted_summarystat_beds/{accession}.bed"
    wildcard_constraints:
        accession = '|'.join(["IMSGC2019"])
    output:
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
    resources:
        mem_mb = 16000
    shell:
        """
        cat <(printf "#Chr\\tstart\\tend\\tP\\tsnpID\\tA1\\tA2\\tP\\tOR\\n") {input} |  bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

def GetPvalColumnAwk(wildcards):
    if wildcards.accession in ["GCST007800", "GCST007799", "GCST004988"]:
        return "$21"
    elif wildcards.accession in ["GCST004131", "GCST004132", "GCST004133"] :
        return "$17"
    else:
        return "$23"

rule CreateBedForGwasCatalogStats:
    """
    ensure that columns1-3 are standard bed format, and column4 is P value.
    Rest of summary stats can be in later columns
    """
    input:
        summarystats = "gwas_summary_stats/full_data/{accession}.tsv.gz"
    wildcard_constraints:
        accession = '|'.join(gwas_df.loc[~gwas_df['FTPPath'].isna()].index)
    log:
        "logs/CreateBedForGwasCatalogStats/{accession}.log"
    output:
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
    resources:
        mem_mb = 16000
    params:
        GetPvalColumnAwk = GetPvalColumnAwk
    shell:
        """
        (zcat {input.summarystats} | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "#Chr", "start", "end", "P", $0 }} NR>1 && $4!="NA" {{print "chr"$3, $4, $4+1, {params.GetPvalColumnAwk}, $0}}' | sort -k 1,1 -k2,2n  | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

rule CreateBedForGwasCatalogStats_RA:
    input:
        summarystats = lambda wildcards: gwas_df.loc[wildcards.accession]['SummaryStatsLocalFilepath']
    wildcard_constraints:
        accession = "TRANSETHNICRA"
    output:
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
    log:
        "logs/CreateBedForGwasCatalogStats/{accession}.log"
    shell:
        """
        (zcat {input.summarystats} | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "#Chr", "start", "end", "P", $0 }} NR>1 {{print "chr"$10, $11, $11+1, $9, $0}}' | sort -k 1,1 -k2,2n  | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

rule CreateBedForGwasCatalogStats_SLE:
    input:
        summarystats = lambda wildcards: gwas_df.loc[wildcards.accession]['SummaryStatsLocalFilepath']
    wildcard_constraints:
        accession = "SLEVYSE"
    output:
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
    log:
        "logs/CreateBedForGwasCatalogStats/{accession}.log"
    shell:
        """
        (zcat {input.summarystats} | awk -F'\\t' -v OFS='\\t' 'NR==1 {{ print "#Chr", "start", "end", "P", $0 }} NR>1 {{print "chr"$24, $25, $25+1, $12, $0}}' | sort -k 1,1 -k2,2n  | bgzip /dev/stdin -c > {output.bed} ) &> {log}
        tabix -p bed {output.bed}
        """

rule GetGWAS_LeadSnpWindows:
    """
    output bed of 1MB window surrounding lead genome-wide signficant autosomal
    snps for each gwas. Exclude blacklistregions (ie MHC). This rule only works
    for summary stats in bed format with Pvalue in column4. Exit status 1 if
    output is empty
    """
    input:
        summarystats = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        blacklistregions = "../data/MHC.hg38.bed",
        chromsizes = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        signif_loci = "gwas_summary_stats/leadSnps/{accession}.bed"
    params:
        PvalThreshold = "5e-8"
    log:
        "logs/GetGWAS_LeadSnpWindows_NonGwasCatalog/{accession}.log"
    wildcard_constraints:
        accession = '|'.join(gwas_df.index)
    resources:
        # mem_mb = much_more_mem_after_first_attempt
        mem_mb = 58000
    shell:
        """
        (python scripts/GetGWASLeadVariantWindowsFromBed.py {input.summarystats} /dev/stdout {params.PvalThreshold} | awk -F'\\t' -v OFS='\\t' '$1~/chr[0-9]+/ {{print $1, $2, $2, $1"_"$2"_N_N_{wildcards.accession}" }}' | bedtools slop -i - -g {input.chromsizes} -b 500000 | bedtools sort -i - | bedtools intersect -a - -b {input.blacklistregions} -wa -sorted -v > {output} ) &> {log}
        [[ -s {output.signif_loci} ]]
        """

rule ConcatGwasLeadSnpWindows:
    input:
        expand(
            "gwas_summary_stats/leadSnps/{accession}.bed", accession=gwas_df.index
        ),
    output:
        "gwas_summary_stats/LeadSnpWindows.bed"
    log:
        "logs/ConcatGwasLeadSnpWindows.log"
    shell:
        """
        cat {input} | bedtools sort -i - > {output}
        """

def GetScriptToOutputStandardizedStats(wildcards):
    if gwas_df.loc[wildcards.accession]['ProcessingMethod'] == 'Custom':
        #processing script A
        if wildcards.accession == "IMSGC2019":
            return "scripts/StandardizeGwasStats_MS.R"
        elif wildcards.accession == "TRANSETHNICRA":
            return "scripts/StandardizeGwasStats_RA.R"
        elif wildcards.accession == "SLEVYSE":
            return "scripts/StandardizeGwasStats_SLE.R"
        elif wildcards.accession in ["GCST004988"]:
            return ""
    elif gwas_df.loc[wildcards.accession]['ProcessingMethod'] == 'GWAS_catalog_harmonised_from_beta_se':
        return "scripts/StandardizeGwasStats_FromGwasCatalog.R"
    elif gwas_df.loc[wildcards.accession]['ProcessingMethod'] == 'GWAS_catalog_harmonised_from_OR_se':
        # ["GCST007799", "GCST007800"]
        return "scripts/StandardizeGwasStats_FromGwasCatalogOR.R"

rule GwasBedStatsAtWindows:
    """
    For coloc, gather summary stats in window centered on lead SNPs
    """
    input:
        signif_loci = "gwas_summary_stats/leadSnps/{accession}.bed",
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
    log:
        "logs/GwasBedStatsAtWindows/{accession}.log"
    output:
        stats = "gwas_summary_stats/StatsForColoc/{accession}.unstandardized.txt.gz"
    wildcard_constraints:
        accession = '|'.join(gwas_df.index)
    shell:
        """
        (tabix -h -R {input.signif_loci} {input.bed} | sort -k1,1 -k2,2n | bedtools intersect -sorted -a - -b {input.signif_loci} -wo -header | awk -F'\\t' -v OFS='\\t' 'NR==1 {{print $0, "chrom", "windowStart", "windowStop", "loci", "NOverlaps"}} NR>1 {{print}}' | gzip - > {output.stats} ) &> {log}
        """


rule GwasBedStatsToStandardizedFormat:
    """
    For coloc, gather summary stats in window centered on lead SNPs, with the
    following columns: loci (chrom_leadSNPPos_A1_A2_gwas), chrom, pos, beta,
    beta_se, A1, and A2 to ensure we can refer to the same allele set as
    molQTL. No need to worry about harmonising the order of A1 and A2 to match
    molQTL. If the GWAS summary stats don't come with A1 and A2 listed, just
    use N and N, which will skip allele checking in later steps when molQTL
    SNPs need to be matched up to the gwas SNPs. Beta effect directions also do
    not need to be harmonised for colocalization. Because some summary stats
    come in different formats, may need to use different scripts to get beta
    and beta_se (eg from Pvals)
    """
    input:
        stats = "gwas_summary_stats/StatsForColoc/{accession}.unstandardized.txt.gz",
    output:
        "gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz"
    wildcard_constraints:
        accession = '|'.join(gwas_df.index)
    params:
        StandardizeScript = GetScriptToOutputStandardizedStats
    resources:
        mem_mb = 16000
    conda:
        "../envs/r_slopes.yml"
    shell:
        """
        Rscript {params.StandardizeScript} {input.stats} {output}
        """

rule GatherGwasStandardizedStats:
    input:
        expand("gwas_summary_stats/StatsForColoc/{accession}.standardized.txt.gz", accession = gwas_df.index),
        "gwas_summary_stats/LeadSnpWindows.bed"

rule GatherAllTopSNPs:
    input:
        MolQTLs = expand("QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz", Phenotype=MyPhenotypes)
    output:
        "QTL_SNP_Enrichment/AllTopSNPs.bed.gz"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/CombinePermutationPassSNPsIntoBed.R {output} {input.MolQTLs}
        """

rule IntersectMolQTLsWithGWAS:
    input:
        bed = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz",
        tbi = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz.tbi",
        MolQTLs = "QTL_SNP_Enrichment/AllTopSNPs.bed.gz"
    output:
        "gwas_summary_stats/MolQTLIntersections/{accession}.bed.gz"
    log:
        "logs/IntersectMolQTLsWithGWAS/{accession}.log"
    shell:
        """
        (zcat {input.bed} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{print $1,$2,$3,$4}}' | bedtools intersect -wo -a - -b {input.MolQTLs} -f 1 -r -sorted | cut --complement -f 3,5,6,7 | gzip - > {output} ) &> {log}
        """

rule GatherIntersectMolQTLsWithGWAS:
    input:
        expand("gwas_summary_stats/MolQTLIntersections/{accession}.bed.gz", accession=gwas_df.index)

rule IntersectTestSNPsWithGWAS_Sample:
    input:
        molQTLs = "QTLs/QTLTools/{Phenotype}/PermutationPass.FDR_Added.txt.gz",
        gwas = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{accession}.bed.gz"
    output:
        "gwas_summary_stats/MolQTLIntersections_ControlSNPs/{accession}/{Phenotype}.txt.gz"
    wildcard_constraints:
        Phenotype="|".join(MyPhenotypes)
    params:
        n = 100000
    shell:
        """
        zcat {input.molQTLs} | awk -F' ' -v OFS='\\t' 'NR>1 {{print $2, $3, $4}}' | bedtools sort -i - | bedtools intersect -a {input.gwas} -sorted -b - -u | awk '{{print $4, "{wildcards.Phenotype}"}}' | shuf -n {params.n} | gzip - > {output}
        """

rule GatherGWAS_TestSNPs_AcrossMolQTLs:
    input:
        expand("gwas_summary_stats/MolQTLIntersections_ControlSNPs/{{accession}}/{Phenotype}.txt.gz", Phenotype=MyPhenotypes)
    output:
        "gwas_summary_stats/MolQTLIntersections_ControlSNPs/{accession}/ALL.txt.gz"
    shell:
        "cat {input} > {output}"

rule GatherGWAS_TestSNPs_AcrossMolQTLs_allTraits:
    input:
        expand("gwas_summary_stats/MolQTLIntersections_ControlSNPs/{accession}/ALL.txt.gz", accession=gwas_df.index)

rule Classify_sQTLs:
    """
    as either p-sQTLs or u-sQTLs
    """
    input:
        QTLs = "pi1/PairwisePi1Traits.P.all.txt.gz",
        IntronAnnotations = "../data/IntronAnnotationsFromYang.Updated.tsv.gz"
    output:
        OneLinePer_Cluster = "SplicingAnalysis/sQTLs_p_and_u.tsv.gz",
        OneLinePer_Intron = "SplicingAnalysis/sQTLs_p_and_u.Full.tsv.gz",
    log:
        "logs/Classify_sQTLs.log"
    conda:
        "../envs/r_2.yaml"
    shell:
        """
        Rscript scripts/Clasify_sQTLs.R {input.QTLs} {output.OneLinePer_Cluster} {output.OneLinePer_Intron} &> {log}
        """

rule MolQTL_GWAS_QQ:
    input:
        sQTLs = "SplicingAnalysis/sQTLs_p_and_u.tsv.gz",
        molQTLs = "gwas_summary_stats/MolQTLIntersections/{gwas}.bed.gz",
        OtherControlSNPs = "gwas_summary_stats/MolQTLIntersections_ControlSNPs/{gwas}/ALL.txt.gz",
        TestSNPS = "gwas_summary_stats/sorted_index_summarystat_hg38beds/{gwas}.bed.gz",
        gwas_table = "config/gwas_table.tsv"
    output:
        pdf = "gwas_summary_stats/qqplots/{gwas}/QQ.pdf",
        data = "gwas_summary_stats/qqplots/{gwas}/QQ.dat.tsv.gz"
    log:
        "logs/MolQTL_GWAS_QQ/{gwas}.log"
    conda:
        "../envs/r_scattermore.yml"
    resources:
        mem_mb = 48000
    shell:
        """
        Rscript scripts/MolQTL_GWAS_QQ.R {output.pdf} {output.data} {input.sQTLs} {input.molQTLs} {input.OtherControlSNPs} {input.TestSNPS} {wildcards.gwas} &> {log}
        """

rule copy_GWAS_qq:
    input:
        "gwas_summary_stats/qqplots/{gwas}/QQ.pdf",
    output:
        "gwas_summary_stats/All_qq/{gwas}.pdf",
    shell:
        """
        cp {input} {output}
        """

rule Gather_GWAS_QQ:
    input:
        expand("gwas_summary_stats/All_qq/{gwas}.pdf", gwas=gwas_df.index)

# rule GetGWASSummaryStatsAtLeadSNPWindows:
#     """
#     a bedtools command to get the summary stats for snps over the lead snp
#     windows, since eventually we will want to read into R to run hyprcoloc, and
#     it will be unweildy and memory intensive to have to read in huge genomewide
#     summary stats and filter for the snps in R
#     """
#     input:
#         signif_loci = "gwas_summary_stats/leadSnps/{accession}.bed",
#         summarystats = "gwas_summary_stats/full_data/{accession}.tsv.gz"
#     output:
#         signif_loci_summarystats = "gwas_summary_stats/leadSnpWindowStats/{accession}.tsv.gz"
#     log:
#         "logs/GetGWASSummaryStatsAtLeadSNPWindows/{accession}.log"
#     shell:
#         """
#         (cat <(zcat {input.summarystats} | head -1 | cut -f1-12,22-23 | awk -F'\\t' -v OFS='\\t' '{{ print $0, "lead_snp" }}' ) <(zcat {input.summarystats} | cut -f1-12,22-23 | awk -F'\\t' -v OFS='\\t' '$1!="NA" && NR>1 {{print "chr"$3, $4, $4+1, $0}}' |  bedtools intersect -a - -wo -b {input.signif_loci} | cut -f4-17,21) | gzip - > {output}) &> {log}
#         """

