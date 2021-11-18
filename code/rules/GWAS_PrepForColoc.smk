rule GetGWAS_LeadSnpWindows:
    """
    output bed of 1MB window surrounding lead genome-wide signficant autosomal
    snps for each gwas. Exclude blacklistregions (ie MHC)
    """
    input:
        summarystats = "gwas_summary_stats/full_data/{accession}.tsv.gz",
        blacklistregions = "../data/MHC.hg38.bed",
        chromsizes = "ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        signif_loci = "gwas_summary_stats/leadSnps/{accession}.bed"
    params:
        PvalThreshold = "5e-8"
    log:
        "logs/GetGWAS_LeadSnpWindows/{accession}.log"
    resources:
        mem_mb = 48000
    shell:
        """
        (python scripts/GetGWASLeadVariantWindows.py {input.summarystats} /dev/stdout {params.PvalThreshold} | awk -F'\\t' -v OFS='\\t' 'NR>1 && $3~/[0-9]+/ {{split($1, a, "_"); print "chr"$3, a[2], a[2], $1"_"$2"_{wildcards.accession}" }}' | bedtools slop -i - -g {input.chromsizes} -b 500000 | bedtools sort -i - | bedtools intersect -a - -b {input.blacklistregions} -wa -sorted -v > {output} ) &> {log}
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

rule GetGWASSummaryStatsAtLeadSNPWindows:
    """
    # TODO: a bedtools command to get the summary stats for snps over the lead
    # snp windows, since eventually we will want to read into R to run
    # hyprcoloc, and it will be unweildy and memory intensive to have to read
    # in huge genomewide summary stats and filter for the snps in R
    """
    input:
        signif_loci = "gwas_summary_stats/leadSnps/{accession}.bed"
    output:
        signif_loci = "gwas_summary_stats/leadSnpWindowStats/{accession}.bed"
    log:
        "logs/GetGWASSummaryStatsAtLeadSNPWindows.log"
    shell:
        """
        """
