# this script will help me make a tsv file (one row per gwas) for storing accession numbers and other relevant info for downloading and processing gwas
# the final tsv file will need at minimum the following columns:
## gwas - unqiue gwas label
## trait - trait name
## either `FTPPath` (for downloading from gwas catalog), or `SummaryStatsLocalFilepath` for summary stats already downloaded
## ProcessingMethod - possible values include "GWAS_catalog_harmonised_from_beta_se", "GWAS_catalog_harmonised_from_OR_p", "Custom"
## Continuous - Either "True" (for continuous traits) or "False" (for binary traits)

# ProcessingMethod specifies scripts/methods to generate "gwas_summary_stats/leadSnpWindowStats/{gwas}.tsv.gz" with summary stats for each 1MB leadSNP window, with at minimum columns labelled "snp_chr", "snp_pos", "ref", "alt", "p", "beta", "se", "lead_snp" (lead.snp.chr_pos_ref_alt_snpID_{gwas})
## Note that gwas catalog hm polarizes effect towards reference allele, which is consistent with molQTL summary stats. Any processing method that does not use hm_beta, must be sure to appropriately polarize betas in a manner that is consistent relative to molQTL betas

# run from code dir


library(tidyverse)

# read the csvs from gwas catalog
csvs_from_gwas_catalog <- list.files("../data/GWAS_catalog_summary_stats_sources", pattern = "*.csv", full.names = T) %>%
  setNames(str_replace(., "../data/GWAS_catalog_summary_stats_sources/(.+?).csv", "\\1"))

dat <- lapply(csvs_from_gwas_catalog, read_csv) %>%
  bind_rows(.id="csv")

# AFTER downloading the summary stats from gwas catalog (invoke snakemake to do this), peak inside the stats to see what columns are present
gwas_catalog_summary_stats_fns <- list.files("gwas_summary_stats/full_data", pattern = "*.tsv.gz", full.names = T) %>%
  setNames(str_replace(., "gwas_summary_stats/full_data/(.+?).tsv.gz", "\\1"))

gwas_catalog_summary_stats_preview <- lapply(gwas_catalog_summary_stats_fns, read.table, nrows=10, sep='\t', header=T) %>%
  bind_rows(.id="accession")

colnames(gwas_catalog_summary_stats_preview)

# from the preview of the summary stats, determine values for required columns for tsv to write out
gwas_catalog_rows <- gwas_catalog_summary_stats_preview %>%
  filter(!is.na(hm_variant_id)) %>%
  distinct(accession, .keep_all=T) %>%
  select(1:9, p_value, standard_error) %>%
  mutate(ProcessingMethod = case_when(
    !is.na(standard_error) & !is.na(hm_beta) ~ "GWAS_catalog_harmonised_from_beta_se",
    !is.na(p_value) & !is.na(hm_odds_ratio) ~ "GWAS_catalog_harmonised_from_OR_p",
    TRUE ~ "Custom"
  )) %>%
  mutate(Continuous = case_when(
    !is.na(hm_odds_ratio) ~ F,
    TRUE ~ T
  )) %>%
  select(accession, ProcessingMethod, Continuous) %>%
  inner_join(dat, by=c("accession"="Study accession")) %>%
  rename(gwas=accession, trait=`Reported trait`, FTPPath=`FTP Path`) %>%
  mutate(SummaryStatsLocalFilepath=NA_character_) %>%
  select(gwas, trait, FTPPath, SummaryStatsLocalFilepath, ProcessingMethod, Continuous, everything(), -`Usage License`)


# add extra rows for other gwas that are not from gwas catalog, and write out final tsv
gwas_catalog_rows %>%
  add_row(
    gwas = "IMSGC2019",
    trait = "Multiple sclerosis",
    SummaryStatsLocalFilepath = "/project2/yangili1/bjf79/gwas_summary_stats/discovery_metav3.0.meta.gz",
    Continuous = F,
    ProcessingMethod = "Custom"
  ) %>%
  write_tsv("config/gwas_table.tsv")

