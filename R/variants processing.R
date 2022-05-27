library(tidyverse)
library(gtsummary)

###################################################################### I ### Load data

path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
selected_gnomad <-
  read.delim(paste0(path, "/gnomAD_raw_data/age variant selected not complete/chr19_age_selected_gnomad.vcf.gz"), 
             header = TRUE, sep = " "
             ) %>% 
  select(c(IDs:nbr_individuals, mann_w))

final_cosmic <-
  read.delim(paste0(path, "/cosmic_raw_data/splitted_data/clean_cosmic_IDs/chr19_cosmic.vcf.gz"), 
             header = TRUE, sep = " "
             )


###################################################################### II ### Binding
selected_cosmic <- inner_join(selected_gnomad, 
                              final_cosmic, 
                              by = c("IDs", "X.CHROM" = "chr", "POS", "REF", "ALT"
                                     ))

write_rds(selected_cosmic, "selected_cosmic.rds")

###################################################################### III ### Extract variant data
CH_variants <- selected_cosmic %>%
  # Generate allele count variables from the INFO var
  mutate(alt_allele_count = str_match(INFO, "AC=(.*?);")[,2]) %>%
  mutate(alt_allele_count_female = str_match(INFO, "AC_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_male = str_match(INFO, "AC_male=(.*?);")[,2]) %>%
  
  mutate(alt_allele_count_afr = str_match(INFO, "AC_afr=(.*?);")[,2]) %>%
  mutate(alt_allele_count_afr_female = str_match(INFO, "AC_afr_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_afr_male = str_match(INFO, "AC_afr_male=(.*?);")[,2]) %>%
  
  mutate(alt_allele_count_amr = str_match(INFO, "AC_amr=(.*?);")[,2]) %>%
  mutate(alt_allele_count_amr_female = str_match(INFO, "AC_amr_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_amr_male = str_match(INFO, "AC_amr_male=(.*?);")[,2]) %>%
  
  # Generate allele freq variables from the INFO var
  mutate(freq_allele_count = str_match(INFO, "AF=(.*?);")[,2]) %>%
  mutate(freq_allele_count_female = str_match(INFO, "AF_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_male = str_match(INFO, "AF_male=(.*?);")[,2]) %>%
  
  mutate(freq_allele_count_afr = str_match(INFO, "AF_afr=(.*?);")[,2]) %>%
  mutate(freq_allele_count_afr_female = str_match(INFO, "AF_afr_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_afr_male = str_match(INFO, "AF_afr_male=(.*?);")[,2]) %>%
  
  mutate(freq_allele_count_amr = str_match(INFO, "AF_amr=(.*?);")[,2]) %>%
  mutate(freq_allele_count_amr_female = str_match(INFO, "AF_amr_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_amr_male = str_match(INFO, "AF_amr_male=(.*?);")[,2]) %>%
  
  
  
  ## NON CANCER
  # Generate allele count variables from the INFO var
  mutate(nc_alt_allele_count = str_match(INFO, "non_cancer_AC=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_female = str_match(INFO, "non_cancer_AC_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_male = str_match(INFO, "non_cancer_AC_male=(.*?);")[,2]) %>%
  
  mutate(nc_alt_allele_count_afr = str_match(INFO, "non_cancer_AC_afr=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_afr_female = str_match(INFO, "non_cancer_AC_afr_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_afr_male = str_match(INFO, "non_cancer_AC_afr_male=(.*?);")[,2]) %>%
  
  mutate(nc_alt_allele_count_amr = str_match(INFO, "non_cancer_AC_amr=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_amr_female = str_match(INFO, "non_cancer_AC_amr_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_amr_male = str_match(INFO, "non_cancer_AC_amr_male=(.*?);")[,2]) %>%
  
  # Generate allele count variables from the INFO var
  mutate(nc_freq_allele_count = str_match(INFO, "non_cancer_AF=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_female = str_match(INFO, "non_cancer_AF_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_male=(.*?);")[,2]) %>%
  
  mutate(nc_freq_allele_count_afr = str_match(INFO, "non_cancer_AF_afr=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_afr_female = str_match(INFO, "non_cancer_AF_afr_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_afr_male = str_match(INFO, "non_cancer_AF_afr_male=(.*?);")[,2]) %>%
  
  mutate(nc_freq_allele_count_amr = str_match(INFO, "non_cancer_AF_amr=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_amr_female = str_match(INFO, "non_cancer_AF_amr_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_amr_male = str_match(INFO, "non_cancer_AF_amr_male=(.*?);")[,2]) %>%
  
  # "Total number of alleles in samples in the non_cancer subset, before removing low-confidence genotypes">
  mutate(nc_freq_allele_count_raw = str_match(INFO, "non_cancer_AF_raw=(.*?);")[,2]) %>% 
  
  # Depth of informative coverage for each sample
  mutate(depth_allele = str_match(INFO, "DP=(.*?);")[,2]) %>%
  
  # allele type
  mutate(allele_type = str_match(INFO, "allele_type=(.*?);")[,2]) %>%
  # variant type
  mutate(variant_type = str_match(INFO, "variant_type=(.*?);")[,2]) %>%
  # Number=A,Type=String,Description="Population with maximum AF in the controls subset">
  mutate(controls_popmax = str_match(INFO, "controls_popmax=(.*?);")[,2]) %>% 
  
  # Extract Consequences
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2])
  
CH_variants <- CH_variants %>%
  separate(col = ens_vep, paste("ens_vep", 1:25, sep=""),
           sep = ",", remove = T, extra = "warn", fill = "right") %>%
  keep(~!all(is.na(.))) %>%
  pivot_longer(cols = starts_with("ens_vep"), names_to = NULL, values_to = "ens_vep") %>%
  drop_na(ens_vep)

ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", 
                       "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", 
                       "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM", 
                       "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE",
                       "HGNC_ID", "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", 
                       "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET", "GMAF",
                       "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF", 
                       "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF",
                       "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF", "CLIN_SIG",
                       "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", 
                       "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter", "LoF_flags", "LoF_info")

CH_variants <- CH_variants %>%
  separate(col = ens_vep, into = ens_vep_var_names,
           sep = "\\|", remove = F, extra = "warn", fill = "right")

# write_rds(CH_variants, "CH_variants.rds")
# write_csv(CH_variants, "CH_variants.csv")

###################################################################### IV ### Analysis by chromosome
# Tables for
# Primary.site Primary.histology Mutation.Description FATHMM.prediction Mutation.somatic.status
# do distinct(IDs)
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(Primary.site, Primary.histology, Mutation.Description,
         FATHMM.prediction, Mutation.somatic.status) %>% 
  tbl_summary() %>% 
  bold_labels()



# Plot Allele count 
# do distinct(IDs)

# Look at final variant product/consequences
# Table for
# Consequence BIOTYPE Existing_variation VARIANT_CLASS CANONICAL
# "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "
# #   SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", "ExAC_FIN_MAF", "
# #   ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF"
# CLIN_SIG, SOMATIC -> Those 2 might be empty
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(Consequence, BIOTYPE, Existing_variation, VARIANT_CLASS,
         CANONICAL, CLIN_SIG, SOMATIC) %>% 
  tbl_summary() %>% 
  bold_labels()

# Need to summarize the mean first then do table...
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SAS_MAF, AA_MAF, EA_MAF, ExAC_MAF, ExAC_Adj_MAF, ExAC_AFR_MAF, 
         ExAC_AMR_MAF, ExAC_EAS_MAF, ExAC_FIN_MAF, 
         ExAC_NFE_MAF, ExAC_OTH_MAF, ExAC_SAS_MAF) %>% 
  tbl_summary() %>% 
  bold_labels()

