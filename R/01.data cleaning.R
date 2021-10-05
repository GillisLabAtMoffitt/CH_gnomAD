# Import Library
library(tidyverse)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research", "Christelle Colin-Leitzinger", "gnomAD")

gnomad <- 
  read.delim(paste0(path, "/data/example_variant_data.txt"))


################################################################################# II ### Data cleaning
gnomad_decoded <- gnomad %>% 
  # Generate allele count variables from the INFO var
  mutate(alt_allele_count = str_match(INFO, "AC=(.*?);")[,2]) %>%
  mutate(alt_allele_count_afr_female = str_match(INFO, "AC_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_male = str_match(INFO, "AC_male=(.*?);")[,2]) %>%
  
  mutate(alt_allele_count_afr = str_match(INFO, "AC_afr=(.*?);")[,2]) %>%
  mutate(alt_allele_count_afr_female = str_match(INFO, "AC_afr_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_male = str_match(INFO, "AC_afr_male=(.*?);")[,2]) %>%
  
  mutate(alt_allele_count_afr = str_match(INFO, "AC_amr=(.*?);")[,2]) %>%
  mutate(alt_allele_count_afr_female = str_match(INFO, "AC_amr_female=(.*?);")[,2]) %>%
  mutate(alt_allele_count_male = str_match(INFO, "AC_amr_male=(.*?);")[,2]) %>%
  
  # Generate allele count variables from the INFO var
  mutate(freq_allele_count = str_match(INFO, "AF=(.*?);")[,2]) %>%
  mutate(freq_allele_count_afr_female = str_match(INFO, "AF_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_male = str_match(INFO, "AF_male=(.*?);")[,2]) %>%
  
  mutate(freq_allele_count_afr = str_match(INFO, "AF_afr=(.*?);")[,2]) %>%
  mutate(freq_allele_count_afr_female = str_match(INFO, "AF_afr_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_male = str_match(INFO, "AF_afr_male=(.*?);")[,2]) %>%
  
  mutate(freq_allele_count_afr = str_match(INFO, "AF_amr=(.*?);")[,2]) %>%
  mutate(freq_allele_count_afr_female = str_match(INFO, "AF_amr_female=(.*?);")[,2]) %>%
  mutate(freq_allele_count_male = str_match(INFO, "AF_amr_male=(.*?);")[,2]) %>%
  
  # Depth of informative coverage for each sample
  mutate(freq_allele_count = str_match(INFO, "DP=(.*?);")[,2]) %>%
  # GQ
  
  
  
  
  
  # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals; 
  # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00">
  
  # age_hist_het_bin_freq,Number=A,Type=String,Description="Histogram of ages of heterozygous individuals; 
  # bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; 
  # total number of individuals of any genotype bin: 2547|3423|4546|8487|10355|12693|11933|10534|8882|5991|4136|1935">
  mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
  separate(col = age_hist_het_bin_freq,
           into = c("≥30.0", "≥35.0", "≥40.0", "≥45.0", "≥50.0", "≥55.0", "≥60.0", "≥65.0", "≥70.0", "75.0 - 80"), 
           sep = "\\|", remove = F, extra = "warn", fill = "right") %>% 
  # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
  mutate(">80" = str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2]) %>%
  # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
  mutate("<30" = str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2]) %>%
  # allele type
  mutate(allele_type = str_match(INFO, "allele_type=(.*?);")[,2]) %>%
  # variant type
  mutate(variant_type = str_match(INFO, "variant_type=(.*?);")[,2]) %>%
  # Number=A,Type=String,Description="Population with maximum AF in the controls subset">
  mutate(controls_popmax = str_match(INFO, "controls_popmax=(.*?);")[,2]) %>%
  
  
  ## NON CANCER
  # Generate allele count variables from the INFO var
  mutate(nc_alt_allele_count = str_match(INFO, "non_cancer_AC=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_afr_female = str_match(INFO, "non_cancer_AC_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_male = str_match(INFO, "non_cancer_AC_male=(.*?);")[,2]) %>%
  
  mutate(nc_alt_allele_count_afr = str_match(INFO, "non_cancer_AC_afr=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_afr_female = str_match(INFO, "non_cancer_AC_afr_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_male = str_match(INFO, "non_cancer_AC_afr_male=(.*?);")[,2]) %>%
  
  mutate(nc_alt_allele_count_afr = str_match(INFO, "non_cancer_AC_amr=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_afr_female = str_match(INFO, "non_cancer_AC_amr_female=(.*?);")[,2]) %>%
  mutate(nc_alt_allele_count_male = str_match(INFO, "non_cancer_AC_amr_male=(.*?);")[,2]) %>%
  
  # Generate allele count variables from the INFO var
  mutate(nc_freq_allele_count = str_match(INFO, "non_cancer_AF=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_afr_female = str_match(INFO, "non_cancer_AF_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_male=(.*?);")[,2]) %>%
  
  mutate(nc_freq_allele_count_afr = str_match(INFO, "non_cancer_AF_afr=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_afr_female = str_match(INFO, "non_cancer_AF_afr_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_afr_male=(.*?);")[,2]) %>%
  
  mutate(nc_freq_allele_count_afr = str_match(INFO, "non_cancer_AF_amr=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_afr_female = str_match(INFO, "non_cancer_AF_amr_female=(.*?);")[,2]) %>%
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_amr_male=(.*?);")[,2]) %>%
  
  # "Total number of alleles in samples in the non_cancer subset, before removing low-confidence genotypes">
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_raw=(.*?);")[,2]) %>%
  
  # Consequence annotations from Ensembl VEP
  # Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|
  # cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|
  # STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|
  # TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|
  # SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|
  # ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|
  # HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2]) #%>% 
  # separate(col = ens_vep, paste("ens_vep", 1:67, sep=""), 
  #          sep = "\\|", remove = F, extra = "warn", fill = "right")




str_count("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|
  cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|
  STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|
  TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|
  SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|
  ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|
  HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info", "\\|")





