## Import library
library(tidyverse)





# CREATE FUNCTION
select_join <- function(gnomad, cosmic){
  selected_variants <- inner_join(gnomad, 
                                cosmic, 
                                by = c("IDs", "X.CHROM", "POS", "REF", "ALT"
                                ))
  
  # Extract allele info
  CH_variants <- selected_variants %>%
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
    mutate(controls_popmax = str_match(INFO, "controls_popmax=(.*?);")[,2])

  write.csv(head(CH_variants, 3), "test.csv")
  
  # filename <- str_match(filename,
  #                       "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/test/(.*?)\\.vcf\\.gz")[,2]
  # print(filename)
  write_delim(head(CH_variants, 3), paste0(i, "CH_variants.vcf.gz"))
  
}

# CH_variants <- select_join(gnomad = gnomad, cosmic = cosmic)


###################################################################### II ### Select gnomAD variant in cosmic, extract var
chr_vec <- paste0("chr", c("X", "Y", seq(1, 22, 1)), "_")
n <-  1
for (i in chr_vec){
  
  ### I ### List and load files
  print(paste("TOUR =", n))
  print(paste("CHR =", i))
  file_list_gnomad <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/age_selected_variants",
                                 pattern = i,
                                 recursive=TRUE,
                                 full.names = TRUE)
  
  # print(file_list_gnomad)
  gnomad <- do.call("rbind",lapply(Sys.glob(file_list_gnomad), read.delim,
                                   header = TRUE, sep = " ")) %>% 
    mutate(X.CHROM = as.character(X.CHROM))
  # print(head(gnomad, 2))
  # print(tail(gnomad, 1))
  file_list_cosmic <- list.files(path = "/share/lab_gillis/Christelle/cosmic_raw_data/splitted_data/clean_cosmic_IDs",
                                 pattern = i,
                                 recursive=TRUE,
                                 full.names = TRUE)
  # print(file_list_cosmic)

  # filename <- file_list_cosmic[[i]]
  cosmic <- read.delim(file_list_cosmic,
                       header = TRUE, sep = " ") %>% 
    mutate(X.CHROM = as.character(chr))
  print(head(cosmic, 2))
  
  ### II ### CALL FUNCTION
  print(head(select_join(gnomad = gnomad, cosmic = cosmic), 2))
  n = n + 1

}


