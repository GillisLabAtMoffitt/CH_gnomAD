## Import library
library(tidyverse)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants/allele_balance/age_distribution/variant_depth/noise_filter",
                        pattern = "*.vcf.gz",
                        recursive = FALSE,
                        full.names = TRUE)


###################################################################### II ### Calculate FDR
gnomad <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " "))
object.size(gnomad)

AB_adjusted <- gnomad %>%
  filter(contaminated_AB == "clean sequencing" &
           black_list == 0 &
           centromeres == 0 &
           segmental_duplication == 0 &
           wm_sdust == 0 &
           dp_het_median >= 10 &
           is_gencode_protein_coding_gene == "Yes") %>%
  select(IDs, allele_balance_pval) %>%
  filter(!is.na(allele_balance_pval)) %>%
  mutate(adjusted_AB_pval = p.adjust(allele_balance_pval, method = "BH", n = length(allele_balance_pval))) %>%
  select(-allele_balance_pval)
object.size(AB_adjusted)

age_adjusted <- gnomad %>%
  filter(contaminated_AB == "clean sequencing" &
           black_list == 0 &
           centromeres == 0 &
           segmental_duplication == 0 &
           wm_sdust == 0 &
           dp_het_median >= 10 &
           is_gencode_protein_coding_gene == "Yes") %>%
  select(IDs, age_mann_pval) %>%
  filter(!is.na(age_mann_pval)) %>%
  mutate(adjusted_age_pval = p.adjust(age_mann_pval, method = "BH", n = length(age_mann_pval))) %>%
  select(-age_mann_pval)
object.size(age_adjusted)

write_delim(AB_adjusted, "AB_adjusted.vcf.gz")
write_delim(age_adjusted, "age_adjusted.vcf.gz")


gnomad <- gnomad %>%
  full_join(., AB_adjusted, by= "IDs") %>%
  full_join(., age_adjusted, by= "IDs")

write_delim(gnomad, "gnomad_long_data.vcf.gz")
            
rm(file_list, AB_adjusted, age_adjusted)
gc()

gnomad_list <- split(gnomad, gnomad$X.CHROM)
object.size(gnomad_list)

lapply(seq_along(gnomad_list), function(i){
  write_delim(gnomad_list[[i]], paste0("chr", gnomad_list[[i]][1,2], "_padjusted.vcf.gz"))
})
