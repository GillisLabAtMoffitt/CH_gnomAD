## Import library
library(tidyverse)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/age_comparison_pval/padjust_age_selected_variants/allele_balance",
                        pattern = "*.vcf.gz",
                        recursive = FALSE,
                        full.names = TRUE)


###################################################################### II ### Calculate FDR
gnomad <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " ")) %>%
  
  select(X.CHROM, everything()) %>%
  mutate(adjusted_BA_pval = p.adjust(allele_balance_pval, method = "BH", n = length(allele_balance_pval)))

gnomad_list <- split(gnomad, gnomad$X.CHROM)

lapply(seq_along(gnomad_list), function(i){
  write_delim(gnomad_list[[i]], paste0("chr", gnomad_list[[i]][1,1], "_BA_padjusted.vcf.gz"))
})
