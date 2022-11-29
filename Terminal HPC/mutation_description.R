## Import library
library(tidyverse)
library(gtsummary)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants",
                        pattern = "*.vcf.gz",
                        recursive=FALSE,
                        full.names = TRUE)


###################################################################### II ### Table
gnomad <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " "))

gnomad %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  select(Mutation.Description) %>% 
  tbl_summary() %>% 
  bold_labels() %>% 
  modify_header(
    label = '**Mutation Characteristics for Cosmic Patients/Cell Line**'
  ) %>% as_gt()
gt::gtsave("Cosmic mutation descriptions in all chromosomes.pdf")
