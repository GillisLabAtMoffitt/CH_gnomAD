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

print("Number of variants")
nrow(gnomad)

print("Number of variants selected padj BH ≤ 0.1")
gnomad <- gnomad %>% 
  filter(adjusted_pval <= 0.1)
nrow(gnomad)

print("Number of variants selected Hedges >= small")
gnomad <- gnomad %>% 
  filter(!str_detect(effect_size_cat, "negligible|small"))
nrow(gnomad)

print("Number of homozygous")
gnomad <- gnomad %>% 
  mutate("<30_hom" = as.numeric(str_match(INFO, "age_hist_hom_n_smaller=(.*?);")[,2])) %>%
  mutate(age_hist_hom_bin_freq = str_match(INFO, "age_hist_hom_bin_freq=(.*?);")[,2]) %>%
  separate(col = age_hist_hom_bin_freq,
           into = c("30-35_hom","35-40_hom","40-45_hom","45-50_hom",
                    "50-55_hom","55-60_hom","60-65_hom","65-70_hom",
                    "70-75_hom","75-80_hom"), 
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>% 
  # "Count of age values falling above highest histogram bin edge for homozygous individuals">
  mutate(">80_hom" = as.numeric(str_match(INFO, "age_hist_hom_n_larger=(.*?);")[,2])) %>% 
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>% 
  mutate(nbr_hom_individuals = rowSums(select(.,`<30_hom`:`>80_hom`), na.rm = TRUE))

table(gnomad$nbr_hom_individuals)
str(gnomad$nbr_hom_individuals)

print("as_kable")

tbl <- gnomad %>% 
  mutate(nbr_hom_individuals_cat = as.character(nbr_hom_individuals)) %>% 
  select(nbr_hom_individuals, nbr_hom_individuals_cat) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  modify_header(
    label = '**Number of homozygous individuals per variant in selected potential CH variants**'
  ) 

tbl %>% as_kable()




