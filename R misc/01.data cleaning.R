# Import Library
library(tidyverse)
library(e1071)
library(purrr)
# library(tictoc)

################################################################################# I ### Load data
# path <- fs::path("", "Volumes", "Gillis_Research", "Christelle Colin-Leitzinger", "gnomAD")

# tic("total")
# tic("loading")
# gnomad <-
  # read.delim(paste0(path, "/data/thousands_first_rows.txt"))
# gnomad <- 
  # read.delim(paste0(path, "/data/bigger_example_variant_data.txt"))

path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
gnomad <-
  read.delim(paste0(path, "/gnomAD_raw_data/dnmt3a_data.vcf.gz"), 
             header = FALSE, 
             col.names = c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
# toc()


################################################################################# II ### Data cleaning
# tic("recoding")
gnomad_decoded <- gnomad %>% 
  mutate(X.CHROM = str_remove(X.CHROM, "chr")) %>% 
  unite(IDs, c(X.CHROM, POS, REF, ALT), sep = "-") %>% 
  # mutate(IDs = factor(row_number())) %>% 
  select(IDs, everything()) %>% 
  
  # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals; 
  # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00">
  
  # age_hist_het_bin_freq,Number=A,Type=String,Description="Histogram of ages of heterozygous individuals; 
  # bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; 
  # total number of individuals of any genotype bin: 2547|3423|4546|8487|10355|12693|11933|10534|8882|5991|4136|1935">
  mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
  
  # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
  mutate("<30" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
  # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
  mutate(">80" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>% 
  
  add_row("IDs"="All individuals", "ID"="All individuals",
          "age_hist_het_bin_freq"= "3423|4546|8487|10355|12693|11933|10534|8882|5991|4136",
          "<30" = 2547, ">80" = 1935) %>%
  
  separate(col = age_hist_het_bin_freq,
           into = c("30-35","35-40","40-45","45-50","50-55","55-60","60-65","65-70","70-75","75-80"), 
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>% 
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>% 
  
  # Filter variant found in more than 10 individuals
  mutate(nbr_individuals = rowSums(select(.,`30-35`:`>80`), na.rm = TRUE)) %>% 
  filter(nbr_individuals > 10) %>% 
  
  pivot_longer(cols = c(grep("[[:digit:]]", colnames(.))), 
               names_to = "age_bin", values_to = "sample_count") %>% 
  mutate(age_bin = factor(age_bin, 
                          levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
                                     "55-60","60-65","65-70","70-75","75-80", ">80"))) %>%
  mutate(center = case_when(age_bin == "<30" ~ 22.5,
                            age_bin == "30-35" ~ 32.5,
                            age_bin == "35-40" ~ 47.5,
                            age_bin == "40-45" ~ 42.5,
                            age_bin == "45-50" ~ 57.5,
                            age_bin == "50-55" ~ 52.5,
                            age_bin == "55-60" ~ 67.5,
                            age_bin == "60-65" ~ 62.5,
                            age_bin == "65-70" ~ 77.5,
                            age_bin == "70-75" ~ 72.5,
                            age_bin == "75-80" ~ 87.5,
                            TRUE ~ 90
  )) %>%
  mutate(left = case_when(age_bin == "<30" ~ 15,
                          age_bin == "30-35" ~ 30,
                          age_bin == "35-40" ~ 35,
                          age_bin == "40-45" ~ 40,
                          age_bin == "45-50" ~ 45,
                          age_bin == "50-55" ~ 50,
                          age_bin == "55-60" ~ 55,
                          age_bin == "60-65" ~ 60,
                          age_bin == "65-70" ~ 65,
                          age_bin == "70-75" ~ 70,
                          age_bin == "75-80" ~ 75,
                          TRUE ~ 80
  )) %>%
  mutate(right = case_when(age_bin == "<30" ~ 30,
                           age_bin == "30-35" ~ 35,
                           age_bin == "35-40" ~ 40,
                           age_bin == "40-45" ~ 45,
                           age_bin == "45-50" ~ 50,
                           age_bin == "50-55" ~ 55,
                           age_bin == "55-60" ~ 60,
                           age_bin == "60-65" ~ 65,
                           age_bin == "65-70" ~ 70,
                           age_bin == "70-75" ~ 75,
                           age_bin == "75-80" ~ 80,
                           TRUE ~ 100
  )) %>%
  mutate(sample_frequency = sample_count / nbr_individuals) %>% 
  mutate(sample_density = sample_count / (right - left)) %>% 
  mutate(sample_density0 = (sample_count / (right - left)) / nbr_individuals) %>% 
  mutate(sample_density1 = sample_frequency / (right - left))


################################################################################# III ### Goodness of Fit / Skewness : Select variant NOT consistent with the reference distribution 
# https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
goodness_fit1 <- data.frame(matrix(nrow=1, ncol=0)) 
goodness_fit2 <- data.frame(matrix(nrow=1, ncol=0)) 
skew <- data.frame(matrix(nrow=1, ncol=0)) 

for(i in unique(gnomad_decoded$IDs)) {
  
  All <- gnomad_decoded %>% filter(IDs == "All individuals") %>% select(sample_count)
  
  a <- bind_cols(gnomad_decoded %>% filter(IDs == i) %>% select(sample_count) , All) %>% 
    `colnames<-`(c(i, "All"))
  result <- chisq.test(a[,i], p=a$All, rescale.p=TRUE)$p.value
  goodness_fit1 <- rbind(goodness_fit1, result)
  
  result <- ks.test(a[,i], a$All)$p.value
  goodness_fit2 <- rbind(goodness_fit2, result)
  
  result <- skewness(a[[i]])
  skew <- rbind(skew, result)
  
}
df <- bind_cols(goodness_fit1, goodness_fit2, skew) %>% `colnames<-`(c("chisq", "kolmogorov", "skewness"))
df$IDs <- c(unique(gnomad_decoded$IDs))
# Bind the distribution characteristic values
gnomad_decoded1 <- full_join(gnomad_decoded, df, by = "IDs")# %>%
  # filter the significatant different
  # filter((chisq <= 0.05 | kolmogorov <= 0.05) & skewness > 0)


################################################################################# t-test option

data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_frequency, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)
data_plot %>% 
  ggplot(aes(x = center, y= sample_frequency, color= IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess", span = 0.6)

test_mean <- data.frame(matrix(nrow=1, ncol=0)) 

for(i in unique(gnomad_decoded$IDs)) {
  
  All <- gnomad_decoded %>% filter(IDs == "All individuals") %>% mutate(IDs = "Reference")
  data_plot <- gnomad_decoded %>% 
    filter(IDs == i) %>% 
    bind_rows(., All) %>% 
    select("IDs", "age_bin", sample_frequency, center)
  
  p <- data_plot %>% 
    ggplot(aes(x = center, y= sample_frequency, color= IDs))+
    geom_smooth(alpha= 0.5, se = FALSE, method = "loess", span = 0.6)
  p_ <- layer_data(p, 1)
  p_
  
  # a <- bind_cols(gnomad_decoded %>% filter(IDs == i) %>% select(sample_count) , All) %>% 
  #   `colnames<-`(c(i, "All"))
  # result <- chisq.test(a[,i], p=a$All, rescale.p=TRUE)$p.value
  # goodness_fit1 <- rbind(goodness_fit1, result)
  # 
  # result <- ks.test(a[,i], a$All)$p.value
  # goodness_fit2 <- rbind(goodness_fit2, result)
  
  result <- t.test(p_ %>% filter(group == 1) %>% select(y), p_ %>% filter(group == 2) %>% select(y))$p.value
  test_mean <- rbind(test_mean, result)
  print(test_mean)
}
test_mean <- test_mean %>% `colnames<-`(c("t_test")) %>% 
  mutate(IDs = c(unique(gnomad_decoded$IDs)))

# Bind the distribution characteristic values
gnomad_decoded1 <- full_join(gnomad_decoded1, test_mean, by = "IDs") %>%
  # filter the significatant different
  filter(t_test > 0)




################################################################################# VI ### Filtering Consequence

gnomad_decoded2 <- gnomad_decoded1 %>% 
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2]) %>% 
  separate(col = ens_vep, paste("ens_vep", 1:25, sep=""),
           sep = ",", remove = T, extra = "warn", fill = "right") %>% 
  keep(~!all(is.na(.))) %>% 
  pivot_longer(cols = starts_with("ens_vep"), names_to = NULL, values_to = "ens_vep") %>% 
  drop_na(ens_vep)
# Consequence annotations from Ensembl VEP
# Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|
# cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|
# STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|
# TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|
# SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|
# ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|
# HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info
ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "
  cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "
  STRAND", "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "
  TREMBL", "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET", "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "
  SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", "ExAC_FIN_MAF", "
  ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "
  HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter", "LoF_flags", "LoF_info")

gnomad_decoded3 <- gnomad_decoded2 %>% 
  separate(col = ens_vep, into = ens_vep_var_names,
           # paste("ens_vep", 1:100, sep=""),
           sep = "\\|", remove = F, extra = "warn", fill = "right")

table(gnomad_decoded3$Consequence)
table(gnomad_decoded3$SYMBOL)

gnomad_decoded4 <- gnomad_decoded3 %>% 
  filter(!str_detect(Consequence, "3_prime_UTR_variant|5_prime_UTR_variant"), SYMBOL == "DNMT3A")

# toc() # 2 sec elapsed
# toc() # 4.437 sec elapsed 21 hrs for all rows


################################################################################# VI ### Plotting

data_plot <- gnomad_decoded4 %>% 
  distinct(IDs, age_bin, .keep_all = TRUE) %>% 
  select("IDs", "nbr_individuals", "age_bin", "sample_count", sample_frequency, skewness) %>% 
  mutate(age_bin = factor(age_bin, 
                          levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
                                     "55-60","60-65","65-70","70-75","75-80", ">80")))

# data_plot <- data_plot[1:60,]

# data_plot <- data_plot %>% 
#   group_by(ID) %>% 
#   pivot_longer(cols = c("<30", "30-35","35-40","40-45","45-50","50-55","55-60",
#                         "60-65","65-70","70-75","75-80", ">80"), 
#                names_to = "age_bin", values_to = "count_in_bin") %>% 
#   ungroup() %>% 
#   mutate(age_bin = factor(age_bin, 
#                           levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
#                                      "55-60","60-65","65-70","70-75","75-80", ">80"))) %>% 
#   mutate(left = case_when(age_bin == "<30" ~ 0,
#                           age_bin == "30-35" ~ 30,
#                           age_bin == "35-40" ~ 35,
#                           age_bin == "40-45" ~ 40,
#                           age_bin == "45-50" ~ 45,
#                           age_bin == "50-55" ~ 50,
#                           age_bin == "55-60" ~ 55,
#                           age_bin == "60-65" ~ 60,
#                           age_bin == "65-70" ~ 65,
#                           age_bin == "70-75" ~ 70,
#                           age_bin == "75-80" ~ 75,
#                           TRUE ~ 80
#   )) %>% 
#   mutate(right = case_when(age_bin == "<30" ~ 30,
#                            age_bin == "30-35" ~ 35,
#                            age_bin == "35-40" ~ 40,
#                            age_bin == "40-45" ~ 45,
#                            age_bin == "45-50" ~ 50,
#                            age_bin == "50-55" ~ 55,
#                            age_bin == "55-60" ~ 60,
#                            age_bin == "60-65" ~ 65,
#                            age_bin == "65-70" ~ 70,
#                            age_bin == "70-75" ~ 75,
#                            age_bin == "75-80" ~ 80,
#                            TRUE ~ 90
#   )) %>%
#   mutate(center = case_when(age_bin == "<30" ~ 15,
#                             age_bin == "30-35" ~ 32.5,
#                             age_bin == "35-40" ~ 47.5,
#                             age_bin == "40-45" ~ 42.5,
#                             age_bin == "45-50" ~ 57.5,
#                             age_bin == "50-55" ~ 52.5,
#                             age_bin == "55-60" ~ 67.5,
#                             age_bin == "60-65" ~ 62.5,
#                             age_bin == "65-70" ~ 77.5,
#                             age_bin == "70-75" ~ 72.5,
#                             age_bin == "75-80" ~ 87.5,
#                             TRUE ~ 90
#   ))
# filter(count_in_bin != 0) #%>% 
# select(ID, age_bin, count_in_bin, everything())

ref <- tibble(age_bin = c("<30", "30-35","35-40","40-45","45-50","50-55","55-60",
                          "60-65","65-70","70-75","75-80", ">80"),
              y = c(2547, 3423, 4546, 8487, 10355, 12693, 11933, 10534, 8882, 5991, 4136, 1935)) %>% 
  mutate(age_bin = factor(age_bin, 
                          levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
                                     "55-60","60-65","65-70","70-75","75-80", ">80"))) %>% 
  mutate(sum = sum(y)) %>% 
  mutate(density = y / sum)
sum <- unique(ref$sum)

data_plot <-  full_join(data_plot, ref, by = "age_bin")
ggplot(data = ref, aes(x =age_bin, y = density))+
  geom_col(fill = "yellow", alpha = 0.5)+
  theme_minimal()



ggplot() + 
  geom_col(data =ref, aes(x =age_bin, y = density), position = "identity", fill = "yellow", alpha = 0.5)+
  geom_col(data = data_plot ,aes(x= age_bin, y= sample_frequency), position = "identity", alpha = 0.5)+
  theme_minimal()


ggplot() + 
  geom_col(data = data_plot ,aes(x= age_bin, y= sample_frequency), position = "identity", alpha = 0.6)+
  theme_minimal()+
  facet_wrap(.~ IDs, scales = "free_y")+
  geom_col(data =ref, aes(x =age_bin, y = density), position = "identity", fill = "yellow", alpha = 0.5)+ 
  geom_text(
    data    = data_plot,
    mapping = aes(x = -Inf, y = -Inf, label = skewness),
    hjust   = 0,
    vjust   = 0
  )
  


# ggplot() + 
#   geom_col(data =ref, aes(x =age_bin, y = density), position = "identity", fill = "yellow", alpha = 0.5)+
#   geom_col(data = data_plot ,aes(x= age_bin, y= sample_count), position = "identity")+
#   
#   scale_y_continuous("Sample Counts", sec.axis = sec_axis( ~ .  * sum, name = "Temperature")) +
#   
#   theme_minimal()+
#   facet_wrap(.~ IDs, scales = "free_y")

data_plot %>% 
  ggplot(aes(x= age_bin, y= sample_count))+
  geom_col()+
  theme_minimal()+
  facet_wrap(.~ IDs, scales = "free_y")




gnomad_decoded2 <- gnomad_decoded1 %>% 
  pivot_wider(names_from = "age_bin", values_from = "sample_count")

################################################################################# IV ### Finishing extracting variables on the subsets of variant selected
gnomad_decoded3 <- gnomad_decoded2 %>% 
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
  mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_raw=(.*?);")[,2])
  
  








