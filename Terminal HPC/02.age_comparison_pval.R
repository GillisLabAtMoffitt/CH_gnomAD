## Import library
library(tidyverse)
library(effsize)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data",
                        pattern = "*.vcf.gz",
                        recursive = FALSE,
                        full.names = TRUE)


# CREATE FUNCTION
age_mannw <- function(gnomad){
  
  gnomad <- gnomad %>%
    mutate(X.CHROM = str_remove(X.CHROM, "chr")) %>%
    unite(IDs, c(X.CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>%
    
    # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals;
    # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00">
    
    # 1. Extract count per age bin for heterozygous ind--
    # age_hist_het_bin_freq,Number=A,Type=String,Description="Histogram of ages of heterozygous individuals;
    # bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0;
    # total number of individuals of any genotype bin: 2547|3423|4546|8487|10355|12693|11933|10534|8882|5991|4136|1935">
    mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
    # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
    mutate("<30_het" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
    # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
    mutate(">80_het" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>%
    separate(col = age_hist_het_bin_freq,
             into = c("30-35_het","35-40_het","40-45_het","45-50_het",
                      "50-55_het","55-60_het","60-65_het","65-70_het",
                      "70-75_het","75-80_het"),
             sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
    mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
    
    # 2. Extract count per age bin for heterozygous ind--
    mutate(age_hist_hom_bin_freq = str_match(INFO, "age_hist_hom_bin_freq=(.*?);")[,2]) %>%
    # "Count of age values falling below lowest histogram bin edge for homozygous individuals">
    mutate("<30_hom" = as.numeric(str_match(INFO, "age_hist_hom_n_smaller=(.*?);")[,2])) %>%
    # "Count of age values falling above highest histogram bin edge for homozygous individuals">
    mutate(">80_hom" = as.numeric(str_match(INFO, "age_hist_hom_n_larger=(.*?);")[,2])) %>%
    separate(col = age_hist_hom_bin_freq,
             into = c("30-35_hom","35-40_hom","40-45_hom","45-50_hom",
                      "50-55_hom","55-60_hom","60-65_hom","65-70_hom",
                      "70-75_hom","75-80_hom"),
             sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
    mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
    
    # 3. Combine count for Het and hom individuals
    mutate("<30" = rowSums(select(.,`30-35_het`,`30-35_hom`), na.rm = TRUE),
           "30-35" = rowSums(select(.,`30-35_het`,`30-35_hom`), na.rm = TRUE),
           "35-40" = rowSums(select(.,`35-40_het`,`35-40_hom`), na.rm = TRUE),
           "40-45" = rowSums(select(.,`40-45_het`,`40-45_hom`), na.rm = TRUE),
           "45-50" = rowSums(select(.,`45-50_het`,`45-50_hom`), na.rm = TRUE),
           "50-55" = rowSums(select(.,`50-55_het`,`50-55_hom`), na.rm = TRUE),
           "55-60" = rowSums(select(.,`55-60_het`,`55-60_hom`), na.rm = TRUE),
           "60-65" = rowSums(select(.,`60-65_het`,`60-65_hom`), na.rm = TRUE),
           "65-70" = rowSums(select(.,`65-70_het`,`65-70_hom`), na.rm = TRUE),
           "70-75" = rowSums(select(.,`70-75_het`,`70-75_hom`), na.rm = TRUE),
           "75-80" = rowSums(select(.,`75-80_het`,`75-80_hom`), na.rm = TRUE),
           ">80" = rowSums(select(.,`>80_het`,`>80_hom`), na.rm = TRUE)
    ) %>%
    add_row("IDs"="All individuals", "ID"="All individuals of any genotype bin",
            "<30" = 2547, "30-35" = 3423,
            "35-40" = 4546,"40-45" = 8487,
            "45-50" = 10355,"50-55" = 12693,
            "55-60" = 11933,"60-65" = 10534,
            "65-70" = 8882,"70-75" = 5991,
            "75-80" = 4136,
            ">80" = 1935) %>%
    mutate(nbr_individuals = rowSums(select(.,`<30`:`>80`), na.rm = TRUE)) %>%
    
    select(-c(ends_with("_het"), ends_with("_hom"))) %>%
    # Filter variants found in more than 10 individuals
    filter(nbr_individuals > 10) %>%
    # Generate AN variable from the INFO var
    # mutate(AN = str_match(INFO, "AN=(.*?);")[,2]) %>%
    # mutate(AN_samples_count = (as.numeric(AN) / 2)
    # ) %>%
    select(IDs, everything())
  
  gnomad <- gnomad %>%
    pivot_longer(cols = c(grep("[[:digit:]]", colnames(.))),
                 names_to = "age_bin", values_to = "sample_count") %>%
    mutate(age_bin = factor(age_bin,
                            levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
                                       "55-60","60-65","65-70","70-75","75-80", ">80"
                            ))) %>%
    mutate(bin_age = case_when(
      # The average age at diagnosis is 8 overall (ages 0 to 19),
      # 5 years old for children (aged 0 to 14), and 17 years old
      # for adolescents (aged 15 to 19), while adults' average
      # age for cancer diagnosis is 65
      str_detect(age_bin, "<30")              ~ 15,
      str_detect(age_bin, "30-35")            ~ 32.5,
      str_detect(age_bin, "35-40")            ~ 37.5,
      str_detect(age_bin, "40-45")            ~ 42.5,
      str_detect(age_bin, "45-50")            ~ 47.5,
      str_detect(age_bin, "50-55")            ~ 52.5,
      str_detect(age_bin, "55-60")            ~ 57.5,
      str_detect(age_bin, "60-65")            ~ 62.5,
      str_detect(age_bin, "65-70")            ~ 67.5,
      str_detect(age_bin, "70-75")            ~ 72.5,
      str_detect(age_bin, "75-80")            ~ 77.5,
      TRUE ~ 90
    )) %>%
    mutate(sample_frequency = sample_count / nbr_individuals) #%>%
  # mutate(variant_frequency = sample_count / AN_samples_count)
  
  
  ## Step 1 : Limit variants with a skewed age distribution
  mann_w <-  data.frame(matrix(nrow=1, ncol=0))
  mann_pval <-  data.frame(matrix(nrow=1, ncol=0))
  variant_mean <-  data.frame(matrix(nrow=1, ncol=0))
  effect_size <-  data.frame(matrix(nrow=1, ncol=0)) 
  effect_size_cat <-  data.frame(matrix(nrow=1, ncol=0))
  set.seed(1234)
  for(i in unique(gnomad$IDs)) {
    
    All <- gnomad %>%
      filter(IDs == "All individuals"
      ) %>%
      mutate(IDs = "Reference")
    data_plot <- gnomad %>%
      filter(IDs == i
      ) %>%
      bind_rows(., All) %>%
      select("IDs", "age_bin", sample_frequency, bin_age)
    
    p <- data_plot %>%
      ggplot(aes(x = bin_age, y= sample_frequency, color= IDs))+
      # ggplot(aes(x = age_bin, y= sample_frequency, color= IDs))+
      # geom_bar(stat = "identity", aes(alpha= IDs))
      geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
                  # method.args=list(family=quasibinomial)#,
                  span = 0.6
      )+
      ylim(0, max(data_plot$sample_frequency))
    p_ <- layer_data(p, 1) %>%
      mutate(colour =  as.factor(colour)) %>%
      mutate_if(is.numeric, ~replace(., is.na(.), 0))
    
    variant <- p_ %>% filter(group == 1)
    variant <- sample(variant$x, 180, prob = variant$y, replace = TRUE)
    ref <- p_ %>% filter(group == 2)
    ref <- sample(ref$x, 180, prob = ref$y, replace = TRUE)
    
    result <- wilcox.test(variant, ref, 
                          alternative = "greater", 
                          paired=FALSE)
    mann_pval <- rbind(mann_pval, result$p.value)
    mann_w <- rbind(mann_w, result$statistic)
    var_mean <- t.test(variant)$estimate
    variant_mean <- rbind(variant_mean, var_mean)
    
    result2 <- cohen.d(variant, ref, hedges.correction=TRUE)$estimate
    effect_size <- rbind(effect_size, result2)
    result3 <- cohen.d(variant, ref,hedges.correction=TRUE)$magnitude[1] %>% 
      as_tibble() %>%
      mutate(value = case_when(
        value == "negligible"   ~ "negligible",
        value == "small"        ~ "small",
        value == "medium"       ~ "medium",
        value == "large"        ~ "large"
      ))
    effect_size_cat <- rbind(effect_size_cat, result3)
    
  }
  mann_w <- mann_w %>% 
    bind_cols(mann_pval) %>% 
    bind_cols(variant_mean) %>% 
    bind_cols(effect_size) %>% 
    bind_cols(effect_size_cat) %>% 
    `colnames<-`(c("mann_w", "mann_pval", "distribution_mean", "effect_size", "effect_size_cat")) %>%
    mutate(IDs = c(unique(gnomad$IDs)))
  
  # Bind the distribution characteristic values
  gnomad <- gnomad %>%
    full_join(., mann_w, by = "IDs") %>%
    distinct(IDs, .keep_all = TRUE) %>% 
    filter(IDs != "All individuals")
  
}


###################################################################### II ### Reading and Cleaning
for (i in seq_along(file_list)){
  filename <- file_list[[i]]

  gnomad <- read.delim(filename,
                       header = FALSE,
                       col.names = c("X.CHROM", "POS", "ID", "REF", "ALT", "INFO")
                       )
  
  gnomad <- age_mannw(gnomad = gnomad)
  
  filename <- str_match(filename, 
                        "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/(.*?)\\.vcf\\.gz")[,2]
  
  write_delim(gnomad, paste0(filename, "_gnomad_age_mannw.vcf.gz"))

}



