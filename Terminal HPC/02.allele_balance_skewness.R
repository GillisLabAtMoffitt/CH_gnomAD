## Import library
library(tidyverse)
library(effsize)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data",
                        pattern = "*.vcf.gz",
                        recursive = FALSE,
                        full.names = TRUE)


# CREATE FUNCTION
allele_balance <- function(gnomad){
  
  gnomad <- gnomad %>%
    mutate(X.CHROM = str_remove(X.CHROM, "chr")) %>%
    unite(IDs, c(X.CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>%
    # 1. Extract count per AB bin for heterozygous ind--
    # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals; 
    # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00>
    mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
    separate(col = ab_hist_alt_bin_freq,
             into = c("0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
                      "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
                      "0.40-0.45","0.45-0.50","0.50-0.55","0.55-0.60",
                      "0.60-0.65","0.65-0.70","0.70-0.75","0.75-0.80",
                      "0.80-0.85","0.85-0.90","0.90-0.95","0.95-1.00"),
             sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
    mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
    add_row("IDs"="Reference", "ID"="DNMT3A 2-25247044-C-T",
            "0.00-0.05" = 0, "0.05-0.10"  = 1, "0.10-0.15" = 22, 
            "0.15-0.20" = 53, "0.20-0.25" = 147, "0.25-0.30" = 545, 
            "0.30-0.35" = 1592, "0.35-0.40" = 3438, "0.40-0.45" = 8408,
            "0.45-0.50" = 9816, "0.50-0.55" = 12487, "0.55-0.60" = 7145,
            "0.60-0.65" = 3974, "0.65-0.70" = 1399, "0.70-0.75" = 419,
            "0.75-0.80" = 208, "0.80-0.85" = 80, "0.85-0.90" = 84,
            "0.90-0.95" = 92, "0.95-1.00" = 0) %>%
    mutate(total_allele_balance = rowSums(select(.,`0.00-0.05`:`0.95-1.00`), na.rm = TRUE)) %>%
    filter(total_allele_balance > 5) %>%
  
    select(IDs, everything())
  
  gnomad <- gnomad %>%
    pivot_longer(cols = c(grep("[[:digit:]]", colnames(.))),
                 names_to = "allele_bin", values_to = "AB_count") %>%
    mutate(allele_bin = factor(allele_bin,
                            levels = c("0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
                                       "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
                                       "0.40-0.45","0.45-0.50","0.50-0.55","0.55-0.60",
                                       "0.60-0.65","0.65-0.70","0.70-0.75","0.75-0.80",
                                       "0.80-0.85","0.85-0.90","0.90-0.95","0.95-1.00"
                            ))) %>%
    mutate(AB_bin = case_when(
      allele_bin == "0.00-0.05"            ~ 0.025,
      allele_bin == "0.05-0.10"            ~ 0.075,
      allele_bin == "0.10-0.15"            ~ 0.125,
      allele_bin == "0.15-0.20"            ~ 0.175,
      allele_bin == "0.20-0.25"            ~ 0.225,
      allele_bin == "0.25-0.30"            ~ 0.275,
      allele_bin == "0.30-0.35"            ~ 0.325,
      allele_bin == "0.35-0.40"            ~ 0.375,
      allele_bin == "0.40-0.45"            ~ 0.425,
      allele_bin == "0.45-0.50"            ~ 0.475,
      allele_bin == "0.50-0.55"            ~ 0.525,
      allele_bin == "0.55-0.60"            ~ 0.575,
      allele_bin == "0.60-0.65"            ~ 0.625,
      allele_bin == "0.65-0.70"            ~ 0.675,
      allele_bin == "0.70-0.75"            ~ 0.725,
      allele_bin == "0.75-0.80"            ~ 0.775,
      allele_bin == "0.80-0.85"            ~ 0.825,
      allele_bin == "0.85-0.90"            ~ 0.875,
      allele_bin == "0.90-0.95"            ~ 0.925,
      allele_bin == "0.95-1.00"            ~ 0.975
    )) %>%
    mutate(AB_frequency = AB_count / total_allele_balance)
  
  
  ## Step 1 : Limit variants with a skewed age distribution
  allele_balance_pval <-  data.frame(matrix(nrow=length(unique(gnomad$IDs)), ncol=1,
                                  dimnames = list(c(unique(gnomad$IDs)),
                                                  c("value")) ))
  variant_mean <-  data.frame(matrix(nrow=length(unique(gnomad$IDs)), ncol=1,
                                     dimnames = list(c(unique(gnomad$IDs)),
                                                     c("value")) ))
  effect_size <-  data.frame(matrix(nrow=length(unique(gnomad$IDs)), ncol=1,
                                    dimnames = list(c(unique(gnomad$IDs)),
                                                    c("value")) ))
  effect_size_cat <-  data.frame(matrix(nrow=length(unique(gnomad$IDs)), ncol=1,
                                        dimnames = list(c(unique(gnomad$IDs)),
                                                        c("value")) ))
  set.seed(1234)
  for(i in unique(gnomad$IDs)) {
    
    All <- gnomad %>%
      filter(IDs == "Reference"
      ) %>% 
      mutate(IDs = "Ref")
    data_plot <- gnomad %>%
      filter(IDs == i
      ) %>%
      bind_rows(., All) %>%
      select("IDs", "allele_bin", AB_frequency, AB_bin)
    
    p <- data_plot %>%
      ggplot(aes(x = AB_bin, y= AB_frequency, color= IDs))+
      geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
                  span = 0.6
      )+
      ylim(0, max(data_plot$AB_frequency))
    p_ <- layer_data(p, 1) %>%
      mutate(colour =  as.factor(colour)) %>%
      mutate_if(is.numeric, ~replace(., is.na(.), 0))
    
    variant <- p_ %>% filter(group == 1)
    variant <- sample(variant$x, 180, prob = variant$y, replace = TRUE)
    ref <- p_ %>% filter(group == 2)
    ref <- sample(ref$x, 180, prob = ref$y, replace = TRUE)
    
    result <- wilcox.test(variant, ref, 
                          alternative = "less", 
                          paired=FALSE)
    allele_balance_pval[i,] <- result$p.value
    variant_mean[i,] <- t.test(variant)$estimate
    
    result2 <- cohen.d(variant, ref, hedges.correction=TRUE)$estimate
    effect_size[i,] <- result2
    result3 <- cohen.d(variant, ref,hedges.correction=TRUE)$magnitude[1] %>% 
      as_tibble() %>%
      mutate(value = case_when(
        value == "negligible"   ~ "negligible",
        value == "small"        ~ "small",
        value == "medium"       ~ "medium",
        value == "large"        ~ "large"
      ))
    effect_size_cat[i,] <- result3
    
  }
  allele_balance_pval <- allele_balance_pval %>% 
    bind_cols(variant_mean) %>% 
    bind_cols(effect_size) %>% 
    bind_cols(effect_size_cat) %>% 
    `colnames<-`(c("allele_balance_pval", "AB_distribution_mean", "AB_effect_size", "AB_effect_size_cat")) %>%
    mutate(IDs = c(unique(gnomad$IDs)))
  
  # Bind the distribution characteristic values
  gnomad <- gnomad %>%
    full_join(., allele_balance_pval, by = "IDs") %>%
    distinct(IDs, .keep_all = TRUE) %>% 
    filter(IDs != "Reference") %>% 
    select(-c(AB_bin, allele_bin))
  
}


###################################################################### II ### Reading and Cleaning
for (i in seq_along(file_list)){
  filename <- file_list[[i]]
  
  gnomad <- read.delim(filename,
                       header = FALSE,
                       col.names = c("X.CHROM", "POS", "ID", "REF", "ALT", "INFO"))
  
  gnomad <- allele_balance(gnomad = gnomad)

  filename <- str_match(filename,
                        "/share/lab_gillis/Christelle/gnomAD_raw_data/splitted_data/(.*?)\\.vcf\\.gz")[,2]

  write_delim(gnomad, paste0(filename, "_allele_balance.vcf.gz"))
  
}

