## Import library

library(tidyverse)

###################################################################### I ### Load data

path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
gnomad <-
  read.delim(paste0(path, "/gnomAD_raw_data/dnmt3a_data.vcf.gz"), 
             header = FALSE, 
             col.names = c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

cosmic <-
  read.delim(paste0(path, "/cosmic_raw_data/cosmic_DNMT3_data.vcf.gz"), 
             header = FALSE, 
             col.names = c("Gene.name", "Accession.Number", "Gene.CDS.length", "HGNC.ID", "Sample.name", 
                           "ID_sample", "ID_tumour", "Primary.site", "Site.subtype.1", "Site.subtype.2", 
                           "Site.subtype.3", "Primary.histology", "Histology.subtype.1", "Histology.subtype.2", 
                           "Histology.subtype.3", "Genome.wide.screen", "GENOMIC_MUTATION_ID", 
                           "LEGACY_MUTATION_ID", "MUTATION_ID", "Mutation.CDS", "Mutation.AA", 
                           "Mutation.Description", "Mutation.zygosity", "LOH", "GRCh", 
                           "Mutation.genome.position", "Mutation.strand", "Resistance.Mutation", 
                           "FATHMM.prediction", "FATHMM.score", "Mutation.somatic.status", "Pubmed_PMID", 
                           "ID_STUDY", "Sample.Type", "Tumour.origin", "Age", "HGVSP", "HGVSC", "HGVSG"))
## First test HPC
head(gnomad)


###################################################################### II ### Cleaning
gnomad_decoded <- gnomad %>% 
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
  mutate(AN = str_match(INFO, "AN=(.*?);")[,2]) %>% 
  mutate(AN_samples_count = (as.numeric(AN) / 2)
  ) %>% 
  select(IDs, everything())

gnomad_decoded <- gnomad_decoded %>% 
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
  mutate(sample_frequency = sample_count / nbr_individuals) %>% 
  mutate(variant_frequency = sample_count / AN_samples_count)


## Step 1 : Limit variants with a skewed age distribution
mann_w <-  data.frame(matrix(nrow=1, ncol=0)) 
set.seed(1234)
for(i in unique(gnomad_decoded$IDs)) {
  
  All <- gnomad_decoded %>% 
    filter(IDs == "All individuals"
    ) %>% 
    mutate(IDs = "Reference")
  data_plot <- gnomad_decoded %>% 
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
  
  result <- wilcox.test(variant, ref, alternative = "greater")$p.value
  mann_w <- rbind(mann_w, result)
  
}
mann_w <- mann_w %>% `colnames<-`(c("mann_w")) %>% 
  mutate(IDs = c(unique(gnomad_decoded$IDs)))

# Bind the distribution characteristic values
selected_variants <- gnomad_decoded %>%
  full_join(., mann_w, by = "IDs") %>%
  filter(mann_w <= 0.1) %>% 
  distinct(IDs, .keep_all = TRUE)



# gnomad_decoded4 <- gnomad_decoded1

# gnomad_decoded4 <- gnomad_decoded4 %>% 
#   mutate(gene_in_gnomad = str_match(gnomad_decoded4$INFO, "(DNMT3A|TET2|TP53|ASXL1|PPM1D)(.*?)$")[,2])

# data_p <- gnomad_decoded4 %>% 
#   bind_rows(., gnomad_decoded1 %>% filter(IDs == "All individuals"))



## Step 2 : Limit variants present in COSMIC
cosmic <-
  cosmic %>% 
  # separate(HGVSG, into = c("chr", "chr_arm", "POS", "ALT"), sep = ":|\\.|>|delins|dup|del", remove = FALSE) %>% 
  # mutate(REF = str_match(POS, "([0-9]*)([A-Z]*)")[,3]) %>% 
  # mutate(POS = str_match(POS, "([0-9]*)([A-Z]*)")[,2],
  #        POS = as.numeric(POS)) %>% 
  # # separate(Mutation.genome.position, into = c("chr", "start", "end"), sep = ":|-", remove = FALSE) %>% 
  # unite(IDs, c(chr, POS, REF, ALT), sep = "-", remove = FALSE) %>% #select(IDs, everything())
  # # separate(HGVSG, into = c("ENSP", "p_prot_modif"), sep = ":", remove = FALSE) %>% 
  # mutate(IDs = case_when(
  #   IDs == "-NA-NA-NA"                 ~ NA_character_,
  #   TRUE                               ~ IDs
  # )) %>% 
  # separate(`Gene.name`, into = c("gene_name_cosmic"), sep = "_", remove = FALSE) %>% 
  # # 1 tumor can have multiple sample
  # distinct(IDs, ID_tumour, .keep_all = TRUE) %>% 
  # arrange(IDs) %>% 
  # select(IDs, "Gene.name", gene_name_cosmic, Gene.name, "Accession.Number", "chr", "POS", REF, "ALT", HGVSG, everything())

final_cosmic <- cosmic %>% 
  group_by(IDs, gene_name_cosmic, chr, POS, REF, ALT) %>% 
  mutate(occurrence_in_cosmic = n())

CH_variants_cpra <- inner_join(selected_variants %>% 
                                 select(-c(AN:variant_frequency)), 
                               final_cosmic, 
                               by = c("IDs", "X.CHROM" = "chr", "POS", "REF", "ALT"))

write_rds(CH_variants_cpra, "CH_variants_cpra.rds")


###################################################################### III ### Extract variant data
# Consequences
# CH_variants <- CH_variants_cpra %>%
#   mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2])
# 
# CH_variants <- CH_variants %>%
#   separate(col = ens_vep, paste("ens_vep", 1:25, sep=""),
#            sep = ",", remove = T, extra = "warn", fill = "right") %>%
#   keep(~!all(is.na(.))) %>%
#   pivot_longer(cols = starts_with("ens_vep"), names_to = NULL, values_to = "ens_vep") %>%
#   drop_na(ens_vep)
# 
# ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "
#   cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "
#   STRAND", "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "
#   TREMBL", "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET", "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "
#   SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF", "ExAC_EAS_MAF", "ExAC_FIN_MAF", "
#   ExAC_NFE_MAF", "ExAC_OTH_MAF", "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "
#   HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter", "LoF_flags", "LoF_info")
# 
# CH_variants <- CH_variants %>%
#   separate(col = ens_vep, into = ens_vep_var_names,
#            # paste("ens_vep", 1:100, sep=""),
#            sep = "\\|", remove = F, extra = "warn", fill = "right")


# Allele info
CH_variants <- CH_variants_cpra %>%
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
  
  # Depth of informative coverage for each sample
  mutate(depth_allele = str_match(INFO, "DP=(.*?);")[,2]) %>%
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
  mutate(nc_freq_allele_count_raw = str_match(INFO, "non_cancer_AF_raw=(.*?);")[,2])

write_rds(CH_variants, "CH_variants.rds")












