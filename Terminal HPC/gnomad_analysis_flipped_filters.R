## Import library
library(tidyverse)
library(gtsummary)
library(ggbreak)
library(gtsummary)
library(ggconsort)
library(viridis)
library(ggridges)
library(ggforce)

theme_gtsummary_compact()
theme_set(theme_classic(base_size = 20))
options(bitmapType='cairo')

###################################################################### I ### List files
file_list <- list.files(
  path =
    "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants/allele_balance/age_distribution/variant_depth/noise_filter/padjust/potential_ch_variants",
  pattern = "*.vcf.gz",
  recursive=FALSE,
  full.names = TRUE)

gnomad_variants <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " "))

vep_bick_genes <-
  readxl::read_xlsx(
    paste0(here::here(), 
           "/Christelle/gnomAD_skweness/analysis/data/gnomad_clonal_mosaicism_variants_list_w_all_annotation_L_CHIP_09.12.2023.xlsx"), 
    na = c(".", "NA"))

considered_CH_variants <- 
  readxl::read_xlsx(paste0(here::here(), "/Christelle/gnomAD_skweness/analysis/data/gene_ids.xlsx")) %>% 
  select(ids, Gene, Variant) %>% 
  tidyr::unite(col = "considered_variants", c(Gene : Variant), sep = " ")


## Cleaning
vep_bick_genes <- vep_bick_genes %>% 
  mutate(mutations_in_BickWHO = case_when(
    In_Bick_WHO == "TRUE"              ~ "Yes",
    TRUE                               ~ NA_character_
  )) %>% 
  mutate(L_CHIP = case_when(
    L_CHIP == "TRUE"                   ~ "Yes",
    TRUE                               ~ NA_character_
  )) %>% 
  mutate(Polyphen2_HDIV_pred = case_when(
    Polyphen2_HDIV_pred == "D"     ~ "Probably damaging (>=0.957)",
    Polyphen2_HDIV_pred == "P"     ~ "Possibly damaging (0.453<=pp2_hdiv<=0.956)",
    Polyphen2_HDIV_pred == "B"     ~ "Benign (pp2_hdiv<=0.452)"
  )) %>% 
  mutate(SIFT_pred = case_when(
    SIFT_pred == "D"               ~ "Deleterious",
    SIFT_pred == "T"               ~ "Tolerated"
  )) %>% 
  mutate(FATHMM_pred = case_when(
    FATHMM_pred == "D"             ~ "Deleterious",
    FATHMM_pred == "T"             ~ "Tolerated"
  )) %>% 
  mutate(CADD_phred_cat = case_when(
    CADD_phred >= 20               ~ "Likely pathogenic",
    CADD_phred < 20                ~ "< 20"
  )) %>% 
  mutate(new_HGVSc_2 = str_match(AAChange.refGene, "c.(.)([:alnum:]*)(.):p.")[,2]) %>% 
  mutate(new_HGVSc_4 = str_match(AAChange.refGene, "c.(.)([:alnum:]*)(.):p.")[,4]) %>% 
  unite(new_HGVSc, c(new_HGVSc_2, new_HGVSc_4), sep = ">", na.rm = TRUE) %>%
  mutate(new_HGVSc = case_when(
    new_HGVSc == ""                ~ NA_character_,
    TRUE                           ~ new_HGVSc
  )) %>% 
  select(-c(In_Bick_WHO, X.CHROM))

## Merge

gnomad_variants <- gnomad_variants %>% 
  # Create filter for removing noise in AB
  # mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
  # separate(col = ab_hist_alt_bin_freq,
  #          into = c("0.00-0.05"),
  #          sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  # mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
  # mutate(contaminated_AB = case_when(
  #   `0.00-0.05` == total_allele_balance               ~ "sequencing noise",
  #   TRUE                                              ~ "clean sequencing"
  # )) %>% 
  # select(-c("0.00-0.05")) %>% 
  # Add VEP and Bick&WHO
  left_join(., vep_bick_genes,
            by= "IDs") %>% 
  mutate(mutations_in_BickWHO = case_when(
    mutations_in_BickWHO == "Yes"              ~ "Yes",
    TRUE                                       ~ "No"
  ), mutations_in_BickWHO = factor(mutations_in_BickWHO, levels = c("No", "Yes"))) %>% 
  mutate(L_CHIP = case_when(
    L_CHIP == "Yes"                            ~ "Yes",
    TRUE                                       ~ "No"
  ), L_CHIP = factor(L_CHIP, levels = c("No", "Yes"))) 

gnomad_variants <- gnomad_variants %>% 
  left_join(., considered_CH_variants,
            by = c("IDs" = "ids")) #%>% 
  # mutate(L_CHIP = case_when(
  #   SYMBOL %in% c("DUSP22", "FAT1", "KMT2D",
  #                 "SYNE1", "ATM", "KMT2C",
  #                 "PCLO", "SPEN", "ARID1A", "NEB",
  #                 "MGA", "RP1L1", "SMARCA4", "FAT2",
  #                 "IGLL5", "CIITA", "DLC1", "MTOR",
  #                 "MYD88", "PTPRD", "SOCS1", "ATR",
  #                 "FAT4", "HIST1H2BK", "NOTCH2")       ~ "Yes",
  #   TRUE                                               ~ "No"
  # ))

# cleaning
rm(considered_CH_variants,
   vep_bick_genes)

object.size(gnomad_variants)
gc()
###################################################################### II ### Apply filter----
print(paste0("The number of missing age in gnomad is ", 
             sum(is.na(gnomad_variants %>% distinct(IDs, .keep_all = TRUE) %>% 
                         select(adjusted_age_pval)))))
print("Number of variants present in cosmic")
nrow(gnomad_variants %>% distinct(IDs, .keep_all = TRUE) %>% filter(variant_in_cosmic == "Yes"))

print("Number of cosmic samples present in gnomad")
nrow(gnomad_variants %>% filter(variant_in_cosmic == "Yes") %>% distinct(Sample.name, .keep_all = TRUE))

# u2af1 <- gnomad_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   filter(IDs == "21-43104346-G-A")
# write_csv(u2af1, "u2af1.csv")
# asxl1 <- gnomad_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   filter(SYMBOL == "ASXL1")
# write_csv(asxl1, "asxl1.csv")
print("table considered_variants in gnomad")
table(gnomad_variants$considered_variants)

clean_AB_variants <- gnomad_variants %>%
  filter(contaminated_AB == "clean sequencing") %>%
  filter(black_list == 0) %>%
  filter(centromeres == 0) %>%
  filter(segmental_duplication == 0) %>%
  filter(wm_sdust == 0)

print("table considered_variants in clean_AB_variants")
table(clean_AB_variants$considered_variants)
allele_depth <- 10
depth_variants <- clean_AB_variants %>%
  filter(dp_het_median >= allele_depth)

print("table considered_variants in depth_variants")
table(depth_variants$considered_variants)
prot_coding <- depth_variants %>%
  filter(is_gencode_protein_coding_gene == "Yes")

print("table considered_variants in prot_coding")
table(prot_coding$considered_variants)
pval_selection <- 0.05
AB_distribution_variants <- prot_coding %>% 
  filter(adjusted_AB_pval < pval_selection)
AB_hedges_selection <- "large"
clonal_mosaicism <- AB_distribution_variants %>%
  # filter(effect_size_cat != "negligible")
  filter(AB_effect_size_cat == AB_hedges_selection)
print(paste0("The number of missing age in CM is ", 
             sum(is.na(clonal_mosaicism %>% distinct(IDs, .keep_all = TRUE) %>% ######################## NEW
                         select(adjusted_age_pval)))))

print("known Bick+WHO variants in clonal_mosaicism")
clonal_mosaicism %>% filter(mutations_in_BickWHO == "Yes") %>%
  select(IDs, considered_variants) %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(considered_variants))
write_csv(clonal_mosaicism %>% filter(mutations_in_BickWHO == "Yes") %>% 
            select(IDs, considered_variants, SYMBOL, 
                   adjusted_age_pval, age_effect_size_cat, 
                   AAChange.refGene) %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(!is.na(considered_variants)),
          "known Bick+WHO variants in clonal_mosaicism.csv")
print("table considered_variants in clonal_mosaicism")
clonal_mosaicism %>% distinct(IDs, .keep_all = TRUE) %>%
  filter(!is.na(considered_variants)) %>% 
  select(considered_variants) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**considered_variants in SM**"
  ) %>% as_kable()
# table(clonal_mosaicism$considered_variants)
print("SM gene that are L-CHIP")
a <- clonal_mosaicism %>% filter(L_CHIP == "Yes") %>% distinct(IDs, .keep_all = TRUE)
a %>% select(SYMBOL) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**L-CHIP in SM**"
  ) %>% as_kable()

write_csv(a %>% 
  select(SYMBOL, alt_allele_count_nfe, 
         alt_allele_count_afr, 
         alt_allele_count_amr,
         alt_allele_count_eas,
         alt_allele_count_sas,
         alt_allele_count,
         nc_alt_allele_count_nfe, 
         nc_alt_allele_count_afr, 
         nc_alt_allele_count_amr,
         nc_alt_allele_count_eas,
         nc_alt_allele_count_sas,
         nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count)),
  "SM gene that are L-CHIP.csv")

write_csv(a %>% 
            select(SYMBOL, nc_freq_allele_count_nfe, 
                   nc_freq_allele_count_afr, 
                   nc_freq_allele_count_amr,
                   nc_freq_allele_count_eas,
                   nc_freq_allele_count_sas,
                   nc_freq_allele_count) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            arrange(desc(nc_freq_allele_count)),
          "SM cn_freq of gene that are L-CHIP.csv")

age_distribution_variants <- clonal_mosaicism %>% 
  filter(adjusted_age_pval < pval_selection)
age_hedges_selection <- "negligible|small"
CH_variants <- age_distribution_variants %>%
  filter(!str_detect(age_effect_size_cat, age_hedges_selection))
print("known Bick+WHO variants in CH_variants")
CH_variants %>% filter(mutations_in_BickWHO == "Yes") %>%
  select(IDs, considered_variants) %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(considered_variants))
write_csv(CH_variants %>% filter(mutations_in_BickWHO == "Yes") %>% 
            select(IDs, considered_variants, SYMBOL, AAChange.refGene) %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(!is.na(considered_variants)),
          "known Bick+WHO variants in CH_variants.csv")
print("table considered_variants in CH_variants")
CH_variants %>% distinct(IDs, .keep_all = TRUE) %>%
  filter(!is.na(considered_variants)) %>% 
  select(considered_variants) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**considered_variants in SM**"
  ) %>% as_kable()
# table(CH_variants$considered_variants)
print("CH gene that are L-CHIP")
a <- CH_variants %>% filter(L_CHIP == "Yes") %>% distinct(IDs, .keep_all = TRUE)
table(as.character(a$SYMBOL))
a %>% select(SYMBOL) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**L-CHIP in CH**"
  ) %>% as_kable()

write_csv(a %>% 
            select(SYMBOL, alt_allele_count_nfe, 
                   alt_allele_count_afr, 
                   alt_allele_count_amr,
                   alt_allele_count_sas,
                   alt_allele_count_eas,
                   alt_allele_count,
                   nc_alt_allele_count_nfe, 
                   nc_alt_allele_count_afr, 
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            arrange(desc(nc_alt_allele_count)),
          "CH gene that are L-CHIP.csv")

write_csv(a %>% 
            select(SYMBOL, nc_freq_allele_count_nfe, 
                   nc_freq_allele_count_afr, 
                   nc_freq_allele_count_amr,
                   nc_freq_allele_count_eas, 
                   nc_freq_allele_count_sas,
                   nc_freq_allele_count) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            arrange(desc(nc_freq_allele_count)),
          "CH cn_freq of gene that are L-CHIP.csv")

cosmic_var <- CH_variants %>% 
  filter(variant_in_cosmic == "Yes")
print("table considered_variants in cosmic_var")
table(cosmic_var$considered_variants)

print(paste("Final number of variants in all chromosome are",
            nrow(cosmic_var %>% distinct(IDs))))
# print(tail(gnomad_variants,2))
# write_delim(CH_variants %>% distinct(IDs, .keep_all = TRUE) %>% select(IDs, "X.CHROM"),
#       "gnomad_CH_variants_list.txt")
# write_delim(cosmic_var %>% distinct(IDs, .keep_all = TRUE) %>% select(IDs, "X.CHROM"),
#           "gnomad_variants_list_after_filtering_without_age_filter.txt")
# write_delim(clonal_mosaicism %>% distinct(IDs, .keep_all = TRUE) %>% select(IDs, "X.CHROM"),
#             "gnomad_clonal_mosaicism_variants_list.txt")

###################################################################### II ### Analysis----
print("FIGURE 1_A")
# # CONSORT FIPPED----
# study_cohorts <-
#   gnomad_variants %>%
#   distinct(IDs, .keep_all = TRUE) %>%
#   cohort_start("**Variants considered**") %>%
#   # Define cohorts using named expressions
#   # Notice that you can use previously defined cohorts in subsequent steps
#   cohort_define(
#     
#     clean_seq = .full %>%
#       filter(contaminated_AB == "clean sequencing"),
#     black_clean = clean_seq %>%
#       filter(black_list == 0),
#     centro_clean = black_clean %>%
#       filter(centromeres == 0),
#     dup_clean = centro_clean %>%
#       filter(segmental_duplication == 0),
#     wm_sdust_clean = dup_clean %>%
#       filter(wm_sdust == 0),
#     clean_AB_variants = wm_sdust_clean %>%
#       filter(dp_het_median >= allele_depth),
#     
#     protcoding_variants = clean_AB_variants  %>%
#       filter(is_gencode_protein_coding_gene == "Yes"),
#     
#     AB_BH = protcoding_variants %>%
#       filter(adjusted_AB_pval <= pval_selection),
#     clonal_mosaicism = AB_BH %>%
#       filter(AB_effect_size_cat == AB_hedges_selection),
#     
#     age_BH = clonal_mosaicism  %>%
#       filter(adjusted_age_pval <= pval_selection),
#     CH_variants = age_BH %>%
#       filter(!str_detect(age_effect_size_cat, age_hedges_selection)),
#     
#     cosmic_variants = CH_variants  %>%
#       filter(variant_in_cosmic == "Yes"),
#     
#     CH_variants_end = cosmic_variants,
#     
#     # anti_join is useful for counting exclusions
#     excluded2 = anti_join(.full, clean_AB_variants, by = "IDs"),
#     excluded2_1 = anti_join(.full, clean_seq, by = "IDs"),
#     excluded2_2 = anti_join(clean_seq, black_clean, by = "IDs"),
#     excluded2_3 = anti_join(black_clean, centro_clean, by = "IDs"),
#     excluded2_4 = anti_join(centro_clean, dup_clean, by = "IDs"),
#     excluded2_5 = anti_join(dup_clean, wm_sdust_clean, by = "IDs"),
#     excluded2_6 = anti_join(wm_sdust_clean, clean_AB_variants, by = "IDs"),
#     
#     excluded3 = anti_join(clean_AB_variants, protcoding_variants, by = "IDs"),
#     
#     excluded1 = anti_join(protcoding_variants, clonal_mosaicism, by = "IDs"),
#     excluded1_1 = anti_join(protcoding_variants, AB_BH, by = "IDs"),
#     excluded1_2 = anti_join(AB_BH, clonal_mosaicism, by = "IDs"),
#     
#     excluded4 = anti_join(clonal_mosaicism, CH_variants, by = "IDs"),
#     excluded4_1 = anti_join(clonal_mosaicism, age_BH, by = "IDs"),
#     excluded4_2 = anti_join(age_BH, CH_variants, by = "IDs"),
#     
#     excluded5 = anti_join(CH_variants, cosmic_variants, by = "IDs")
#   ) %>%
#   # Provide text labels for cohorts
#   cohort_label(
#     excluded2 = "<span style='color:darkblue'>Exclude variants with sequencing artifacts</span>",
#     excluded2_1 = "Exclude variants with sequencing artifacts",
#     excluded2_2 = "Exclude 'curated “blacklist” developed by ENCODE (73), and excluded mutations in the blacklisted regions'",
#     excluded2_3 = "Exclude 'Centromeres (“Centromeres” in the “Mapping and Sequencing” group)'",
#     excluded2_4 = "Exclude 'Segmental duplications (“segmental dups”)'",
#     excluded2_5 = "Exclude 'Low complexity regions (“WM + SDust”) track'",
#     excluded2_6 = "Exclude variants with depth inferieur to 10",
#     excluded3 = "<span style='color:darkblue'>Select for protein coding mutations</span>",
#     excluded1 = "<span style='color:darkblue'>Exclude germline variants</span>",
#     excluded1_1 = "BH pval < 0.05",
#     excluded1_2 = "Hedges g == 'large'",
#     excluded5 = "<span style='color:darkblue'>Select variants present in Cosmic</span>",
#     excluded4 = "<span style='color:darkblue'>Exclude non-age-skewed mutations</span>",
#     excluded4_1 = "BH pval < 0.05",
#     excluded4_2 = "Hedges g == 'large' or 'medium'",
#     clean_AB_variants = "**Variants after quality control filtering**",
#     protcoding_variants = "**Protein-coding somatic mutations**",
#     clonal_mosaicism = "**Peripheral blood somatic mosaicism (SM) variants (clonal_mosaicism)**",
#     cosmic_variants = "**CH variants observed in tumors**",
#     CH_variants = "**CH variants**",
#     CH_variants_end = "**?**"
#   )
# 
# study_consort <- study_cohorts %>%
#   consort_box_add(
#     # top eligibility box at 40 height
#     "full", 0.05, 50, cohort_count_adorn(study_cohorts, .full)
#   ) %>%
#   # first top excluded box at 0.2 (right from 0.05 eligibility box)
#   consort_box_add(
#     "exc_contamination", 0.2, 45, glue::glue(
#       "{cohort_count_adorn(study_cohorts, excluded2)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_1)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_2)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_3)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_4)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_5)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded2_6)}
#       ")
#   ) %>%
#   consort_box_add(
#     "clean_AB_variants", 0.05, 40, cohort_count_adorn(study_cohorts, clean_AB_variants)
#   ) %>%
#   consort_box_add(
#     "exc_non_protein_coding", 0.2, 35, glue::glue(
#       "{cohort_count_adorn(study_cohorts, excluded3)}<br>
#       • include variants present in a protein coding gene
#       ")
#   ) %>%
#   consort_box_add(
#     "protcoding_variants", 0.05, 30, cohort_count_adorn(study_cohorts, protcoding_variants)
#   ) %>%
#   consort_box_add(
#     "exc_germline", 0.2, 25, glue::glue(
#       "{cohort_count_adorn(study_cohorts, excluded1)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded1_1)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded1_2)}
#       ")
#   ) %>%
#   consort_box_add(
#     "clonal_mosaicism", 0.05, 20, cohort_count_adorn(study_cohorts, clonal_mosaicism)
#   ) %>%
#   consort_box_add(
#     "exc_non_age_skew", 0.2, 15, glue::glue(
#       "{cohort_count_adorn(study_cohorts, excluded4)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded4_1)}<br>
#       • {cohort_count_adorn(study_cohorts, excluded4_2)}
#       ")
#   ) %>%
#   # Add bottom box
#   consort_box_add(
#     "CH_variants", 0.05, 10, cohort_count_adorn(study_cohorts, CH_variants)
#   ) %>%
#   consort_box_add(
#     "exc_non_cosmic", 0.2, 5, glue::glue(
#       "{cohort_count_adorn(study_cohorts, excluded5)}<br>
#       • Select variants present in Cosmic
#       ")
#   ) %>%
#   consort_box_add(
#     "cosmic_variants", 0.05, 0, cohort_count_adorn(study_cohorts, cosmic_variants)
#   ) %>%
#   # Add exclusion arrows
#   consort_arrow_add(
#     end = "exc_contamination", end_side = "left", start_x = 0.05, start_y = 45
#   ) %>%
#   consort_arrow_add(
#     end = "exc_non_protein_coding", end_side = "left", start_x = 0.05, start_y = 35
#   ) %>%
#   consort_arrow_add(
#     end = "exc_germline", end_side = "left", start_x = 0.05, start_y = 25
#   ) %>%
#   consort_arrow_add(
#     end = "exc_non_cosmic", end_side = "left", start_x = 0.05, start_y = 5
#   ) %>%
#   consort_arrow_add(
#     end = "exc_non_age_skew", end_side = "left", start_x = 0.05, start_y = 15
#   ) %>%
#   # Add top to bottom arrow
#   consort_arrow_add(
#     end = "cosmic_variants", end_side = "top", start_x = 0.05, start_y = 50
#   )
# 
# jpeg("Figure1A consort diagram all chromosomes_flipped.jpeg", width = 950, height = 650)
# study_consort %>%
#   ggplot() +
#   geom_consort() + #theme_classic()
#   theme_consort(margin_h = c(1,12), # bottom, left
#                 margin_v = c(1,50)) # top, right
# dev.off()

print("FIGURE 1_B")
# Fig1_B geom_area filter----

fig1_b <- gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, AB_distribution_mean) %>%
  mutate(data = "Variants considered") %>%

  bind_rows(., prot_coding %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "Variants in protein-coding genes")) %>%
  bind_rows(., clonal_mosaicism %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "Peripheral blood somatic mosaicism (SM) variants")) %>%
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  bind_rows(., cosmic_var %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "CH variants observed in tumors")) %>%

  mutate(AB_distribution_mean = round(AB_distribution_mean, 2)) %>%
  filter(!is.na(AB_distribution_mean)) %>%
  group_by(AB_distribution_mean, data) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(data = factor(data, levels= c("Variants considered",
                                       "Variants in protein-coding genes",
                                       "Peripheral blood somatic mosaicism (SM) variants",
                                       "Clonal hematopoiesis (CH) variants",
                                       "CH variants observed in tumors"
  )))

pdf("Figure1B_all_chromosomes_part1_flipped.pdf")
fig1_b %>%
  ggplot(aes(x= AB_distribution_mean, y=count, fill= data
  ))+
  geom_area(position = "identity") +
  # scale_fill_viridis(discrete = TRUE, name = NULL, option = "A")+
  scale_fill_manual(values = c("#440154FF", "#3B528BFF",
                               "#21908CFF","#5DC863FF", "#FDE725FF"
                               ),
                    name= NULL)+
  labs(x= paste("Average Allele Balance Distribution"),
       y= "Number of Variants"
  )+
  theme(legend.position = "none")
dev.off()

# extract legend
legend <- fig1_b %>%
  ggplot(aes(x= AB_distribution_mean, y=count, fill= data
  ))+
  geom_area(position = "identity") +
  # scale_fill_viridis(discrete = TRUE, name = NULL, option = "A")+
  scale_fill_manual(values = c("#440154FF", "#3B528BFF",
                               "#21908CFF","#5DC863FF", "#FDE725FF"
  ),
  name= NULL)+
  labs(x= paste("Average Allele Balance Distribution"),
       y= "Number of Variants"
  )+
  theme(legend.position = "bottom")+ 
  guides(fill = guide_legend(nrow=3, byrow=TRUE))
legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
pdf("Figure1B_legend_bottom.pdf", width = 20, height = 3)
ggpubr::as_ggplot(legend)
dev.off()

fig1_b <- clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, AB_distribution_mean) %>%
  mutate(data = "Peripheral blood somatic mosaicism (SM) variants") %>%
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  bind_rows(., cosmic_var %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, AB_distribution_mean) %>%
              mutate(data = "CH variants observed in tumors")) %>%

  mutate(AB_distribution_mean = round(AB_distribution_mean, 2)) %>%
  filter(!is.na(AB_distribution_mean)) %>%
  group_by(AB_distribution_mean, data) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(data = factor(data, levels= c("Peripheral blood somatic mosaicism (SM) variants",
                                       "Clonal hematopoiesis (CH) variants",
                                       "CH variants observed in tumors"
  )))

a <- fig1_b %>%
  ggplot(aes(x= AB_distribution_mean, y=count, fill= data
  ))+
  geom_area(position = "identity") +
  # scale_fill_viridis(discrete = TRUE, name = NULL, option = "A")+
  scale_fill_manual(values = c(#"#440154FF", "#3B528BFF",
                               "#21908CFF","#5DC863FF", "#FDE725FF"),
                    name= NULL)+
  labs(x= paste("Average Allele Balance Distribution"),
       y= "Number of variants"
  )+
  theme(legend.position = "none",
        axis.title.y = element_blank())
pdf("Figure1B_all_chromosomes_part2_cutaxis.pdf", width = 4)
a
dev.off()

# Figure1_C----
print("FIGURE S1")
pdf("FigureS1_packed_circle.pdf")
library(ggraph)
library(igraph)

cat_to_cat <- tibble(from= c("gnomad_variants", "prot_coding", "clonal_mosaicism", "CH_variants",
                             "gnomad_variants", "prot_coding", "clonal_mosaicism", "CH_variants"),
                     to= c("prot_coding", "clonal_mosaicism", "CH_variants", "cosmic_var", 
                           "Lost_incleaning", "lost_inAB", "lost_inAge", "lost_inCosmic"))

size_plot <- tibble(name= c("gnomad_variants", "prot_coding", "clonal_mosaicism", "CH_variants", "cosmic_var",
                            "Lost_incleaning", "lost_inAB", "lost_inAge", "lost_inCosmic"),
                    size= c(3934541, 2684332, 503703, 89361, 3834,
                            1250209, 2180629, 414342, 85527),
                    shortName= c("Variants considered", "Variants in protein-coding genes",
                                 "Peripheral blood somatic mosaicism (SM) variants ", 
                                 "Clonal hematopoiesis (CH) variants",
                                 "CH variants observed in tumors",
                                 "Lost_incleaning", "lost_inAB", "lost_inAge", "lost_inCosmic"))
mygraph1 <- graph_from_data_frame(cat_to_cat, vertices=size_plot)

ggraph(mygraph1, layout = 'circlepack', weight=size) +
  geom_node_circle(data= . %>% filter(shortName %in% c("Variants considered", "Variants in protein-coding genes",
                                                       "Peripheral blood somatic mosaicism (SM) variants ",
                                                       "Clonal hematopoiesis (CH) variants", "CH variants observed in tumors")),
                   aes(#x0= 0, 
                     fill= factor(shortName,
                                  levels= c("Variants considered", "Variants in protein-coding genes",
                                            "Peripheral blood somatic mosaicism (SM) variants ",
                                            "Clonal hematopoiesis (CH) variants", "CH variants observed in tumors"
                                  )
                     ))) +
  scale_fill_manual(values = c("#440154FF",  "#3B528BFF",  "#21908CFF",  "#5DC863FF",  "#FDE725FF"),
                    name = NULL)+
  theme(legend.position = "none")
dev.off()

print("FIGURE 2")
# Fig2_A AB / Age distribution in CH variants----
extract_AB <- CH_variants %>%
  filter(!is.na(considered_variants)) %>%
  # filter(IDs %in% c("2-25234373-C-T", "2-25235726-A-G",
  #                   "2-197402110-T-C")) %>%
  distinct(IDs, .keep_all = TRUE)

print("Adjusted FDR for CH considered_variants")
CH_variants %>%
  filter(!is.na(considered_variants)) %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(considered_variants, adjusted_AB_pval, adjusted_age_pval)
write_csv(clonal_mosaicism %>%
            filter(!is.na(considered_variants)) %>%
            distinct(IDs, .keep_all = TRUE), "considered_variants in SM.csv")

gnomad <- extract_AB %>%
  filter(IDs != "All individuals") %>%
  select(-AB_count) %>%
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
          "0.90-0.95" = 92, "0.95-1.00" = 0, .before = 1
  ) %>%
  mutate(total_allele_balance = rowSums(select(.,`0.00-0.05`:`0.95-1.00`), na.rm = TRUE)) %>%

  select(IDs, everything())

gnomad <- gnomad %>%
  pivot_longer(cols = c(grep("^[[:digit:]]", colnames(.))),
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

All <- gnomad %>%
  filter(IDs == "Reference"
  ) %>%
  mutate(IDs = "Ref")

pdf("Fig3_reference_AB_distribution.pdf")
All %>%
  ggplot(aes(x = AB_bin, y= AB_frequency))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              span = 0.6, color= "#3333FF"
  )+
  labs(title = All$ID)+
  labs(x= "Allele balance", y= "Allele balance frequency")+
  geom_col(data = All ,
           aes(x= AB_bin, y= AB_frequency),
           position = "identity",
           fill = "#6699FF",
           color= "#3333FF",
           alpha = 0.1)+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 24))+
  ylim(0, 0.4)
dev.off()

data_plot <- gnomad %>%
  filter(IDs == "2-25234373-C-T"
  ) %>%
  bind_rows(., All) %>%
  distinct(IDs, AB_bin, .keep_all = TRUE) %>% 
  mutate(IDs = case_when(
    IDs == "2-25234373-C-T"     ~ "Observed distribution",
    TRUE                        ~ "Expected distribution"
  ))

p <- data_plot %>%
  ggplot(aes(x = AB_bin, y= AB_frequency, color= IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              span = 0.6
  )+
  scale_color_manual(values = c("#3333FF", "#F8766D"))+
  geom_col(data = data_plot,
           aes(x= AB_bin, y= AB_frequency, fill = IDs),
           position = "identity",
           alpha = 0.1)+
  scale_fill_manual(values = c("#6699FF", "red"))+
  theme(legend.position = "bottom",
        legend.title = element_blank())

legend <- ggpubr::get_legend(p)
# Convert to a ggplot and print
pdf("Figure2_legend_bottom.pdf", width = 20, height = 3)
ggpubr::as_ggplot(legend)
dev.off()

set.seed(1234)
for(i in unique(gnomad$IDs)) {

  if (i != "Reference") {
    # All <- gnomad %>%
    #   filter(IDs == "Reference"
    #   ) %>%
    #   mutate(IDs = "Ref")
    data_plot <- gnomad %>%
      filter(IDs == i
      ) %>%
      bind_rows(., All) %>%
      distinct(IDs, AB_bin, .keep_all = TRUE)

    plot_title <- data_plot %>%
      filter(IDs == i) %>%
      select(considered_variants) %>% distinct() %>% as.character()
    # n <- n + 1
    p <- data_plot %>%
      ggplot(aes(x = AB_bin, y= AB_frequency, color= IDs))+
      geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
                  span = 0.6
      )+
      scale_color_manual(values = c("#F8766D", "#3333FF"))+
      labs(title = plot_title)+
      labs(x= "Allele balance", y= "Allele balance frequency")+
      geom_col(data = data_plot %>%
                 filter(IDs == i) ,
               aes(x= AB_bin, y= AB_frequency),
               position = "identity",
               fill = "red",
               alpha = 0.1)+
      theme(axis.text = element_text(size = 24),
            axis.title = element_text(size = 26),
            plot.title = element_text(size = 24))+
      ylim(0, 0.4)+
      theme(legend.position = "none")

    p <- p+
      geom_text(
        data = data_plot %>%
          filter(IDs == i),
        mapping = aes(x = 0.60, y = 0.4,
                      label = paste0("FDR q = ", round(
                        adjusted_AB_pval, 5)
                        )),
        colour = "black", size = 20/.pt,
        hjust   = 0
        # vjust   = -45 # -28 would be lower
      )+
      geom_text(
        data = data_plot %>%
          filter(IDs == i), #base_size = text.size,
        mapping = aes(x = 0.60, y = 0.37,
                      label = paste0("Hedges g = ", round(AB_effect_size, 2))),
        colour = "black", size = 20/.pt, family = "serif",
        hjust   = 0
        # vjust   = -43 # -28 would be lower
      )

    pdf(paste0("Figure3_AB_for_", i, ".pdf"))
    print(p)
    dev.off()

  }
}

extract_Age <- CH_variants %>%
  filter(!is.na(considered_variants)) %>%
  # filter(IDs %in% c("2-25234373-C-T", "2-25235726-A-G",
  #                   "2-197402110-T-C")) %>%
  distinct(IDs, .keep_all = TRUE)

extract_Age <- extract_Age %>%

  select(-sample_frequency, -sample_count, -nbr_individuals) %>%

  # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
  mutate("<30" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
  mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
  separate(col = age_hist_het_bin_freq,
           into = c("30-35","35-40","40-45","45-50",
                    "50-55","55-60","60-65","65-70",
                    "70-75","75-80"),
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
  mutate(">80" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>%
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%

  add_row("IDs"="All individuals", "ID"="All age data in gnomAD (exons)",
          "<30" = 2547, "30-35" = 3423,
          "35-40" = 4546,"40-45" = 8487,
          "45-50" = 10355,"50-55" = 12693,
          "55-60" = 11933,"60-65" = 10534,
          "65-70" = 8882,"70-75" = 5991,
          "75-80" = 4136,
          ">80" = 1935, .before = 1
  ) %>%
  mutate(nbr_individuals = rowSums(select(.,`<30`:`>80`), na.rm = TRUE))

extract_Age <- extract_Age %>%
  pivot_longer(cols = c(grep("^([[:digit:]]|<|>)", colnames(.))),
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
  mutate(age_frequency = sample_count / nbr_individuals)

All <- extract_Age %>%
  filter(IDs == "All individuals"
  )

pdf("Fig3_reference_Age_distribution.pdf")
All %>%
  ggplot(aes(x = bin_age, y= age_frequency))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              span = 0.6, color= "#3333FF"
  )+
  labs(title = All$ID)+
  labs(x= "Age", y= "Age balance frequency")+
  geom_col(data = All ,
           aes(x= bin_age, y= age_frequency),
           position = "identity",
           fill = "#6699FF",
           color= "#3333FF",
           alpha = 0.1)+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 24))+
  ylim(0, 0.4)
dev.off()

set.seed(1234)
for(i in unique(extract_Age$IDs)) {

  if (i != "All individuals") {
    # All <- extract_Age %>%
    #   filter(IDs == "All individuals"
    # ) %>%
    #   mutate(IDs = paste0("Reference= ", ID))
    data_plot <- extract_Age %>%
      filter(IDs == i
      ) %>%
      bind_rows(., All) %>%
      distinct(IDs, bin_age, .keep_all = TRUE)

    plot_title <- data_plot %>%
      filter(IDs == i) %>%
      select(considered_variants) %>% distinct() %>% as.character()
    # n <- n + 1
    p <- data_plot %>%
      ggplot(aes(x = bin_age, y= age_frequency, color= IDs))+
      geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
                  span = 0.6
      )+
      scale_color_manual(values = c("#F8766D", "#3333FF"))+
      labs(title = plot_title)+
      labs(x= "Age", y= "Age frequency")+
      geom_col(data = data_plot %>%
                 filter(IDs == i) ,
               aes(x= bin_age, y= age_frequency),
               position = "identity",
               fill = "red",
               alpha = 0.1)+
      theme(axis.text = element_text(size = 24),
            axis.title = element_text(size = 26),
            plot.title = element_text(size = 24))+
      ylim(0, 0.4)+
      theme(legend.position = "none")

    p <- p+
      geom_text(
        data    = data_plot %>%
          filter(IDs == i),
        mapping = aes(x = 60, y = 0.4,
                      label = paste0("FDR q = ", round(
                        adjusted_age_pval, 5)
                        )),
        colour = "black", size = 20/.pt,
        hjust   = 0,
        # vjust   = -28
      )+
      geom_text(
        data = data_plot %>%
          filter(IDs == i), #base_size = text.size,
        mapping = aes(x = 60, y = 0.37,
                      label = paste0("Hedges g = ", round(age_effect_size, 2))),
        colour = "black", size = 20/.pt,
        hjust   = 0
        # vjust   = -43 # -28 would be lower
      )

    pdf(paste0("Figure3_Age_for", i, ".pdf"))
    print(p)
    dev.off()
  }

}

# cleaning
rm(extract_Age, extract_AB)
ls()

# # Fig3 AB / Age distribution in CM variants----
# extract_AB <- clonal_mosaicism %>%
#   filter(!is.na(considered_variants)) %>%
#   # filter(IDs %in% c("2-25234373-C-T", "2-25235726-A-G",
#   #                   "2-197402110-T-C")) %>%
#   distinct(IDs, .keep_all = TRUE)
# 
# gnomad <- extract_AB %>%
#   filter(IDs != "All individuals") %>%
#   select(-AB_count) %>%
#   # 1. Extract count per AB bin for heterozygous ind--
#   # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals;
#   # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00>
#   mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
#   separate(col = ab_hist_alt_bin_freq,
#            into = c("0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
#                     "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
#                     "0.40-0.45","0.45-0.50","0.50-0.55","0.55-0.60",
#                     "0.60-0.65","0.65-0.70","0.70-0.75","0.75-0.80",
#                     "0.80-0.85","0.85-0.90","0.90-0.95","0.95-1.00"),
#            sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
#   mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
#   add_row("IDs"="Reference", "ID"="DNMT3A 2-25247044-C-T",
#           "0.00-0.05" = 0, "0.05-0.10"  = 1, "0.10-0.15" = 22,
#           "0.15-0.20" = 53, "0.20-0.25" = 147, "0.25-0.30" = 545,
#           "0.30-0.35" = 1592, "0.35-0.40" = 3438, "0.40-0.45" = 8408,
#           "0.45-0.50" = 9816, "0.50-0.55" = 12487, "0.55-0.60" = 7145,
#           "0.60-0.65" = 3974, "0.65-0.70" = 1399, "0.70-0.75" = 419,
#           "0.75-0.80" = 208, "0.80-0.85" = 80, "0.85-0.90" = 84,
#           "0.90-0.95" = 92, "0.95-1.00" = 0, .before = 1
#   ) %>%
#   mutate(total_allele_balance = rowSums(select(.,`0.00-0.05`:`0.95-1.00`), na.rm = TRUE)) %>%
# 
#   select(IDs, everything())
# 
# gnomad <- gnomad %>%
#   pivot_longer(cols = c(grep("^[[:digit:]]", colnames(.))),
#                names_to = "allele_bin", values_to = "AB_count") %>%
#   mutate(allele_bin = factor(allele_bin,
#                              levels = c("0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
#                                         "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
#                                         "0.40-0.45","0.45-0.50","0.50-0.55","0.55-0.60",
#                                         "0.60-0.65","0.65-0.70","0.70-0.75","0.75-0.80",
#                                         "0.80-0.85","0.85-0.90","0.90-0.95","0.95-1.00"
#                              ))) %>%
#   mutate(AB_bin = case_when(
#     allele_bin == "0.00-0.05"            ~ 0.025,
#     allele_bin == "0.05-0.10"            ~ 0.075,
#     allele_bin == "0.10-0.15"            ~ 0.125,
#     allele_bin == "0.15-0.20"            ~ 0.175,
#     allele_bin == "0.20-0.25"            ~ 0.225,
#     allele_bin == "0.25-0.30"            ~ 0.275,
#     allele_bin == "0.30-0.35"            ~ 0.325,
#     allele_bin == "0.35-0.40"            ~ 0.375,
#     allele_bin == "0.40-0.45"            ~ 0.425,
#     allele_bin == "0.45-0.50"            ~ 0.475,
#     allele_bin == "0.50-0.55"            ~ 0.525,
#     allele_bin == "0.55-0.60"            ~ 0.575,
#     allele_bin == "0.60-0.65"            ~ 0.625,
#     allele_bin == "0.65-0.70"            ~ 0.675,
#     allele_bin == "0.70-0.75"            ~ 0.725,
#     allele_bin == "0.75-0.80"            ~ 0.775,
#     allele_bin == "0.80-0.85"            ~ 0.825,
#     allele_bin == "0.85-0.90"            ~ 0.875,
#     allele_bin == "0.90-0.95"            ~ 0.925,
#     allele_bin == "0.95-1.00"            ~ 0.975
#   )) %>%
#   mutate(AB_frequency = AB_count / total_allele_balance)
# 
# All <- gnomad %>%
#   filter(IDs == "Reference"
#   ) %>%
#   mutate(IDs = "Ref")
# 
# set.seed(1234)
# for(i in unique(gnomad$IDs)) {
# 
#   if (i != "Reference") {
#     # All <- gnomad %>%
#     #   filter(IDs == "Reference"
#     #   ) %>%
#     #   mutate(IDs = "Ref")
#     data_plot <- gnomad %>%
#       filter(IDs == i
#       ) %>%
#       bind_rows(., All) %>%
#       distinct(IDs, AB_bin, .keep_all = TRUE)
# 
#     plot_title <- data_plot %>%
#       filter(IDs == i) %>%
#       select(considered_variants) %>% distinct() %>% as.character()
#     # n <- n + 1
#     p <- data_plot %>%
#       ggplot(aes(x = AB_bin, y= AB_frequency, color= IDs))+
#       geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
#                   span = 0.6
#       )+
#       scale_color_manual(values = c("#F8766D", "#3333FF"))+
#       labs(title = plot_title)+
#       labs(x= "Allele balance", y= "Allele balance frequency")+
#       geom_col(data = data_plot %>%
#                  filter(IDs == i) ,
#                aes(x= AB_bin, y= AB_frequency),
#                position = "identity",
#                fill = "red",
#                alpha = 0.1)+
#       theme(axis.text = element_text(size = 24),
#             axis.title = element_text(size = 26),
#             plot.title = element_text(size = 24))+
#       ylim(0, 0.4)+
#       theme(legend.position = "none")
# 
#     # p <- p+
#     #   geom_text(
#     #     data = data_plot %>%
#     #       filter(IDs == i),
#     #     mapping = aes(x = 0.60, y = 0.4,
#     #                   label = paste0("q = ", round(adjusted_AB_pval, 3))),
#     #     colour = "black", size = 20/.pt,
#     #     hjust   = 0,
#     #     # vjust   = -28
#     #   )+
#     #   geom_text(
#     #     data = data_plot %>%
#     #       filter(IDs == i), #base_size = text.size,
#     #     mapping = aes(x = 0.60, y = 0.37,
#     #                   label = paste0("Hedges g = ", round(AB_effect_size, 3))),
#     #     colour = "black", size = 20/.pt,
#     #     hjust   = 0
#     #     # vjust   = -43 # -28 would be lower
#     #   )
# 
#     pdf(paste0("Figure3_CM_AB_for", i, ".pdf"))
#     print(p)
#     dev.off()
# 
#   }
# }
# 
# extract_Age <- clonal_mosaicism %>%
#   filter(!is.na(considered_variants)) %>%
#   # filter(IDs %in% c("2-25234373-C-T", "2-25235726-A-G",
#   #                   "2-197402110-T-C")) %>%
#   distinct(IDs, .keep_all = TRUE)
# 
# extract_Age <- extract_Age %>%
# 
#   select(-sample_frequency, -sample_count, -nbr_individuals) %>%
# 
#   # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
#   mutate("<30" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
#   mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
#   separate(col = age_hist_het_bin_freq,
#            into = c("30-35","35-40","40-45","45-50",
#                     "50-55","55-60","60-65","65-70",
#                     "70-75","75-80"),
#            sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
#   # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
#   mutate(">80" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>%
#   mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
# 
#   add_row("IDs"="All individuals", "ID"="All individuals of any genotype bin",
#           "<30" = 2547, "30-35" = 3423,
#           "35-40" = 4546,"40-45" = 8487,
#           "45-50" = 10355,"50-55" = 12693,
#           "55-60" = 11933,"60-65" = 10534,
#           "65-70" = 8882,"70-75" = 5991,
#           "75-80" = 4136,
#           ">80" = 1935, .before = 1
#   ) %>%
#   mutate(nbr_individuals = rowSums(select(.,`<30`:`>80`), na.rm = TRUE))
# 
# extract_Age <- extract_Age %>%
#   pivot_longer(cols = c(grep("^([[:digit:]]|<|>)", colnames(.))),
#                names_to = "age_bin", values_to = "sample_count") %>%
#   mutate(age_bin = factor(age_bin,
#                           levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
#                                      "55-60","60-65","65-70","70-75","75-80", ">80"
#                           ))) %>%
#   mutate(bin_age = case_when(
#     # The average age at diagnosis is 8 overall (ages 0 to 19),
#     # 5 years old for children (aged 0 to 14), and 17 years old
#     # for adolescents (aged 15 to 19), while adults' average
#     # age for cancer diagnosis is 65
#     str_detect(age_bin, "<30")              ~ 15,
#     str_detect(age_bin, "30-35")            ~ 32.5,
#     str_detect(age_bin, "35-40")            ~ 37.5,
#     str_detect(age_bin, "40-45")            ~ 42.5,
#     str_detect(age_bin, "45-50")            ~ 47.5,
#     str_detect(age_bin, "50-55")            ~ 52.5,
#     str_detect(age_bin, "55-60")            ~ 57.5,
#     str_detect(age_bin, "60-65")            ~ 62.5,
#     str_detect(age_bin, "65-70")            ~ 67.5,
#     str_detect(age_bin, "70-75")            ~ 72.5,
#     str_detect(age_bin, "75-80")            ~ 77.5,
#     TRUE ~ 90
#   )) %>%
#   mutate(age_frequency = sample_count / nbr_individuals)
# 
# All <- extract_Age %>%
#   filter(IDs == "All individuals"
#   )
# 
# set.seed(1234)
# 
# for(i in unique(extract_Age$IDs)) {
# 
#   if (i != "All individuals") {
#     # All <- extract_Age %>%
#     #   filter(IDs == "All individuals"
#     # ) %>%
#     #   mutate(IDs = paste0("Reference= ", ID))
#     data_plot <- extract_Age %>%
#       filter(IDs == i
#       ) %>%
#       bind_rows(., All) %>%
#       distinct(IDs, bin_age, .keep_all = TRUE)
# 
#     plot_title <- data_plot %>%
#       filter(IDs == i) %>%
#       select(considered_variants) %>% distinct() %>% as.character()
#     # n <- n + 1
#     p <- data_plot %>%
#       ggplot(aes(x = bin_age, y= age_frequency, color= IDs))+
#       geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
#                   span = 0.6
#       )+
#       scale_color_manual(values = c("#F8766D", "#3333FF"))+
#       labs(title = plot_title)+
#       labs(x= "Age", y= "Age frequency")+
#       geom_col(data = data_plot %>%
#                  filter(IDs == i) ,
#                aes(x= bin_age, y= age_frequency),
#                position = "identity",
#                fill = "red",
#                alpha = 0.1)+
#       theme(axis.text = element_text(size = 24),
#             axis.title = element_text(size = 26),
#             plot.title = element_text(size = 24))+
#       ylim(0, max(data_plot$age_frequency))+
#       theme(legend.position = "none")
# 
#     # p <- p+
#     #   geom_text(
#     #     data    = data_plot %>%
#     #       filter(IDs == i),
#     #     mapping = aes(x = 25, y = -Inf,
#     #                   label = paste0("BH= ", round(adjusted_age_pval, 3))),
#     #     colour = "black",
#     #     # hjust   = 0,
#     #     vjust   = -28
#     #   )
# 
#     pdf(paste0("Figure3_CM_Age_for", i, ".pdf"))
#     print(p)
#     dev.off()
#   }
# 
# }
# 
# # cleaning
# rm(extract_Age, extract_AB)
# ls()
# 
# # Fig3_A Table Chromosomes----
# print("Table Chromosome")
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(X.CHROM) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency")) %>% 
#   bold_labels() %>% 
#   modify_header(
#     label = "**Number of CH Variants per Chromosomes**"
#   ) %>% as_kable()
# 
# a <- CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(X.CHROM, SYMBOL, exon_effective_length) %>% 
#   group_by(X.CHROM, SYMBOL) %>% 
#   mutate(exon_effective_length_sum1 = row_number(SYMBOL)) %>% 
#   mutate(exon_effective_length_sum2 = case_when(
#     exon_effective_length_sum1 == 1       ~ exon_effective_length
#   )) %>% 
#   group_by(X.CHROM) %>% 
#   mutate(exon_effective_length_sum = sum(exon_effective_length_sum2, na.rm = TRUE)) %>%
#   select(-c(exon_effective_length_sum1, exon_effective_length_sum2)) %>% 
#   group_by(X.CHROM) %>% 
#   mutate(variants_nbr = n()) %>% 
#   ungroup() %>% 
#   mutate(results = variants_nbr / exon_effective_length_sum) %>% 
#   distinct(X.CHROM, .keep_all = TRUE) %>% 
#   select(-exon_effective_length, -SYMBOL) %>% 
#   mutate(total = sum(results)) %>% 
#   mutate(results_percent = (results / total) * 100) %>% 
#   arrange(desc(results))
# print(a, n = 30)
# clonal_mosaicism %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(X.CHROM) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency")) %>% 
#   bold_labels() %>% 
#   modify_header(
#     label = "**Number of Clonal Mosaicism Variants per Chromosomes**"
#   ) %>% as_kable()
# a <- clonal_mosaicism %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(X.CHROM, SYMBOL, exon_effective_length) %>% 
#   group_by(X.CHROM, SYMBOL) %>% 
#   mutate(exon_effective_length_sum1 = row_number(SYMBOL)) %>% 
#   mutate(exon_effective_length_sum2 = case_when(
#     exon_effective_length_sum1 == 1       ~ exon_effective_length
#   )) %>% 
#   group_by(X.CHROM) %>% 
#   mutate(exon_effective_length_sum = sum(exon_effective_length_sum2, na.rm = TRUE)) %>%
#   select(-c(exon_effective_length_sum1, exon_effective_length_sum2)) %>% 
#   group_by(X.CHROM) %>% 
#   mutate(variants_nbr = n()) %>% 
#   ungroup() %>% 
#   mutate(results = variants_nbr / exon_effective_length_sum) %>% 
#   distinct(X.CHROM, .keep_all = TRUE) %>% 
#   select(-exon_effective_length, -SYMBOL) %>% 
#   mutate(total = sum(results)) %>% 
#   mutate(results_percent = (results / total) * 100) %>% 
#   arrange(desc(results))
# print(a, n = 30)

# print("FIGURE 3_B")
# Fig3_B Compare top 10-20 genes for CH vs CM----
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   filter(SYMBOL =="DNMT3A") %>% 
#   mutate(AN_number = str_match(INFO, "AN=(.*?);")[,2]) %>% 
#   select(IDs, SYMBOL, AN_number, alt_allele_count) %>% 
#   mutate(AN_divided_by_2 = (as.numeric(AN_number) / 2)
#   ) %>% group_by(SYMBOL) %>% 
#   mutate(AN_divided_by_2_sum = sum(AN_divided_by_2)) %>% head()

# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(IDs, SYMBOL) %>% 
#   mutate(data = "CH variants") %>% 
#   bind_rows(., clonal_mosaicism %>% 
#               distinct(IDs, .keep_all = TRUE) %>% 
#               select(IDs, SYMBOL) %>% 
#               mutate(data = "CM variants")) %>% 
#   select(SYMBOL, data) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency")) %>% 
#   bold_labels() %>% 
#   modify_header(
#     label = "**Compare top 10-20 genes for CH vs CM - not prevalence**"
#   ) %>% as_kable()

print("AC per variants in CH")
write_csv(CH_variants %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, nc_alt_allele_count, mutations_in_BickWHO, L_CHIP,
                   alt_allele_count, alt_allele_count_female, alt_allele_count_male,
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_freq_allele_count, nc_freq_allele_count_female, nc_freq_allele_count_male,
                   nc_freq_allele_count_afr,
                   nc_freq_allele_count_amr,
                   nc_freq_allele_count_eas,
                   nc_freq_allele_count_sas,
                   nc_freq_allele_count_nfe),
          "Overall AC in CH.csv")
write_csv(CH_variants %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC over gene in CH.csv")
write_csv(CH_variants %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(mutations_in_BickWHO == "Yes") %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC BW region over gene in CH.csv")
write_csv(CH_variants %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(L_CHIP == "Yes") %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC L-CH region over gene in CH.csv")
write_csv(CH_variants %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, nc_alt_allele_count, mutations_in_BickWHO, L_CHIP,
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe,
                   exon_effective_length) %>% 
            group_by(SYMBOL, exon_effective_length) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            mutate(across(where(is.numeric), ~ . / exon_effective_length)) %>% 
            arrange(desc(nc_alt_allele_count)),
          "Relative AC over exon length in CH.csv")


write_csv(clonal_mosaicism %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, nc_alt_allele_count, mutations_in_BickWHO, L_CHIP,
                   alt_allele_count, alt_allele_count_female, alt_allele_count_male,
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_freq_allele_count, nc_freq_allele_count_female, nc_freq_allele_count_male,
                   nc_freq_allele_count_afr,
                   nc_freq_allele_count_amr,
                   nc_freq_allele_count_eas,
                   nc_freq_allele_count_sas,
                   nc_freq_allele_count_nfe),
          "Overall AC in CM.csv")
write_csv(clonal_mosaicism %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC over gene in SM.csv")
write_csv(clonal_mosaicism %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(mutations_in_BickWHO == "Yes") %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC BW region over gene in SM.csv")
write_csv(clonal_mosaicism %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            filter(L_CHIP == "Yes") %>% 
            select(SYMBOL, 
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe) %>% 
            group_by(SYMBOL) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            arrange(desc(nc_alt_allele_count)),
          "summarized AC L-CH region over gene in SM.csv")
write_csv(clonal_mosaicism %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(SYMBOL, nc_alt_allele_count, mutations_in_BickWHO, L_CHIP,
                   alt_allele_count, 
                   alt_allele_count_afr,
                   alt_allele_count_amr,
                   alt_allele_count_eas,
                   alt_allele_count_sas,
                   alt_allele_count_nfe,
                   nc_alt_allele_count,
                   nc_alt_allele_count_afr,
                   nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas,
                   nc_alt_allele_count_sas,
                   nc_alt_allele_count_nfe,
                   exon_effective_length) %>% 
            group_by(SYMBOL, exon_effective_length) %>% 
            mutate(across(where(is.numeric), ~ sum(.))) %>%
            distinct(SYMBOL, .keep_all = TRUE) %>% 
            ungroup() %>% 
            mutate(across(where(is.numeric), ~ . / exon_effective_length)) %>% 
            arrange(desc(nc_alt_allele_count)),
          "Relative AC over exon length in CM.csv")
print("The first CH 30 genes with higer cumulative alt count")
a <- CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)
print("The first CH 30 genes with higer cumulative non cancer alt count")
a <- CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(nc_alt_allele_count_sum = sum(nc_alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_sum)) %>% head(30)
print(a, n=30)

print("The first SM 30 genes with higer cumulative alt count")
a <- clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)
print("The first SM 30 genes with higer cumulative non cancer alt count")
a <- clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(nc_alt_allele_count_sum = sum(nc_alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_sum)) %>% head(30)
print(a, n=30)

print("The first CH 30 genes-region in BW with higer cumulative alt count")
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)

print("The first SM 30 genes-region in BW with higer cumulative alt count")
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)

print("The first CH 30 full genes-not region in BW with higer cumulative alt count")
a <- CH_variants %>% 
  group_by(SYMBOL) %>% 
  fill(mutations_in_BickWHO, .direction = "updown") %>% 
  ungroup() %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)

print("The first SM 30 full genes-not region in BW with higer cumulative alt count")
a <- clonal_mosaicism %>% 
  group_by(SYMBOL) %>% 
  fill(mutations_in_BickWHO, .direction = "updown") %>% 
  ungroup() %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(alt_allele_count_sum = sum(alt_allele_count)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sum)) %>% head(30)
print(a, n=30)




# Fig2_C Table Yi-Han----
print("FIGURE 3_A")
print("CH VEPs extracted from Yi-Han")
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(Func.refGene, ExonicFunc.refGene, 
         new_HGVSc, SIFT_pred, Polyphen2_HDIV_pred, 
         FATHMM_pred, CADD_phred, CADD_phred_cat, CLNSIG,
         REF_dbGaP_PopFreq, ALT_dbGaP_PopFreq) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**CH VEPs extracted from Yi-Han**"
  ) %>% as_kable()

CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(REF_dbGaP_PopFreq, ALT_dbGaP_PopFreq) %>% 
  tbl_summary(sort = list(everything() ~ "frequency"),
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "** CH dbGaP ranges extracted from Yi-Han**"
  ) %>% as_kable()

# pdf("dbgap boxplot.pdf")
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(REF_dbGaP_PopFreq) %>% 
#   mutate(name = "REF_dbGaP_PopFreq") %>% 
#   # pivot_longer(cols = everything()) %>% 
#   ggplot(aes(x= name, y= REF_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()
# 
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(REF_dbGaP_PopFreq) %>% 
#   mutate(name = "REF_dbGaP_PopFreq") %>% 
#   # pivot_longer(cols = everything()) %>% 
#   ggplot(aes(x= name, y= REF_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()+
#   ylim(0.9965, 1)
# 
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(REF_dbGaP_PopFreq) %>% 
#   mutate(name = "REF_dbGaP_PopFreq") %>% 
#   # pivot_longer(cols = everything()) %>% 
#   ggplot(aes(x= name, y= REF_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()+
#   facet_zoom(ylim = c(0.999, 1))
# 
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(ALT_dbGaP_PopFreq) %>% 
#   mutate(name = "ALT_dbGaP_PopFreq") %>% 
#   ggplot(aes(x= name, y= ALT_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()
# 
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(ALT_dbGaP_PopFreq) %>% 
#   mutate(name = "ALT_dbGaP_PopFreq") %>% 
#   ggplot(aes(x= name, y= ALT_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()+
#   ylim(0, 0.003)
# 
# CH_variants %>% 
#   distinct(IDs, .keep_all = TRUE) %>% 
#   select(ALT_dbGaP_PopFreq) %>% 
#   mutate(name = "ALT_dbGaP_PopFreq") %>% 
#   ggplot(aes(x= name, y= ALT_dbGaP_PopFreq))+
#   geom_boxplot()+
#   geom_jitter()+
#   facet_zoom(ylim = c(0, 0.001))
# dev.off()

pdf("dbgap_distribution.pdf")
clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(ALT_dbGaP_PopFreq) %>%
  mutate(ALT_dbGaP_PopFreq = round(ALT_dbGaP_PopFreq, 2)) %>% 
  group_by(ALT_dbGaP_PopFreq) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  # mutate(name = "ALT_dbGaP_PopFreq") %>%
  ggplot(aes(x= ALT_dbGaP_PopFreq, y = n))+
  geom_area(fill="blue")+
  xlim(0, 0.25)
dev.off()


vep_in_gnomAD <- CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2])
vep_in_gnomAD <- vep_in_gnomAD %>%
  separate(col = ens_vep, paste("ens_vep", 1:25, sep=""),
           sep = ",", remove = T, extra = "warn", fill = "right") %>%
  keep(~!all(is.na(.))) %>%
  pivot_longer(cols = starts_with("ens_vep"), names_to = NULL, values_to = "ens_vep") %>%
  drop_na(ens_vep)

ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
                       "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
                       "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
                       "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
                       "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
                       "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
                       "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
                       "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
                       "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
                       "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
                       "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
                       "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
                       "LoF_flags", "LoF_info")

vep_in_gnomAD <- vep_in_gnomAD %>%
  separate(col = ens_vep, into = ens_vep_var_names,
           sep = "\\|", remove = F, extra = "warn", fill = "right")

vep_in_gnomAD <- vep_in_gnomAD %>%
  mutate(SOMATIC_1 = case_when(
    str_detect(SOMATIC, "1")      ~ 1,
    TRUE                          ~ 0
  )) %>% 
  group_by(IDs) %>% 
  mutate(SOMATIC = max(SOMATIC_1),
         SOMATIC = case_when(
           SOMATIC == 1                  ~ "Yes",
           TRUE                          ~ "No"
         )) %>% 
  ungroup() %>% 
  distinct(IDs, .keep_all = TRUE)

vep_in_gnomAD %>% 
  select(
    SOMATIC) %>% 
  tbl_summary(type = list(c(SOMATIC) ~ "categorical"),
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**VEP in CH variants**"
  ) %>% as_kable()


print("CM VEPs extracted from Yi-Han")
clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(Func.refGene, ExonicFunc.refGene, 
         new_HGVSc, SIFT_pred, Polyphen2_HDIV_pred, 
         FATHMM_pred, CADD_phred, CADD_phred_cat, CLNSIG,
         REF_dbGaP_PopFreq, ALT_dbGaP_PopFreq) %>% 
  tbl_summary(sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**CM VEPs extracted from Yi-Han**"
  ) %>% as_kable()

clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(REF_dbGaP_PopFreq, ALT_dbGaP_PopFreq) %>% 
  tbl_summary(sort = list(everything() ~ "frequency"),
              statistic = list(all_continuous() ~ "{median} ({min}, {max})")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**CM dbGaP ranges extracted from Yi-Han**"
  ) %>% as_kable()

vep_in_gnomAD <- clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2])
vep_in_gnomAD <- vep_in_gnomAD %>%
  separate(col = ens_vep, paste("ens_vep", 1:25, sep=""),
           sep = ",", remove = T, extra = "warn", fill = "right") %>%
  keep(~!all(is.na(.))) %>%
  pivot_longer(cols = starts_with("ens_vep"), names_to = NULL, values_to = "ens_vep") %>%
  drop_na(ens_vep)

ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
                       "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
                       "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
                       "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
                       "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
                       "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
                       "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
                       "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
                       "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
                       "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
                       "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
                       "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
                       "LoF_flags", "LoF_info")

vep_in_gnomAD <- vep_in_gnomAD %>%
  separate(col = ens_vep, into = ens_vep_var_names,
           sep = "\\|", remove = F, extra = "warn", fill = "right")

vep_in_gnomAD <- vep_in_gnomAD %>%
  mutate(SOMATIC_1 = case_when(
    str_detect(SOMATIC, "1")      ~ 1,
    TRUE                          ~ 0
  )) %>% 
  group_by(IDs) %>% 
  mutate(SOMATIC = max(SOMATIC_1),
         SOMATIC = case_when(
           SOMATIC == 1                  ~ "Yes",
           TRUE                          ~ "No"
         )) %>% 
  ungroup() %>% 
  distinct(IDs, .keep_all = TRUE)

vep_in_gnomAD %>% 
  select(
    SOMATIC) %>% 
  tbl_summary(type = list(c(SOMATIC) ~ "categorical"),
              sort = list(everything() ~ "frequency")) %>% 
  bold_labels() %>% 
  modify_header(
    label = "**VEP in Clonal Mosaicism Variants**"
  ) %>% as_kable()

# cleaning
rm(ens_vep_var_names, vep_in_gnomAD)




# Fig4----
print("FIGURE 4")

write_rds(CH_variants %>%
        distinct(IDs, .keep_all = TRUE) %>%
        unite(IDs, c(IDs, SYMBOL)) %>%
        select(IDs,
               nc_freq_allele_count_afr,
               nc_freq_allele_count_amr,
               nc_freq_allele_count_eas,
               nc_freq_allele_count_sas,
               nc_freq_allele_count_nfe
        ), "CH_variants_PAF_data.rds")

write_rds(clonal_mosaicism %>%
        distinct(IDs, .keep_all = TRUE) %>%
        unite(IDs, c(IDs, SYMBOL)) %>%
        select(IDs,
               nc_freq_allele_count_afr,
               nc_freq_allele_count_amr,
               nc_freq_allele_count_eas,
               nc_freq_allele_count_sas,
               nc_freq_allele_count_nfe
        ), "SM_variants_PAF_data.rds")

print("Fig4F_SM") 
# Fig4E_CH ----
pdf("Fig4F_CH_Popultaion_AF_Distribution.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  unite(IDs, c(IDs, SYMBOL)) %>%
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas",
           "sas_male", "sas_female", "sas",
           "amr_male", "amr_female", "amr",
           "afr_male", "afr_female", "afr")
         )) %>%
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>% 
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color=Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in CH Variants")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()


CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})")
              # digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF race/eth**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              # digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF race/eth**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              # digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF race/eth without the 0 PAF - mean**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              # digits = all_continuous() ~ 20,
              # statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF race/eth without the 0 PAF - median**"
  ) %>% as_kable()

print("Pop freq raw value for nfe")
mean(CH_variants$nc_freq_allele_count_nfe, na.rm= TRUE)
median(CH_variants$nc_freq_allele_count_nfe, na.rm= TRUE)
range(CH_variants$nc_freq_allele_count_nfe, na.rm= TRUE)
sd(CH_variants$nc_freq_allele_count_nfe, na.rm= TRUE)
quantile(CH_variants$nc_freq_allele_count_nfe, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop freq raw value for afr")
mean(CH_variants$nc_freq_allele_count_afr, na.rm= TRUE)
median(CH_variants$nc_freq_allele_count_afr, na.rm= TRUE)
range(CH_variants$nc_freq_allele_count_afr, na.rm= TRUE)
sd(CH_variants$nc_freq_allele_count_afr, na.rm= TRUE)
quantile(CH_variants$nc_freq_allele_count_afr, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop freq raw value for eas")
mean(CH_variants$nc_freq_allele_count_eas, na.rm= TRUE)
median(CH_variants$nc_freq_allele_count_eas, na.rm= TRUE)
range(CH_variants$nc_freq_allele_count_eas, na.rm= TRUE)
sd(CH_variants$nc_freq_allele_count_eas, na.rm= TRUE)
quantile(CH_variants$nc_freq_allele_count_eas, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop freq raw value for sas")
mean(CH_variants$nc_freq_allele_count_sas, na.rm= TRUE)
median(CH_variants$nc_freq_allele_count_sas, na.rm= TRUE)
range(CH_variants$nc_freq_allele_count_sas, na.rm= TRUE)
sd(CH_variants$nc_freq_allele_count_sas, na.rm= TRUE)
quantile(CH_variants$nc_freq_allele_count_sas, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop freq raw value for amr")
mean(CH_variants$nc_freq_allele_count_amr, na.rm= TRUE)
median(CH_variants$nc_freq_allele_count_amr, na.rm= TRUE)
range(CH_variants$nc_freq_allele_count_amr, na.rm= TRUE)
sd(CH_variants$nc_freq_allele_count_amr, na.rm= TRUE)
quantile(CH_variants$nc_freq_allele_count_amr, probs = c(.25, .5, 0.75), na.rm= TRUE)

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})")
              # digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF male vs female**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Pop AF male vs female**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CH Pop AF male vs female without the 0 PAF - mean**"
  ) %>% as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              # statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CH Pop AF male vs female without the 0 PAF - median**"
  ) %>% as_kable()

print("Fig4E_SM") 
# Fig4E_SM ----
pdf("Fig4E_SM_Popultaion_AF_Distribution.pdf")
clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  unite(IDs, c(IDs, SYMBOL)) %>%
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas",
           "sas_male", "sas_female", "sas",
           "amr_male", "amr_female", "amr",
           "afr_male", "afr_female", "afr")
         )) %>%
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>% 
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color=Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in SM Variants")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})")
              # digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CM Pop AF race/eth**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CM Pop AF race/eth**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              # digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF race/eth without the 0 PAF - mean**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_afr, nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
         nc_freq_allele_count_amr, nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
         nc_freq_allele_count_eas, nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
         nc_freq_allele_count_sas, nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
         nc_freq_allele_count_nfe, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         # name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           # "male", "female", "overall", 
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>% 
  filter(Ethnicity != "overall") %>% 
  filter(Sex == "overall") %>% 
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              # digits = all_continuous() ~ 20,
              # statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF race/eth without the 0 PAF - median**"
  ) %>% as_kable()

print("Pop SM freq raw value for nfe")
mean(clonal_mosaicism$nc_freq_allele_count_nfe)
median(clonal_mosaicism$nc_freq_allele_count_nfe)
range(clonal_mosaicism$nc_freq_allele_count_nfe)
sd(clonal_mosaicism$nc_freq_allele_count_nfe)
quantile(clonal_mosaicism$nc_freq_allele_count_nfe, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop SM freq raw value for afr")
mean(clonal_mosaicism$nc_freq_allele_count_afr, na.rm= TRUE)
median(clonal_mosaicism$nc_freq_allele_count_afr, na.rm= TRUE)
range(clonal_mosaicism$nc_freq_allele_count_afr, na.rm= TRUE)
sd(clonal_mosaicism$nc_freq_allele_count_afr, na.rm= TRUE)
quantile(clonal_mosaicism$nc_freq_allele_count_afr, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop SM freq raw value for eas")
mean(clonal_mosaicism$nc_freq_allele_count_eas, na.rm= TRUE)
median(clonal_mosaicism$nc_freq_allele_count_eas, na.rm= TRUE)
range(clonal_mosaicism$nc_freq_allele_count_eas, na.rm= TRUE)
sd(clonal_mosaicism$nc_freq_allele_count_eas, na.rm= TRUE)
quantile(clonal_mosaicism$nc_freq_allele_count_eas, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop SM freq raw value for sas")
mean(clonal_mosaicism$nc_freq_allele_count_sas, na.rm= TRUE)
median(clonal_mosaicism$nc_freq_allele_count_sas, na.rm= TRUE)
range(clonal_mosaicism$nc_freq_allele_count_sas, na.rm= TRUE)
sd(clonal_mosaicism$nc_freq_allele_count_sas, na.rm= TRUE)
quantile(clonal_mosaicism$nc_freq_allele_count_sas, probs = c(.25, .5, 0.75), na.rm= TRUE)

print("Pop SM freq raw value for amr")
mean(clonal_mosaicism$nc_freq_allele_count_amr, na.rm= TRUE)
median(clonal_mosaicism$nc_freq_allele_count_amr, na.rm= TRUE)
range(clonal_mosaicism$nc_freq_allele_count_amr, na.rm= TRUE)
sd(clonal_mosaicism$nc_freq_allele_count_amr, na.rm= TRUE)
quantile(clonal_mosaicism$nc_freq_allele_count_amr, probs = c(.25, .5, 0.75), na.rm= TRUE)

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              digits = all_continuous() ~ 20,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})")
              # digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CM Pop AF male vs female**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CM Pop AF male vs female**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF male vs female without the 0 PAF - mean**"
  ) %>% as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs,
         nc_freq_allele_count_female, nc_freq_allele_count_male
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = str_replace(name, "nc_freq_allele_count", "overall"),
         name = factor(name, levels = c(
           "male", "female")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  filter(value != 0) %>% 
  select(`Population Allele Frequency` = value, Sex) %>% 
  tbl_summary(by= Sex,
              # statistic = list(all_continuous() ~ "{mean} ({p25}, {p75})"),
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF male vs female without the 0 PAF - median**"
  ) %>% as_kable()

# Supplemental fig Pop AF----
pdf("Supplemental_figure_S6_Population_AF_per_IDs.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(considered_variants)) %>%
  filter(considered_variants != "DNMT3A W860R") %>%
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs, considered_variants,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -c(IDs, considered_variants)) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas",
           "sas_male", "sas_female", "sas",
           "amr_male", "amr_female", "amr",
           "afr_male", "afr_female", "afr")
         )) %>%
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex, considered_variants) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>% 
  
  ggplot(aes(x=Ethnicity, y= value, color= Ethnicity))+
  geom_point(size= 3)+
  labs(#title = "Variant Allele Frequency in Non-Cancer population", 
       y= "Popultaion Allele Frequency", x= "Frequency Density by Population")+
         
  scale_color_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  theme_classic(base_size = 16)+
  theme(legend.title=element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust=1),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "italic"))+
  facet_wrap(. ~ considered_variants, ncol = 2)+
  coord_flip()
dev.off()

pdf("Supplemental_CM_figure_6F_Population_AF_per_IDs.pdf")
clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(considered_variants)) %>%
  filter(considered_variants != "DNMT3A W860R") %>%
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs, considered_variants,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -c(IDs, considered_variants)) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas",
           "sas_male", "sas_female", "sas",
           "amr_male", "amr_female", "amr",
           "afr_male", "afr_female", "afr")
         )) %>%
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex, considered_variants) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>% 
  
  ggplot(aes(x=Ethnicity, y= value, color= Ethnicity))+
  geom_point(size= 3)+
  labs(#title = "Variant Allele Frequency in Non-Cancer population", 
       y= "Popultaion Allele Frequency", x= "Frequency Density by Population")+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  theme_classic(base_size = 16)+
  theme(legend.title=element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust=1),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "italic"))+
  facet_wrap(. ~ considered_variants, ncol = 2)+
  coord_flip()
dev.off()

# Figure rates of CH by age----
pdf("Fig3C_Rates_of_CH_by_age.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%

  select(-sample_frequency, -sample_count, -nbr_individuals) %>%

  # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
  mutate("<30" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
  mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
  separate(col = age_hist_het_bin_freq,
           into = c("30-35","35-40","40-45","45-50",
                    "50-55","55-60","60-65","65-70",
                    "70-75","75-80"),
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
  mutate(">80" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>%
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%

  select(IDs, SYMBOL, "<30" : ">80") %>%
  summarise_at(c("<30",
                 "30-35","35-40","40-45","45-50",
                 "50-55","55-60","60-65","65-70",
                 "70-75","75-80", ">80"), ~ sum(.)) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name, levels = c("<30",
                                        "30-35","35-40","40-45","45-50",
                                        "50-55","55-60","60-65","65-70",
                                        "70-75","75-80", ">80"))) %>%
  mutate(sum = sum(value),
         perc = (value / sum) * 100) %>% 
  ggplot(aes(x= name, y= perc))+
  geom_bar(stat = "identity", fill= "darkblue")+
  labs(#title = "CH variants", 
       x= "Age",
       y= "Individuals (#)")+
  ylim(0, 20.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf("Fig3C_Rates_of_CM_by_age.pdf")
clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%

  select(-sample_frequency, -sample_count, -nbr_individuals) %>%

  # "Count of age values falling below lowest histogram bin edge for heterozygous individuals">
  mutate("<30" = as.numeric(str_match(INFO, "age_hist_het_n_smaller=(.*?);")[,2])) %>%
  mutate(age_hist_het_bin_freq = str_match(INFO, "age_hist_het_bin_freq=(.*?);")[,2]) %>%
  separate(col = age_hist_het_bin_freq,
           into = c("30-35","35-40","40-45","45-50",
                    "50-55","55-60","60-65","65-70",
                    "70-75","75-80"),
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  # "Count of age values falling above highest histogram bin edge for heterozygous individuals">
  mutate(">80" = as.numeric(str_match(INFO, "age_hist_het_n_larger=(.*?);")[,2])) %>%
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%

  select(IDs, SYMBOL, "<30" : ">80") %>%
  summarise_at(c("<30",
                 "30-35","35-40","40-45","45-50",
                 "50-55","55-60","60-65","65-70",
                 "70-75","75-80", ">80"), ~ sum(.)) %>%
  pivot_longer(cols = everything()) %>%
  mutate(name = factor(name, levels = c("<30",
                                        "30-35","35-40","40-45","45-50",
                                        "50-55","55-60","60-65","65-70",
                                        "70-75","75-80", ">80"))) %>%
  mutate(sum = sum(value),
         perc = (value / sum) * 100) %>% 
  ggplot(aes(x= name, y= perc))+
  geom_bar(stat = "identity", fill= "tomato")+
  labs(#title = "CM variants", 
       x= "Age",
       y= "Individuals (#)")+
  ylim(0, 20.1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


# Fig5A Age in COSMIC----
a <- CH_variants %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Age, Sample.name, Primary.site, mutations_in_BickWHO) %>%
  mutate(data = "CH variants")

pdf("FigS5_Age_in_COSMIC.pdf")
gnomad_variants %>%
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>%
  select(Age) %>%
  mutate(data = "no CH variants") %>%
  bind_rows(., a
  ) %>%
  group_by(data) %>%
  mutate(median_ = median(Age, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(data = factor(data, levels = c("no CH variants",
                                        "CH variants"))) %>% 
  ggplot(aes(x=Age, fill= data))+
  geom_histogram(binwidth = 1)+
  scale_fill_manual(values = c("lightgrey", "purple"), name="")+
  geom_vline(aes(xintercept= median_, color= data))+
  scale_color_manual(values = c("darkgrey", "purple"))+
  labs(y= "Individuals (#)")+
  guides(color = FALSE)
dev.off()

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(variant_in_cosmic == "Yes") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not CH variants") %>% 
  bind_rows(., a
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels()%>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>% as_kable()

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not CH") %>% 
  bind_rows(., a %>% 
              filter(Primary.site == "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic HEME Patients Characteristics**"
  ) %>% as_kable()

gnomad_variants %>% 
  # filter(variant_in_cosmic == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not CH") %>% 
  bind_rows(., a %>% 
              filter(Primary.site != "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic non-HEME Patients Characteristics**"
  ) %>% as_kable()









a <- clonal_mosaicism %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Age, Sample.name, Primary.site, mutations_in_BickWHO) %>%
  mutate(data = "SM variants")

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not SM variants") %>% 
  bind_rows(., a
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels()%>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>% as_kable()

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not SM") %>% 
  bind_rows(., a %>% 
              filter(Primary.site == "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic HEME Patients Characteristics**"
  ) %>% as_kable()

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not SM") %>% 
  bind_rows(., a %>% 
              filter(Primary.site != "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic non-HEME Patients Characteristics**"
  ) %>% as_kable()









# COSMIC Age in bick
a <- CH_variants %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  arrange(Primary.site) %>%
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Age, Sample.name, Primary.site, mutations_in_BickWHO) %>%
  mutate(data = "CH variants")

gnomad_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not CH") %>% 
  bind_rows(., a %>% 
              filter(mutations_in_BickWHO == "Yes") %>%
              filter(Primary.site == "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic HEME Patients Characteristics in Bick**"
  ) %>% as_kable()

gnomad_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Age) %>%
  mutate(data = "not CH") %>% 
  bind_rows(., a %>% 
              filter(mutations_in_BickWHO == "Yes") %>% 
              filter(Primary.site != "haematopoietic_and_lymphoid_tissue")
  ) %>% 
  select(Age, data) %>% 
  tbl_summary(by=data) %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cosmic non-HEME Patients Characteristics in Bick**"
  ) %>% as_kable()

a %>% 
  bind_rows(., gnomad_variants %>% 
              filter(variant_in_cosmic == "Yes") %>% 
              arrange(Primary.site) %>%
              distinct(Sample.name, .keep_all = TRUE) %>%
              # filter(!Sample.name %in% c(a$Sample.name)) %>% 
              mutate(data = "not CH")) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Primary.site, data) %>% 
  tbl_summary(by=data,
              sort = everything() ~ "frequency", percent = "row") %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**TRUE Cancer by row percent**"
  ) %>%
  as_kable()



gnomad_variants %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  mutate(data = "not CH") %>% 
  bind_rows(., a) %>% 
  select(Primary.site, data) %>% 
  tbl_summary(by=data,
              sort = everything() ~ "frequency", percent = "row") %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = "**Cancer by row percent**"
  ) %>%
  as_kable()


# Fig5B Cancers in COSMIC----
a <- CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Primary.site, Primary.histology) %>% 
  mutate(Primary.histology = case_when(
    Primary.site == "haematopoietic_and_lymphoid_tissue"    ~ Primary.histology,
    TRUE    ~ NA_character_
  )) %>% 
  group_by(Primary.site, Primary.histology) %>%
  summarise(sum_individuals_by_primary = n()) %>%
  ungroup()
write_rds(a, "Fig5B_Prevalence_Cancers_inCosmic.rds")

print("FIGURE 5")
# pdf("Fig5B_Prevalence_Cancers_inCosmic.pdf"#, width = 7, height = 5
# )
# CH_variants %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>%
#   select(Primary.site, Primary.histology) %>% 
#   mutate(Primary.histology = case_when(
#     Primary.site == "haematopoietic_and_lymphoid_tissue"    ~ Primary.histology,
#     TRUE    ~ NA_character_
#   ), Primary.histology = factor(Primary.histology, 
#                                 levels = c(
#                                   "lymphoid_neoplasm",
#                                   "haematopoietic_neoplasm",
#                                   "other"))) %>% 
#   group_by(Primary.site, Primary.histology) %>%
#   summarise(sum_individuals_by_primary = n()) %>%
#   ungroup() %>% 
#   ggplot(aes(x= fct_reorder(Primary.site, sum_individuals_by_primary), 
#              y=sum_individuals_by_primary, fill= Primary.histology))+
#   geom_bar(stat = "identity")+ 
#   # scale_fill_manual(breaks = c("haematopoietic_neoplasm", "lymphoid_neoplasm"), 
#   #                   values = c("purple", "green")) +
#   labs(x= "Primary Site")+
#   scale_y_break(c(5000, 44900),
#                 ticklabels=c(45000, 46000, 47000),
#                 expand = FALSE,
#                 scales="fixed") +
#   ylim(0, 47000)+
#   coord_flip()
# dev.off()

# Calculate % CH in each cancers
a <- CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Sample.name, Primary.site, Primary.histology) %>%
  mutate(data = "CH variants")

b <- gnomad_variants %>% 
  filter(variant_in_cosmic == "Yes") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  select(Sample.name, Primary.site, Primary.histology) %>%
  mutate(data = "not CH variants") %>% 
  bind_rows(., a
  ) %>% 
  mutate(Primary.histology = case_when(
    Primary.site == "haematopoietic_and_lymphoid_tissue"    ~ Primary.histology,
    TRUE    ~ NA_character_
  ))

c <- b %>% 
  group_by(data, Primary.site, Primary.histology) %>%
    summarise(sum_individuals_by_primaryby_ch = n()) %>%
    ungroup() %>% 
  arrange(desc(Primary.site)) %>% 
  group_by(Primary.site) %>%
  mutate(sum_individuals_by_primary = sum(sum_individuals_by_primaryby_ch)) %>% 
  ungroup() %>% 
  mutate(perc = sum_individuals_by_primaryby_ch / sum_individuals_by_primary * 100)
write_rds(c, "New_Fig5B_percent_Prevalence_Cancers_inCosmic.rds")

pdf("Fig5B_Prevalence_Cancers_inCosmic.pdf", width = 7, height = 5
)
c %>% 
  filter(data == "CH variants") %>% 
  ggplot(aes(x= fct_reorder(Primary.site, perc),
             y=perc, fill= Primary.histology))+
  geom_bar(stat = "identity")+
  # scale_fill_manual(breaks = c("haematopoietic_neoplasm", "lymphoid_neoplasm"),
  #                   values = c("purple", "green")) +
  labs(x= "Primary Site")+
#   scale_y_break(c(5000, 44900),
#                 ticklabels=c(45000, 46000, 47000),
#                 expand = FALSE,
#                 scales="fixed") +
#   ylim(0, 47000)+
  coord_flip()
dev.off()




CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Primary.site) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Samples, careful if 1 sample shows multiple gene-mutations**") %>% 
  as_kable()

CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  select(Primary.histology) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics <br>for Primary.site == 'haematopoietic_and_lymphoid_tissue'**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Samples, careful if 1 sample shows multiple gene-mutations**") %>% 
  as_kable()

CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Mutation.Description,
         FATHMM.prediction, Mutation.somatic.status) %>% 
  mutate(Mutation.somatic.status = case_when(
    str_detect(Mutation.somatic.status, "somatic")        ~ "reported or confirmed somatic",
    TRUE                                                  ~ Mutation.somatic.status
  )) %>% 
  mutate(Mutation.Description = case_when(
    Mutation.Description == ""           ~ NA_character_,
    Mutation.Description == "Unknown"    ~ NA_character_,
    TRUE                                 ~ Mutation.Description
  )) %>% 
  mutate(FATHMM.prediction = case_when(
    FATHMM.prediction == ""              ~ NA_character_,
    FATHMM.prediction == "Unknown"       ~ NA_character_,
    TRUE                                 ~ FATHMM.prediction
  )) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Samples, careful if 1 sample shows multiple gene-mutations**") %>% 
  as_kable()


# Figure 5D Compare top 10-20 genes for CH vs CM----
# print("FIGURE 5 top 10-20 genes for CH vs CM")
# CH_variants %>% 
#   distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
#   select(SYMBOL) %>% 
#   mutate(data = "CH variants") %>% 
#   bind_rows(., clonal_mosaicism %>% 
#               filter(variant_in_cosmic == "Yes") %>% 
#               distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
#               select(SYMBOL) %>% 
#               mutate(data = "CM variants present in cosmic")) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency")) %>% 
#   bold_labels() %>% 
#   modify_header(
#     label = "**Compare top 10-20 genes for CH vs CM - Prevalence - individual in Cosmic**"
#   ) %>% 
#   as_kable()

print("AC for variants in bick per gene in CH")
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count_nfe, 
         # alt_allele_count_afr, alt_allele_count_amr, 
         alt_allele_count_eas, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_nfe)) %>% head(30)
print(a, n=30)
print("NC_AC for variants in bick per gene in CH")
write_csv(CH_variants %>% 
            filter(mutations_in_BickWHO == "Yes") %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(IDs, SYMBOL, nc_alt_allele_count_nfe, 
                   nc_alt_allele_count_afr, nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas, nc_alt_allele_count_sas,
                   nc_alt_allele_count,
                   alt_allele_count_nfe, 
                   alt_allele_count_afr, alt_allele_count_amr,
                   alt_allele_count_eas, alt_allele_count_sas,
                   alt_allele_count),
                   "NC_AC for variants in bick per gene in CH.csv")
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, nc_alt_allele_count_nfe, 
         # nc_alt_allele_count_afr, nc_alt_allele_count_amr, 
         nc_alt_allele_count_eas, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_nfe)) %>% head(30)
print(a, n=30)
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, nc_alt_allele_count_nfe, 
         # nc_alt_allele_count_afr, nc_alt_allele_count_amr, 
         nc_alt_allele_count_sas, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_sas)) %>% head(30)
print(a, n=30)
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, nc_alt_allele_count_nfe, 
         nc_alt_allele_count_afr, nc_alt_allele_count_amr) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_afr)) %>% head(30)
print(a, n=30)
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #nc_alt_allele_count_nfe, 
         # nc_alt_allele_count_afr, nc_alt_allele_count_amr,
         nc_alt_allele_count_eas, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_eas)) %>% head(30)
print(a, n=30)
a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #nc_alt_allele_count_nfe, 
         # nc_alt_allele_count_afr, 
         nc_alt_allele_count_amr,
         nc_alt_allele_count_afr, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_amr)) %>% head(30)
print(a, n=30)

print("AC for variants in bick per gene in CM")
write_csv(clonal_mosaicism %>% 
            filter(mutations_in_BickWHO == "Yes") %>% 
            distinct(IDs, .keep_all = TRUE) %>% 
            select(IDs, SYMBOL, nc_alt_allele_count_nfe, 
                   nc_alt_allele_count_afr, nc_alt_allele_count_amr,
                   nc_alt_allele_count_eas, nc_alt_allele_count_sas,
                   nc_alt_allele_count,
                   alt_allele_count_nfe, 
                   alt_allele_count_afr, alt_allele_count_amr,
                   alt_allele_count_eas, alt_allele_count_sas,
                   alt_allele_count),
          "NC_AC for variants in bick per gene in CM.csv")
print("NC_AC for variants in bick per gene in CM")
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, alt_allele_count_nfe, 
         # alt_allele_count_afr, alt_allele_count_amr, 
         alt_allele_count_eas, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_nfe)) %>% head(30)
print(a, n=30)
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #alt_allele_count_nfe, 
         # alt_allele_count_afr, alt_allele_count_amr, 
         alt_allele_count_sas, alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(alt_allele_count_sas)) %>% head(30)
print(a, n=30)
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #nc_alt_allele_count_nfe, 
         nc_alt_allele_count_afr, nc_alt_allele_count_amr) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_afr)) %>% head(30)
print(a, n=30)
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #nc_alt_allele_count_nfe, 
         # alt_allele_count_afr, alt_allele_count_amr,
         nc_alt_allele_count_eas, nc_alt_allele_count) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_eas)) %>% head(30)
print(a, n=30)
a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(SYMBOL, #nc_alt_allele_count_nfe, 
         nc_alt_allele_count_amr) %>% 
  group_by(SYMBOL) %>% 
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  arrange(desc(nc_alt_allele_count_amr)) %>% head(30)
print(a, n=30)

print("AC for genes in bick per gene in CH")
a <- CH_variants %>%
  group_by(SYMBOL) %>%
  fill(mutations_in_BickWHO, .direction = "updown") %>%
  ungroup() %>%
  filter(mutations_in_BickWHO == "Yes") %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(SYMBOL, alt_allele_count_nfe,
         # alt_allele_count_afr, alt_allele_count_amr,
         alt_allele_count_eas, alt_allele_count) %>%
  group_by(SYMBOL) %>%
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  arrange(desc(alt_allele_count_nfe)) %>% head(30)
print(a, n=30)

print("AC for genes in bick per gene in CM")
a <- clonal_mosaicism %>%
  group_by(SYMBOL) %>%
  fill(mutations_in_BickWHO, .direction = "updown") %>%
  ungroup() %>%
  filter(mutations_in_BickWHO == "Yes") %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(SYMBOL, alt_allele_count_nfe,
         # alt_allele_count_afr, alt_allele_count_amr,
         alt_allele_count_eas, alt_allele_count) %>%
  group_by(SYMBOL) %>%
  mutate(across(where(is.numeric), ~ sum(.))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  arrange(desc(alt_allele_count_nfe)) %>% head(30)
print(a, n=30)

# VennD genes----
print("data with Venndiagram for genes in pop")
a <- CH_variants %>% distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(SYMBOL))

print("NC_freq for nfe == 0")
a %>% filter(nc_freq_allele_count_nfe == 0 & 
               nc_freq_allele_count_afr > 0 & 
               nc_freq_allele_count_eas > 0 & 
               nc_freq_allele_count_amr > 0 & 
               nc_freq_allele_count_sas > 0) %>% 
  select(SYMBOL, nc_freq_allele_count_nfe,
         nc_freq_allele_count_afr, nc_freq_allele_count_eas,
         nc_freq_allele_count_amr, nc_freq_allele_count_sas) %>% 
  arrange(desc(nc_freq_allele_count_afr))

print("NC_freq for afr == 0")
a %>% filter(nc_freq_allele_count_nfe > 0 & 
               nc_freq_allele_count_afr == 0 & 
               nc_freq_allele_count_eas > 0 & 
               nc_freq_allele_count_amr > 0 & 
               nc_freq_allele_count_sas > 0) %>% 
  select(SYMBOL, nc_freq_allele_count_nfe,
         nc_freq_allele_count_afr, nc_freq_allele_count_eas,
         nc_freq_allele_count_amr, nc_freq_allele_count_sas) %>% 
  arrange(desc(nc_freq_allele_count_nfe))

print("NC_freq for eas == 0")
a %>% filter(nc_freq_allele_count_nfe > 0 & 
               nc_freq_allele_count_afr > 0 & 
               nc_freq_allele_count_eas == 0 & 
               nc_freq_allele_count_amr > 0 & 
               nc_freq_allele_count_sas > 0) %>% 
  select(SYMBOL, nc_freq_allele_count_nfe,
         nc_freq_allele_count_afr, nc_freq_allele_count_eas,
         nc_freq_allele_count_amr, nc_freq_allele_count_sas) %>% 
  arrange(desc(nc_freq_allele_count_afr))

print("NC_freq for amr == 0")
a %>% filter(nc_freq_allele_count_nfe > 0 & 
               nc_freq_allele_count_afr > 0 & 
               nc_freq_allele_count_eas > 0 & 
               nc_freq_allele_count_amr == 0 & 
               nc_freq_allele_count_sas > 0) %>% 
  select(SYMBOL, nc_freq_allele_count_nfe,
         nc_freq_allele_count_afr, nc_freq_allele_count_eas,
         nc_freq_allele_count_amr, nc_freq_allele_count_sas) %>% 
  arrange(desc(nc_freq_allele_count_afr))

print("NC_freq for sas == 0")
a %>% filter(nc_freq_allele_count_nfe > 0 & 
               nc_freq_allele_count_afr > 0 & 
               nc_freq_allele_count_eas > 0 & 
               nc_freq_allele_count_amr > 0 & 
               nc_freq_allele_count_sas == 0) %>% 
  select(SYMBOL, nc_freq_allele_count_nfe,
         nc_freq_allele_count_afr, nc_freq_allele_count_eas,
         nc_freq_allele_count_amr, nc_freq_allele_count_sas) %>% 
  arrange(desc(nc_freq_allele_count_afr))





# Fig5C COSMIC focus on 5 first cancer----
print("FIGURE focus on 5 first cancer")
CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site %in% 
           c("haematopoietic_and_lymphoid_tissue", 
             "large_intestine", "breast", "skin", "lung", "endometrium")) %>% 
  mutate(Primary.site = factor(
    Primary.site, levels = 
      c("haematopoietic_and_lymphoid_tissue", 
        "large_intestine", "breast", "skin", "lung", "endometrium"))) %>% 
  select(Primary.site, Age, Mutation.somatic.status,# SYMBOL
         ) %>% 
  mutate(Mutation.somatic.status = case_when(
    str_detect(Mutation.somatic.status, "somatic")        ~ "reported or confirmed somatic",
    TRUE                                                  ~ Mutation.somatic.status
  )) %>% 
  tbl_summary(by= Primary.site,
              sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% add_overall() %>% 
  add_p() %>% 
  as_kable()



a <- CH_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site %in% 
           c("haematopoietic_and_lymphoid_tissue", 
             "large_intestine", "endometrium", "stomach", "bone",
             "skin", "pancreas", "thyroid")) %>% 
  select(Primary.site, Age) %>% 
  mutate(data = "CH variants")

gnomad_variants %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>% 
  filter(Primary.site %in% 
           c("haematopoietic_and_lymphoid_tissue", 
             "large_intestine", "endometrium", "stomach", "bone",
             "skin", "pancreas", "thyroid")) %>% 
  select(Primary.site, Age) %>% 
  mutate(data = "Not CH variants") %>% 
  bind_rows(., a
  ) %>% 
  mutate(Primary.site = factor(
    Primary.site, levels = 
      c("haematopoietic_and_lymphoid_tissue", 
        "large_intestine", "endometrium", "stomach", "bone",
        "skin", "pancreas", "thyroid"))) %>% 
  tbl_strata(strata = Primary.site,
             .tbl_fun =
               ~ .x %>%
               tbl_summary(by = data) %>%
               add_n() %>% 
               add_p() %>% bold_p(t=.05)
  ) %>% 
  bold_labels() %>% #add_overall() %>% 
  # add_p() %>% bold_p(t=.05) %>% 
  modify_header(
    label = '**Cosmic Patients Characteristics** "haematopoietic_and_lymphoid_tissue", 
             "large_intestine", "endometrium", "stomach", "bone",
             "skin", "pancreas", "thyroid"'
  ) %>% 
  as_kable()

# Focus on HEM
a <- CH_variants %>%
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>%
  select(Age, Sample.name) %>%
  mutate(data = "CH variants")

# pdf("Fig5A_Age_in_HEM_COSMIC.pdf")
# gnomad_variants %>%
#   distinct(Sample.name, .keep_all = TRUE) %>%
#   filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>%
#   filter(!Sample.name %in% c(a$Sample.name)) %>%
#   select(Age) %>%
#   mutate(data = "not CH selected") %>%
#   bind_rows(., a
#   ) %>%
#   group_by(data) %>%
#   mutate(median_ = median(Age, na.rm = TRUE)) %>%
#   ggplot(aes(x=Age, fill= data))+
#   geom_histogram(binwidth = 1)+
#   scale_fill_manual(values = c("lightgrey", "purple"), name="")+
#   geom_vline(aes(xintercept= median_, color= data))+
#   scale_color_manual(values = c("#999999", "#990066"))+
#   guides(color = FALSE)
# dev.off()

gnomad_variants %>%
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>%
  filter(!Sample.name %in% c(a$Sample.name)) %>%
  select(Age) %>%
  mutate(data = "not selected") %>%
  bind_rows(., a
  ) %>%
  select(Age, data) %>%
  tbl_summary(by=data) %>%
  bold_labels() %>% add_stat_label() %>% add_overall() %>%
  add_p() %>% bold_p(t=.05) %>%
  modify_header(
    label = "**Cosmic HEM Patients Characteristics**"
  ) %>%
  as_kable()


# Fig6 ABick & WHO AB distribution----
print("Figure 6")
pdf("Fig6B_Average_AB_Density_Distribution_BickWho.pdf")
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(AB_distribution_mean, age_distribution_mean,
         mutations_in_BickWHO, variant_in_cosmic) %>% 
  bind_rows(., 
            CH_variants %>% 
              distinct(IDs, .keep_all = TRUE) %>% 
              select(AB_distribution_mean, age_distribution_mean,
                     variant_in_cosmic) %>% 
              mutate(mutations_in_BickWHO = "Overall")
  ) %>% 
  mutate(mutations_in_BickWHO = case_when(
    mutations_in_BickWHO == "Overall"   ~ "All CH",
    mutations_in_BickWHO == "Yes"       ~ "Curated M-CHIP list",
    mutations_in_BickWHO == "No"        ~ "Not Curated M-CHIP"
  )) %>% 
  mutate(mutations_in_BickWHO = factor(mutations_in_BickWHO, 
                                       levels= c("All CH", "Curated M-CHIP list", "Not Curated M-CHIP"))) %>% 
  ggplot(aes(x=AB_distribution_mean, y=mutations_in_BickWHO, fill=mutations_in_BickWHO)) +
  geom_density_ridges(alpha=0.5)+
  # ggtitle("Average Allele Balance Density Distribution")+
  labs(x= "Allele Balance Distribution",
       y= "Frequency Density", 
       fill= "Curated M-CH variant")+
  scale_fill_discrete(limits = c("Not Curated M-CHIP",
                                 "Curated M-CHIP list",
                                 "All CH"))+
  scale_x_continuous(limits=c(0, 0.5),
                     breaks=seq(0, 0.5, 0.1))+
  theme(axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()

pdf("Fig6D_Average_AB_Density_Distribution_L_CHIP.pdf")
CH_variants %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(AB_distribution_mean, age_distribution_mean,
         L_CHIP, variant_in_cosmic) %>% 
  bind_rows(., 
            CH_variants %>% 
              distinct(IDs, .keep_all = TRUE) %>% 
              select(AB_distribution_mean, age_distribution_mean,
                     variant_in_cosmic) %>% 
              mutate(L_CHIP = "Overall")
  ) %>% 
  mutate(L_CHIP = case_when(
    L_CHIP == "Overall"   ~ "All CH",
    L_CHIP == "Yes"       ~ "Curated L-CHIP list",
    L_CHIP == "No"        ~ "Not curated L-CHIP"
  )) %>% 
  mutate(L_CHIP = factor(L_CHIP, 
                         levels= c("All CH", "Curated L-CHIP list", "Not curated L-CHIP"))) %>% 
  ggplot(aes(x=AB_distribution_mean, y=L_CHIP, fill=L_CHIP)) +
  geom_density_ridges(alpha=0.5)+
  # ggtitle("Average Allele Balance Density Distribution")+
  labs(x= "Allele Balance Distribution",
       y= "Frequency Density", 
       fill= "Curated L-CH variant")+
  scale_fill_discrete(limits = c("Not curated L-CHIP",
                                 "Curated L-CHIP list",
                                 "All CH"))+
  scale_x_continuous(limits=c(0, 0.5),
                     breaks=seq(0, 0.5, 0.1),)+
  theme(axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()


pdf("FigS6B_SM_Average_AB_Density_Distribution_BickWho.pdf")
clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(AB_distribution_mean, age_distribution_mean,
         mutations_in_BickWHO, variant_in_cosmic) %>% 
  bind_rows(., 
            clonal_mosaicism %>% 
              distinct(IDs, .keep_all = TRUE) %>% 
              select(AB_distribution_mean, age_distribution_mean,
                     variant_in_cosmic) %>% 
              mutate(mutations_in_BickWHO = "Overall")
  ) %>% 
  mutate(mutations_in_BickWHO = case_when(
    mutations_in_BickWHO == "Overall"   ~ "All CH",
    mutations_in_BickWHO == "Yes"       ~ "Curated M-CHIP list",
    mutations_in_BickWHO == "No"        ~ "Not Curated M-CHIP"
  )) %>% 
  mutate(mutations_in_BickWHO = factor(mutations_in_BickWHO, 
                                       levels= c("All CH", "Curated M-CHIP list", "Not Curated M-CHIP"))) %>% 
  ggplot(aes(x=AB_distribution_mean, y=mutations_in_BickWHO, fill=mutations_in_BickWHO)) +
  geom_density_ridges(alpha=0.5)+
  # ggtitle("Average Allele Balance Density Distribution")+
  labs(x= "Allele Balance Distribution",
       y= "Frequency Density", 
       fill= "Curated M-CH variant")+
  scale_fill_discrete(limits = c("Not Curated M-CHIP",
                                 "Curated M-CHIP list",
                                 "All CH"))+
  scale_x_continuous(limits=c(0, 0.5),
                     breaks=seq(0, 0.5, 0.1))+
  theme(axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()

pdf("FigS6D_SM_Average_AB_Density_Distribution_L_CHIP.pdf")
clonal_mosaicism %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(AB_distribution_mean, age_distribution_mean,
         L_CHIP, variant_in_cosmic) %>% 
  bind_rows(., 
            clonal_mosaicism %>% 
              distinct(IDs, .keep_all = TRUE) %>% 
              select(AB_distribution_mean, age_distribution_mean,
                     variant_in_cosmic) %>% 
              mutate(L_CHIP = "Overall")
  ) %>% 
  mutate(L_CHIP = case_when(
    L_CHIP == "Overall"   ~ "All CH",
    L_CHIP == "Yes"       ~ "Curated L-CHIP list",
    L_CHIP == "No"        ~ "Not curated L-CHIP"
  )) %>% 
  mutate(L_CHIP = factor(L_CHIP, 
                         levels= c("All CH", "Curated L-CHIP list", "Not curated L-CHIP"))) %>% 
  ggplot(aes(x=AB_distribution_mean, y=L_CHIP, fill=L_CHIP)) +
  geom_density_ridges(alpha=0.5)+
  # ggtitle("Average Allele Balance Density Distribution")+
  labs(x= "Allele Balance Distribution",
       y= "Frequency Density", 
       fill= "Curated L-CH variant")+
  scale_fill_discrete(limits = c("Not curated L-CHIP",
                                 "Curated L-CHIP list",
                                 "All CH"))+
  scale_x_continuous(limits=c(0, 0.5),
                     breaks=seq(0, 0.5, 0.1),)+
  theme(axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()



CH_variants %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_variants %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_variants %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_variants %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()

CH_variants %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()

CH_variants %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with SM variants in Bick & WHO**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with SM variants in Bick & WHO**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with SM variants in Bick & WHO**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with SM variants in L_CHIP**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with SM variants in L_CHIP**"
  ) %>% 
  as_kable()

clonal_mosaicism %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(L_CHIP)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**SM In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with SM variants in L_CHIP**"
  ) %>% 
  as_kable()

# Fig 6C----
pdf("Fig6C_CH_Popultaion_AF_Density_Distribution_in_BickWHO.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>%
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color= Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in restricetd CH Variants list")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()

# Fig 6F----
pdf("Fig6F_CH_Popultaion_AF_Density_Distribution_in_L_CHIP.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(L_CHIP == "Yes") %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>%
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color= Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in restricetd CH Variants list")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()

legend <- CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>%
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color= Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  ggtitle("Popultaion Allele Frequency Distribution in restricetd CH Variants list")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  theme(legend.title=element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")+
  xlim(-0.00001,0.00015)
legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
pdf("Figure4_legend.pdf")
ggpubr::as_ggplot(legend)
dev.off()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas, nc_freq_allele_count_sas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In Bick+WHO?**") %>% 
  modify_header(
    label = "**CH Popultaion Allele Frequency in Bick+WHO**"
  ) %>% 
  as_kable()

CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas, nc_freq_allele_count_sas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         L_CHIP) %>% 
  tbl_summary(by= L_CHIP,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In L_CHIP?**") %>% 
  modify_header(
    label = "**CH Popultaion Allele Frequency in L_CHIP**"
  ) %>% 
  as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas, nc_freq_allele_count_sas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In Bick+WHO?**") %>% 
  modify_header(
    label = "**SM Popultaion Allele Frequency in Bick+WHO**"
  ) %>% 
  as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas, nc_freq_allele_count_sas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         L_CHIP) %>% 
  tbl_summary(by= L_CHIP,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In L_CHIP?**") %>% 
  modify_header(
    label = "**SM Popultaion Allele Frequency in L_CHIP**"
  ) %>% 
  as_kable()

# Fig 6B----
pdf("Fig6B_CM_Popultaion_AF_Density_Distribution_in_Bick_WHO.pdf")
clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>%
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color= Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in restricetd SM Variants list")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()

# Fig 6E----
pdf("Fig6E_CM_Popultaion_AF_Density_Distribution_in_L_CHIP.pdf")
clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(L_CHIP == "Yes") %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs,
         nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe
  ) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>%
  ggplot(aes(x=value, y=Ethnicity, fill=Ethnicity, color= Ethnicity)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in restricetd SM Variants list")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("Black",
                                 "East Asian",
                                 "Hispanic",
                                 "South Asian",
                                 "White"))+
  scale_color_discrete(limits = c("Black",
                                  "East Asian",
                                  "Hispanic",
                                  "South Asian",
                                  "White"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 45, 
        #                            vjust = 1,
        #                            hjust=1),
        axis.ticks.y = element_blank())
dev.off()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In Bick+WHO?**") %>% 
  modify_header(
    label = "**CM Popultaion Allele Frequency in Bick+WHO**"
  ) %>% 
  as_kable()

clonal_mosaicism %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(nc_freq_allele_count_afr,
         nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_female, nc_freq_allele_count_male,
         L_CHIP) %>% 
  tbl_summary(by= L_CHIP,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(0.05) %>% add_overall() %>% 
  modify_spanning_header(c(stat_1, stat_2) ~ "**In L_CHIP?**") %>% 
  modify_header(
    label = "**CM Popultaion Allele Frequency in L_CHIP**"
  ) %>% 
  as_kable()

# pdf("Fig6C_Popultaion_AF_Density_Distribution_by_gene.pdf")
# CH_variants %>%
#   distinct(IDs, .keep_all = TRUE) %>%
#   filter(mutations_in_BickWHO == "Yes") %>%
#   # filter(SYMBOL %in% c("DNMT3A", "SF3B1")) %>%
#   select(IDs, SYMBOL,
#          # nc_freq_allele_count, nc_freq_allele_count_female, nc_freq_allele_count_male,
#          nc_freq_allele_count_afr,# nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
#          nc_freq_allele_count_amr,# nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
#          nc_freq_allele_count_eas,# nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
#          nc_freq_allele_count_sas,# nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
#          nc_freq_allele_count_nfe#, nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
#   ) %>%
#   pivot_longer(cols = -c(IDs, SYMBOL)) %>%
#   mutate(value = as.numeric(value)) %>%
#   mutate(name = str_remove(name, "nc_freq_allele_count_"),
#          # name = str_replace(name, "nc_freq_allele_count", "overall"),
#          name = factor(name, levels = c(
#            # "male", "female", "overall",
#            "nfe_male", "nfe_female", "nfe",
#            "eas_male", "eas_female", "eas", 
#            "sas_male", "sas_female", "sas", 
#            "amr_male", "amr_female", "amr", 
#            "afr_male", "afr_female", "afr")
#          )) %>% 
#   mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
#     is.na(Sex)             ~ "overall",
#     TRUE                   ~ Sex
#   )) %>%
#   mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
#     is.na(Ethnicity)       ~ "overall",
#     TRUE                   ~ Ethnicity
#   )) %>%
#   filter(!is.na(value)) %>% 
#   # ggplot(aes(x=value, y=name, fill=Ethnicity, linetype= Sex)) +
#   # geom_density_ridges(alpha=0.5)+
#   ggplot(aes(x=value, y=name, fill=SYMBOL, color= SYMBOL)) +
#   geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
#   ggtitle("Popultaion Allele Frequency Distribution - xlim(-0.00001,0.00015)")+
#   labs(x= "Popultaion Allele Frequency", y= NULL)+
#   xlim(-0.00001,0.00015)
#   # facet_wrap(. ~ SYMBOL, scales = "free_x")
# dev.off()

# CH_variants %>%
#   distinct(IDs, .keep_all = TRUE) %>%
#   filter(mutations_in_BickWHO == "Yes") %>% 
#   # filter(SYMBOL %in% c("DNMT3A", "SF3B1")) %>% 
#   select(nc_freq_allele_count_afr,
#          nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
#          nc_freq_allele_count_eas, nc_freq_allele_count_sas,
#          nc_freq_allele_count_female, nc_freq_allele_count_male,
#          SYMBOL) %>% 
#   tbl_summary(by = SYMBOL,
#               type = list(everything() ~ "continuous"),
#               digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
#   ) %>% 
#   bold_labels() %>% 
#   add_p() %>% bold_p(0.05) %>% add_overall() %>% 
#   modify_spanning_header(c(stat_1, stat_2) ~ "**Gene in Bick+WHO**") %>% 
#   modify_header(
#     label = "**Popultaion Allele Frequency by gene in Bick+WHO**"
#   ) %>% 
#   as_kable()

# pdf("Fig6C_CM_Popultaion_AF_Density_Distribution_by_gene.pdf")
# clonal_mosaicism %>%
#   distinct(IDs, .keep_all = TRUE) %>%
#   filter(mutations_in_BickWHO == "Yes") %>%
#   # filter(SYMBOL %in% c("DNMT3A", "SF3B1")) %>%
#   select(IDs, SYMBOL,
#          # nc_freq_allele_count, nc_freq_allele_count_female, nc_freq_allele_count_male,
#          nc_freq_allele_count_afr,# nc_freq_allele_count_afr_female, nc_freq_allele_count_afr_male,
#          nc_freq_allele_count_amr,# nc_freq_allele_count_amr_female, nc_freq_allele_count_amr_male,
#          nc_freq_allele_count_eas,# nc_freq_allele_count_eas_female, nc_freq_allele_count_eas_male,
#          nc_freq_allele_count_sas,# nc_freq_allele_count_sas_female, nc_freq_allele_count_sas_male,
#          nc_freq_allele_count_nfe,# nc_freq_allele_count_nfe_female, nc_freq_allele_count_nfe_male 
#   ) %>%
#   pivot_longer(cols = -c(IDs, SYMBOL)) %>%
#   mutate(value = as.numeric(value)) %>%
#   mutate(name = str_remove(name, "nc_freq_allele_count_"),
#          # name = str_replace(name, "nc_freq_allele_count", "overall"),
#          name = factor(name, levels = c(
#            # "male", "female", "overall",
#            "nfe_male", "nfe_female", "nfe",
#            "eas_male", "eas_female", "eas", 
#            "sas_male", "sas_female", "sas", 
#            "amr_male", "amr_female", "amr", 
#            "afr_male", "afr_female", "afr")
#          )) %>% 
#   mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
#     is.na(Sex)             ~ "overall",
#     TRUE                   ~ Sex
#   )) %>%
#   mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
#     is.na(Ethnicity)       ~ "overall",
#     TRUE                   ~ Ethnicity
#   )) %>%
#   # ggplot(aes(x=value, y=name, fill=Ethnicity, linetype= Sex)) +
#   # geom_density_ridges(alpha=0.5)+
#   ggplot(aes(x=value, y=name, fill=SYMBOL, color= SYMBOL)) +
#   geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
#   ggtitle("Popultaion Allele Frequency Distribution - xlim(-0.00001,0.00015)")+
#   labs(x= "Popultaion Allele Frequency", y= NULL)+
#   xlim(-0.00001,0.00015)
#   # facet_wrap(. ~ SYMBOL, scales = "free_x")
# dev.off()

# clonal_mosaicism %>%
#   distinct(IDs, .keep_all = TRUE) %>%
#   filter(mutations_in_BickWHO == "Yes") %>% 
#   # filter(SYMBOL %in% c("DNMT3A", "SF3B1")) %>% 
#   select(nc_freq_allele_count_afr,
#          nc_freq_allele_count_amr, nc_freq_allele_count_nfe,
#          nc_freq_allele_count_eas, nc_freq_allele_count_sas,
#          nc_freq_allele_count_female, nc_freq_allele_count_male,
#          SYMBOL) %>% 
#   tbl_summary(by = SYMBOL,
#               type = list(everything() ~ "continuous"),
#               digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
#   ) %>% 
#   bold_labels() %>% 
#   add_p() %>% bold_p(0.05) %>% add_overall() %>% 
#   modify_spanning_header(c(stat_1, stat_2) ~ "**Gene in Bick+WHO**") %>% 
#   modify_header(
#     label = "**CM Popultaion Allele Frequency by gene in Bick+WHO**"
#   ) %>% 
#   as_kable()

# Fig 7----
# pdf("Fig7A_Number_of_overall_COSMIC_samples_per_gene.pdf")
# CH_variants %>% 
#   distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
#   ggplot(aes(x= fct_infreq(SYMBOL), fill= SYMBOL))+
#   geom_bar()+
#   ggtitle("Number of overall COSMIC samples per gene In Bick+WHO",
#           subtitle = "Warning patient counted multiple time if had both gene")+
#   scale_fill_discrete(name = NULL)+
#   labs(x= NULL, y="Number of samples")+
#   coord_flip()
# dev.off()
# 
# pdf("Fig7B_Number_of_overall_COSMIC_samples_per_gene_In_BickWHO.pdf")
# CH_variants %>% 
#   filter(mutations_in_BickWHO == "Yes") %>% 
#   distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
#   select(Sample.name, SYMBOL) %>% 
#   group_by(Sample.name) %>% 
#   mutate(n = n()) %>% 
#   ggplot(aes(x= fct_infreq(SYMBOL), fill= SYMBOL))+
#   geom_bar()+
#   ggtitle("Number of overall COSMIC samples per gene In Bick+WHO",
#           subtitle = "Warning patient counted multiple time if had both gene")+
#   scale_fill_discrete(name = NULL)+
#   labs(x= NULL, y="Number of samples")+
#   coord_flip()
# dev.off()

# pdf("Fig7C_CM_Number_of_overall_COSMIC_samples_per_gene.pdf")
# clonal_mosaicism %>% 
#   distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
#   ggplot(aes(x= fct_infreq(SYMBOL), fill= SYMBOL))+
#   geom_bar()+
#   ggtitle("CM Number of overall COSMIC samples per gene In Bick+WHO",
#           subtitle = "Warning patient counted multiple time if had both gene")+
#   scale_fill_discrete(name = NULL)+
#   labs(x= NULL, y="Number of samples")+
#   coord_flip()
# dev.off()

a <- CH_variants %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
  select(Sample.name, SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  mutate(n = n()) %>% 
  arrange(desc(n), Sample.name)

write_delim(a, "bickwho_ch_variants_cosmic_samples_gene_list.txt")

write_delim(a %>% distinct(SYMBOL), "bickwho_ch_variants_gene_list.txt")

a <- clonal_mosaicism %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  distinct(Sample.name, SYMBOL, .keep_all = TRUE) %>% 
  select(Sample.name, SYMBOL) %>% 
  group_by(SYMBOL) %>% 
  mutate(n = n()) %>% 
  arrange(desc(n), Sample.name)

write_delim(a, "bickwho_SM_variants_cosmic_samples_gene_list.txt")

write_delim(a %>% distinct(SYMBOL), "bickwho_SM_variants_gene_list.txt")

# VAF cutoff

pdf("Popultaion Allele Frequency density gnomAD.pdf")
gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs, nc_freq_allele_count
  ) %>%
  ggplot(aes(x=nc_freq_allele_count)) +
  geom_density(color = "#F8766D")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+

  scale_x_continuous(limits = c(-0.00001,0.0005), 
                     labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf("Popultaion Allele Frequency density CH.pdf")
CH_variants %>%
  distinct(IDs, .keep_all = TRUE) %>% 
  unite(IDs, c(IDs, SYMBOL)) %>% 
  select(IDs, nc_freq_allele_count
  ) %>%
  ggplot(aes(x=nc_freq_allele_count)) +
  geom_density(color = "#00BFC4")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
  
  scale_x_continuous(limits = c(-0.00001,0.0005), 
                     labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf("Popultaion Allele Frequency count ylim10000.pdf")
gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, nc_freq_allele_count) %>%
  mutate(data = "Variants considered") %>%
  
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, nc_freq_allele_count) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  mutate(data = factor(data, levels = c(
    "Variants considered",
    "Clonal hematopoiesis (CH) variants"
    ))) %>% 
  
  mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
  filter(!is.na(nc_freq_allele_count)) %>%
  group_by(nc_freq_allele_count, data) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  
  ggplot(aes(x= nc_freq_allele_count, y=count, fill= data
  ))+
  geom_area(position = "identity") +
  scale_x_continuous(limits = c(0,0.02), labels = function(x) format(x, scientific = TRUE)
                     )+
  coord_cartesian(ylim = c(0, 10000))+
  theme(legend.position = "bottom")
dev.off()

pdf("Popultaion Allele Frequency count zoomed.pdf")
gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, nc_freq_allele_count) %>%
  mutate(data = "Variants considered") %>%
  
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, nc_freq_allele_count) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  mutate(data = factor(data, levels = c(
    "Variants considered",
    "Clonal hematopoiesis (CH) variants"
  ))) %>% 
  
  mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
  filter(!is.na(nc_freq_allele_count)) %>%
  group_by(nc_freq_allele_count, data) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  
  ggplot(aes(x= nc_freq_allele_count, y=count, fill= data
  ))+
  geom_area(position = "identity") +
  scale_x_continuous(limits = c(0,0.02), labels = function(x) format(x, scientific = TRUE)
  )+
  theme(legend.position = "bottom")+
  facet_zoom(ylim = c(0, 10000))
dev.off()  

pdf("Popultaion Allele Frequency density zoomed.pdf")
gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, nc_freq_allele_count) %>%
  mutate(data = "Variants considered") %>%
  
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, nc_freq_allele_count) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  mutate(data = factor(data, levels = c(
    "Variants considered",
    "Clonal hematopoiesis (CH) variants"
  ))) %>% 
  
  mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
  filter(!is.na(nc_freq_allele_count)) %>%
  group_by(nc_freq_allele_count, data) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  
  ggplot(aes(x= nc_freq_allele_count, fill= data
  ))+
  # ggplot(aes(x=nc_freq_allele_count)) +
  geom_density(alpha = 0.5)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
  
  scale_x_continuous(limits = c(-0.00001,1),# labels = function(x) format(x, scientific = TRUE)
                     breaks = c(seq(0, 1, by= 0.05), 0.01, 0.02)
                     )+
  # scale_y_discrete(expand = c(0, 0.2))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  facet_zoom(xlim = c(-0.00001, 0.2))
dev.off()

a <- gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, nc_freq_allele_count) %>%
  mutate(data = "Variants considered") %>%
  
  bind_rows(., CH_variants %>%
              distinct(IDs, .keep_all = TRUE) %>%
              select(IDs, nc_freq_allele_count) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  mutate(data = factor(data, levels = c(
    "Variants considered",
    "Clonal hematopoiesis (CH) variants"
  )))
write_rds(a, "population allele frequency data for VAF cut-point.rds")


# For age

a <- gnomad_variants %>%
  select(IDs, Sample.name, Age, Primary.site, mutations_in_BickWHO, L_CHIP, variant_in_cosmic) %>%
  mutate(data = "Variants considered") %>%
  
  bind_rows(., CH_variants %>%
              select(IDs, Sample.name, Age, Primary.site, mutations_in_BickWHO, L_CHIP, variant_in_cosmic) %>%
              mutate(data = "Clonal hematopoiesis (CH) variants")) %>%
  mutate(data = factor(data, levels = c(
    "Variants considered",
    "Clonal hematopoiesis (CH) variants"
  )))
write_rds(a, "age.rds")

