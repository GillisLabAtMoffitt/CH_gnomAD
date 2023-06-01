## Import library
library(tidyverse)
library(gtsummary)
library(GenomicRanges)
library(gtsummary)
library(ggconsort)
library(viridis)
library(ggridges)

theme_gtsummary_compact()
theme_set(theme_classic())


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants",
                        pattern = "*.vcf.gz",
                        recursive=FALSE,
                        full.names = TRUE)

black_list <- 
  read.delim(paste0(here::here(), "/blacklist_hg38_ENCFF356LFX.bed.gz"), header=FALSE) %>% 
  `colnames<-`(c("seqnames", "chromStart", "chromEnd"))

dat <- gnomad_variants %>% select(X.CHROM, chromStart = POS) %>% 
  mutate(chromEnd = chromStart + 0) %>% 
  mutate(seqnames = paste0("chr", X.CHROM))

dat <- makeGRangesFromDataFrame(dat)
df <- makeGRangesFromDataFrame(black_list,
                               start.field="chromStart",
                               end.field="chromEnd")
black_list <- countOverlaps(dat, df) %>% as_tibble()

centromeres <- 
  read.delim(paste0(here::here(), "/Centromeres_hg38.gz"))
df <- makeGRangesFromDataFrame(centromeres, 
                               seqnames.field = "chrom",
                               start.field="chromStart",
                               end.field="chromEnd")
centromeres <- countOverlaps(dat, df) %>% as_tibble()
segmental_dup <- 
  read.delim(paste0(here::here(), "/Segmental_Dups_hg38.gz"))
df <- makeGRangesFromDataFrame(segmental_dup,
                               seqnames.field = "chrom",
                               start.field="chromStart",
                               end.field="chromEnd")
segmental_dup <- countOverlaps(dat, df) %>% as_tibble()

wm_sdust <- 
  read.delim(paste0(here::here(), "/WM_SDust_hg38.gz"))
df <- makeGRangesFromDataFrame(wm_sdust,
                               seqnames.field = "chrom",
                               start.field="chromStart",
                               end.field="chromEnd")
wm_sdust <- countOverlaps(dat, df) %>% as_tibble()

gnomad_variants <- gnomad_variants %>% 
  bind_cols(., black_list %>% rename("value" = "black_list"),
            centromeres %>% rename("value" = "centromeres"),
            segmental_dup %>% rename("value" = "segmental_duplication"),
            wm_sdust %>% rename("value" = "wm_sdust"))

bick_genes <-
  readxl::read_xlsx(paste0(here::here(), "/Chr2 w Bick & WHO CH filter.xlsx"))
bick_genes <- bick_genes %>% 
  filter(`In Bick & WHO` == "TRUE") %>% 
  select(IDs) %>% 
  mutate(CH_bick_genes = "IDs in Bick+WHO")
known_CH_gene_list <- c(unique(bick_genes$IDs))
vep_info <-
  readxl::read_xlsx(paste0(here::here(), "/Chr2 w Bick & WHO CH filter w annotation.xlsx"), 
                    na = ".")

considered_CH_variants <- 
  readxl::read_xlsx(paste0(here::here(), "/gene_ids.xlsx")) %>% 
  select(ids, Gene, Variant) %>% 
  tidyr::unite(col = "considered_variants", c(Gene : Variant), sep = " ")

gnomad_variants <- gnomad_variants %>% 
  left_join(., considered_CH_variants,
            by = c("IDs" = "ids"))

vep_info <- vep_info %>% 
  distinct() %>% 
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
  unite(new_HGVSc, c(new_HGVSc_2, new_HGVSc_4), sep = ">") %>% 
  select(-c(SYMBOL, `In Bick & WHO`))

gnomad_variants <- gnomad_variants %>% 
  mutate(mutations_in_BickWHO = case_when(
    IDs %in% known_CH_gene_list        ~ "Yes",
    TRUE                               ~ "No"
  ), mutations_in_BickWHO = factor(mutations_in_BickWHO, levels = c("No", "Yes"))) %>% 
  mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
  separate(col = ab_hist_alt_bin_freq,
           into = c("0.00-0.05"),
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
  mutate(contaminated_AB = case_when(
    `0.00-0.05` == total_allele_balance               ~ "sequencing noise",
    TRUE                                              ~ "clean sequencing"
  )) %>% 
  select(-c("0.00-0.05")) %>% 
  left_join(., vep_info,
            by= "IDs")

# cleaning
rm(black_list, centromeres, segmental_dup, wm_sdust, 
   dat, df, considered_CH_variants,
   known_CH_gene_list, bick_genes, vep_info)

###################################################################### II ### Apply filter----
clean_AB_variants <- gnomad_variants %>% 
  filter(contaminated_AB == "clean sequencing") %>% 
  filter(black_list == 0) %>% 
  filter(centromeres == 0) %>% 
  filter(segmental_duplication == 0) %>% 
  filter(wm_sdust == 0)
allele_depth <- 10
depth_variants <- clean_AB_variants %>% 
  filter(dp_het_median >= allele_depth)
pval_selection <- 0.05
AB_distribution_variants <- depth_variants %>% 
  filter(adjusted_AB_pval < pval_selection)
AB_hedges_selection <- "large"
hedges_variants <- AB_distribution_variants %>%
  filter(AB_effect_size_cat == AB_hedges_selection)
prot_coding <- hedges_variants %>%
  filter(is_gencode_protein_coding_gene == "Yes")
cosmic_var <- prot_coding %>% 
  filter(variant_in_cosmic == "Yes")
age_distribution_variants <- cosmic_var %>% 
  filter(adjusted_age_pval < pval_selection)
age_hedges_selection <- "negligible|small"
CH_variants <- age_distribution_variants %>%
  filter(!str_detect(age_effect_size_cat, age_hedges_selection))

print(paste("Final number of variants in all chromosome are", 
            nrow(CH_variants %>% distinct(IDs))))
write_csv(CH_variants %>% distinct(IDs, .keep_all = TRUE) %>% select(IDs, X.CHROM),
      "gnomad_variants_list_after_filtering.csv")


###################################################################### II ### Analysis----
# Fig1_A CONSORT----
study_cohorts <-
  gnomad_variants %>%
  distinct(IDs, .keep_all = TRUE) %>%
  cohort_start("**Variants considered**") %>%
  # Define cohorts using named expressions
  # Notice that you can use previously defined cohorts in subsequent steps
  cohort_define(

    clean_seq = .full %>%
      filter(contaminated_AB == "clean sequencing"),
    black_clean = clean_seq %>%
      filter(black_list == 0),
    centro_clean = black_clean %>%
      filter(centromeres == 0),
    dup_clean = centro_clean %>%
      filter(segmental_duplication == 0),
    wm_sdust_clean = dup_clean %>%
      filter(wm_sdust == 0),
    clean_AB_variants = wm_sdust_clean %>%
      filter(dp_het_median >= allele_depth),

    AB_BH = clean_AB_variants %>%
      filter(adjusted_AB_pval <= pval_selection),
    AB_variants = AB_BH %>%
      filter(AB_effect_size_cat == AB_hedges_selection),

    protcoding_variants = AB_variants  %>%
      filter(is_gencode_protein_coding_gene == "Yes"),

    cosmic_variants = protcoding_variants  %>%
      filter(variant_in_cosmic == "Yes"),

    age_BH = cosmic_variants  %>%
      filter(adjusted_age_pval <= pval_selection),
    CH_variants = age_BH %>%
      filter(!str_detect(age_effect_size_cat, age_hedges_selection)),

    CH_variants_end = CH_variants,

    # anti_join is useful for counting exclusions
    excluded2 = anti_join(.full, clean_AB_variants, by = "IDs"),
    excluded2_1 = anti_join(.full, clean_seq, by = "IDs"),
    excluded2_2 = anti_join(clean_seq, black_clean, by = "IDs"),
    excluded2_3 = anti_join(black_clean, centro_clean, by = "IDs"),
    excluded2_4 = anti_join(centro_clean, dup_clean, by = "IDs"),
    excluded2_5 = anti_join(dup_clean, wm_sdust_clean, by = "IDs"),
    excluded2_6 = anti_join(wm_sdust_clean, clean_AB_variants, by = "IDs"),

    excluded1 = anti_join(clean_AB_variants, AB_variants, by = "IDs"),
    excluded1_1 = anti_join(clean_AB_variants, AB_BH, by = "IDs"),
    excluded1_2 = anti_join(AB_BH, AB_variants, by = "IDs"),

    excluded3 = anti_join(AB_variants, protcoding_variants, by = "IDs"),
    excluded5 = anti_join(protcoding_variants, cosmic_variants, by = "IDs"),

    excluded4 = anti_join(cosmic_variants, CH_variants, by = "IDs"),
    excluded4_1 = anti_join(cosmic_variants, age_BH, by = "IDs"),
    excluded4_2 = anti_join(age_BH, CH_variants, by = "IDs")
  ) %>%
  # Provide text labels for cohorts
  cohort_label(
    excluded2 = "<span style='color:darkblue'>Exclude variants with sequencing artifacts</span>",
    excluded2_1 = "Exclude variants with sequencing artifacts",
    excluded2_2 = "Exclude 'curated “blacklist” developed by ENCODE (73), and excluded mutations in the blacklisted regions'",
    excluded2_3 = "Exclude 'Centromeres (“Centromeres” in the “Mapping and Sequencing” group)'",
    excluded2_4 = "Exclude 'Segmental duplications (“segmental dups”)'",
    excluded2_5 = "Exclude 'Low complexity regions (“WM + SDust”) track'",
    excluded2_6 = "Exclude variants with depth inferieur to 10",
    excluded1 = "<span style='color:darkblue'>Exclude germline variants</span>",
    excluded1_1 = "BH pval < 0.05",
    excluded1_2 = "Hedges g == 'large'",
    excluded3 = "<span style='color:darkblue'>Select for protein coding mutations</span>",
    excluded5 = "<span style='color:darkblue'>Select variants present in Cosmic</span>",
    excluded4 = "<span style='color:darkblue'>Exclude non-age-skewed mutations</span>",
    excluded4_1 = "BH pval < 0.05",
    excluded4_2 = "Hedges g == 'large' or 'medium'",
    clean_AB_variants = "**Variants after quality control filtering**",
    AB_variants = "**Peripheral blood somatic mutations**",
    protcoding_variants = "**Protein-coding somatic mutations**",
    cosmic_variants = "**Candidate variants in COSMIC**",
    CH_variants_end = "**Age-dependent peripheral blood somatic mutations**"
  )
study_consort <- study_cohorts %>%
  consort_box_add(
    # top eligibility box at 40 height
    "full", 0.05, 50, cohort_count_adorn(study_cohorts, .full)
  ) %>%
    # first top excluded box at 0.2 (right from 0.05 eligibility box)
  consort_box_add(
    "exc_contamination", 0.2, 45, glue::glue(
      "{cohort_count_adorn(study_cohorts, excluded2)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_1)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_2)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_3)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_4)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_5)}<br>
      • {cohort_count_adorn(study_cohorts, excluded2_6)}
      ")
  ) %>%
  consort_box_add(
    "clean_AB_variants", 0.05, 40, cohort_count_adorn(study_cohorts, clean_AB_variants)
  ) %>%
    consort_box_add(
    "exc_germline", 0.2, 35, glue::glue(
      "{cohort_count_adorn(study_cohorts, excluded1)}<br>
      • {cohort_count_adorn(study_cohorts, excluded1_1)}<br>
      • {cohort_count_adorn(study_cohorts, excluded1_2)}
      ")
  ) %>%
  consort_box_add(
    "AB_variants", 0.05, 30, cohort_count_adorn(study_cohorts, AB_variants)
  ) %>%
  consort_box_add(
    "exc_non_protein_coding", 0.2, 25, glue::glue(
      "{cohort_count_adorn(study_cohorts, excluded3)}<br>
      • include variants present in a protein coding gene
      ")
  ) %>%
  consort_box_add(
    "protcoding_variants", 0.05, 20, cohort_count_adorn(study_cohorts, protcoding_variants)
  ) %>%
  consort_box_add(
    "exc_non_cosmic", 0.2, 15, glue::glue(
      "{cohort_count_adorn(study_cohorts, excluded5)}<br>
      • Select variants present in Cosmic
      ")
  ) %>%
  consort_box_add(
    "cosmic_variants", 0.05, 10, cohort_count_adorn(study_cohorts, cosmic_variants)
  ) %>%
  consort_box_add(
    "exc_non_age_skew", 0.2, 5, glue::glue(
      "{cohort_count_adorn(study_cohorts, excluded4)}<br>
      • {cohort_count_adorn(study_cohorts, excluded4_1)}<br>
      • {cohort_count_adorn(study_cohorts, excluded4_2)}
      ")
  ) %>%
  # Add bottom box
  consort_box_add(
    "CH_variants", 0.05, 0, cohort_count_adorn(study_cohorts, CH_variants_end)
  ) %>%
  # Add exclusion arrows
  consort_arrow_add(
    end = "exc_contamination", end_side = "left", start_x = 0.05, start_y = 45
  ) %>%
  consort_arrow_add(
    end = "exc_germline", end_side = "left", start_x = 0.05, start_y = 35
  ) %>%
  consort_arrow_add(
    end = "exc_non_protein_coding", end_side = "left", start_x = 0.05, start_y = 25
  ) %>%
  consort_arrow_add(
    end = "exc_non_cosmic", end_side = "left", start_x = 0.05, start_y = 15
  ) %>%
  consort_arrow_add(
    end = "exc_non_age_skew", end_side = "left", start_x = 0.05, start_y = 5
  ) %>%
  # Add top to bottom arrow
  consort_arrow_add(
    end = "CH_variants", end_side = "top", start_x = 0.05, start_y = 50
  )
jpeg("consort diagram.jpeg", width = 950, height = 650)
study_consort %>%
  ggplot() +
  geom_consort() + #theme_classic()
  theme_consort(margin_h = c(1,12), # bottom, left
                margin_v = c(1,40)) # top, right
dev.off()

pdf("consort diagram all chromosome.pdf")
study_consort %>%
  ggplot() +
  geom_consort() + #theme_classic()
  theme_consort(margin_h = c(1,12), # bottom, left
                margin_v = c(1,40)) # top, right
dev.off()



