# Import library
library(tidyverse)
library(GenomicRanges)

# Load data
# consequence_ranking <- 
#   read_csv(
#     "/share/lab_gillis/Christelle/gnomAD_skweness/other_data/gnomAD_mut_consequence_ranking.csv")
gencode <- 
  read.delim(
    "/share/lab_gillis/Christelle/gnomAD_skweness/other_data/gencode_filtered.vcf.gz",
    header = TRUE, sep = " ")

black_list <- 
  read.delim(paste0(here::here(), "/Christelle/gnomAD_skweness/analysis/data/blacklist_hg38_ENCFF356LFX.bed.gz"), header=FALSE) %>% 
  `colnames<-`(c("seqnames", "chromStart", "chromEnd"))
centromeres <- 
  read.delim(paste0(here::here(), "/Christelle/gnomAD_skweness/analysis/data/Centromeres_hg38.gz"))
segmental_dup <- 
  read.delim(paste0(here::here(), "/Christelle/gnomAD_skweness/analysis/data/Segmental_Dups_hg38.gz"))
wm_sdust <- 
  read.delim(paste0(here::here(), "/Christelle/gnomAD_skweness/analysis/data/WM_SDust_hg38.gz"))


# CREATE FUNCTION
noise_filter <- function(gnomad){
  
  print(unique(gnomad$X.CHROM))
  
  gnomad_vep <- gnomad %>%
  # Extract SYMBOL from gnomAD VEP
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2]) %>%
  # Evaluate the number of vep in each vep string
  mutate(number_of_vep = sapply(strsplit(ens_vep, ","), length)) %>%
  # Separate each vep
  separate(col = ens_vep, paste("ens_vep", 1:max(.$number_of_vep), sep="_"),
           sep = ",", remove = T, extra = "warn", fill = "right") %>%
  # pivot longer VEPs
  pivot_longer(cols = starts_with("ens_vep_"), names_to = "ENSVEP", values_to = "ens_vep") %>%
  drop_na(ens_vep) %>% select(-ENSVEP)
  
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

  gnomad_vep <- gnomad_vep %>%
    # Extract Consequence and SYMBOL
    separate(col = ens_vep, into = ens_vep_var_names,
             sep = "\\|", remove = F, extra = "warn", fill = "right") %>%
    select(IDs, SYMBOL) %>% 
    filter(!is.na(SYMBOL)) %>% 
    distinct(IDs, .keep_all = TRUE)

  gnomad_vep <- gnomad_vep %>%
    # Add rank to variants Ids
    # left_join(., consequence_ranking,
    #           by = c("Consequence" = "mut_consequence")) %>%
    # arrange(IDs, Rank) %>%
    # distinct(IDs, .keep_all = TRUE) %>%
    # select(IDs, SYMBOL, Consequence, consequence_rank = Rank, consequence_group) %>%
    # Add protein coding gene var
    left_join(., gencode %>%
                dplyr::rename(gencode_gene_type = gene_type),
              by = c("SYMBOL" = "gene_name")) %>%
    mutate(is_gencode_protein_coding_gene = case_when(
      is.na(gencode_gene_type)     ~ "No",
      !is.na(gencode_gene_type)    ~ "Yes"
    ))
  
  # Merge with initial gnomAD
  gnomad <- left_join(gnomad,
                      gnomad_vep,
                      by = "IDs")
  # Cleaning
  rm(gnomad_vep)
  
  # Add cleaning filter
  dat <- gnomad %>% select(X.CHROM, chromStart = POS) %>% 
    mutate(chromEnd = chromStart + 0) %>% 
    mutate(seqnames = paste0("chr", X.CHROM))
  dat <- makeGRangesFromDataFrame(dat)
  df <- makeGRangesFromDataFrame(black_list,
                                 start.field="chromStart",
                                 end.field="chromEnd")
  black_list <- countOverlaps(dat, df) %>% as_tibble()

  df <- makeGRangesFromDataFrame(centromeres, 
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  centromeres <- countOverlaps(dat, df) %>% as_tibble()
  
  df <- makeGRangesFromDataFrame(segmental_dup,
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  segmental_dup <- countOverlaps(dat, df) %>% as_tibble()
  
  df <- makeGRangesFromDataFrame(wm_sdust,
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  wm_sdust <- countOverlaps(dat, df) %>% as_tibble()
  
  gnomad <- gnomad %>% 
    bind_cols(., black_list %>% dplyr::rename(black_list = value),
              centromeres %>% dplyr::rename(centromeres = value),
              segmental_dup %>% dplyr::rename(segmental_duplication = value),
              wm_sdust %>% dplyr::rename(wm_sdust = value))
  
  # # Join gnomad and cosmic for each chromosome file
  # selected_variants <- left_join(gnomad, 
  #                               cosmic, 
  #                               by = c("IDs", "X.CHROM", "POS", "REF", "ALT"
  #                               )) %>% 
  #   mutate(variant_in_cosmic = case_when(
  #     is.na(variant_in_cosmic)     ~ "No",
  #     variant_in_cosmic == "Yes"   ~ "Yes"
  #   ))
  # 
  # # Extract allele info
  # CH_variants <- selected_variants %>%
  #   # Generate allele count variables from the INFO var
  #   mutate(alt_allele_count = str_match(INFO, "AC=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_female = str_match(INFO, "AC_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_male = str_match(INFO, "AC_male=(.*?);")[,2]) %>%
  # 
  #   mutate(alt_allele_count_afr = str_match(INFO, "AC_afr=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_afr_female = str_match(INFO, "AC_afr_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_afr_male = str_match(INFO, "AC_afr_male=(.*?);")[,2]) %>%
  # 
  #   mutate(alt_allele_count_amr = str_match(INFO, "AC_amr=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_amr_female = str_match(INFO, "AC_amr_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_amr_male = str_match(INFO, "AC_amr_male=(.*?);")[,2]) %>%
  #   
  #   mutate(alt_allele_count_eas = str_match(INFO, "AC_eas=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_eas_female = str_match(INFO, "AC_eas_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_eas_male = str_match(INFO, "AC_eas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(alt_allele_count_sas = str_match(INFO, "AC_sas=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_sas_female = str_match(INFO, "AC_sas_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_sas_male = str_match(INFO, "AC_sas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(alt_allele_count_nfe = str_match(INFO, "AC_nfe=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_nfe_female = str_match(INFO, "AC_nfe_female=(.*?);")[,2]) %>%
  #   mutate(alt_allele_count_nfe_male = str_match(INFO, "AC_nfe_male=(.*?);")[,2]) %>%
  #   
  #   mutate(alt_allele_count_nfe_nwe = str_match(INFO, "AC_nfe_nwe=(.*?);")[,2]) %>%
  # 
  #   # Generate allele freq variables from the INFO var
  #   mutate(freq_allele_count = str_match(INFO, "AF=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_female = str_match(INFO, "AF_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_male = str_match(INFO, "AF_male=(.*?);")[,2]) %>%
  # 
  #   mutate(freq_allele_count_afr = str_match(INFO, "AF_afr=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_afr_female = str_match(INFO, "AF_afr_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_afr_male = str_match(INFO, "AF_afr_male=(.*?);")[,2]) %>%
  # 
  #   mutate(freq_allele_count_amr = str_match(INFO, "AF_amr=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_amr_female = str_match(INFO, "AF_amr_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_amr_male = str_match(INFO, "AF_amr_male=(.*?);")[,2]) %>%
  #   
  #   mutate(freq_allele_count_eas = str_match(INFO, "AF_eas=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_eas_female = str_match(INFO, "AF_eas_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_eas_male = str_match(INFO, "AF_eas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(freq_allele_count_sas = str_match(INFO, "AF_sas=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_sas_female = str_match(INFO, "AF_sas_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_sas_male = str_match(INFO, "AF_sas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(freq_allele_count_nfe = str_match(INFO, "AF_nfe=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_nfe_female = str_match(INFO, "AF_nfe_female=(.*?);")[,2]) %>%
  #   mutate(freq_allele_count_nfe_male = str_match(INFO, "AF_nfe_male=(.*?);")[,2]) %>%
  #   
  #   mutate(freq_allele_count_nfe_nwe = str_match(INFO, "AF_nfe_nwe=(.*?);")[,2]) %>%
  # 
  # 
  # 
  #   ## NON CANCER
  #   # Generate allele count variables from the INFO var
  #   mutate(nc_alt_allele_count = str_match(INFO, "non_cancer_AC=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_female = str_match(INFO, "non_cancer_AC_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_male = str_match(INFO, "non_cancer_AC_male=(.*?);")[,2]) %>%
  # 
  #   mutate(nc_alt_allele_count_afr = str_match(INFO, "non_cancer_AC_afr=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_afr_female = str_match(INFO, "non_cancer_AC_afr_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_afr_male = str_match(INFO, "non_cancer_AC_afr_male=(.*?);")[,2]) %>%
  # 
  #   mutate(nc_alt_allele_count_amr = str_match(INFO, "non_cancer_AC_amr=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_amr_female = str_match(INFO, "non_cancer_AC_amr_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_amr_male = str_match(INFO, "non_cancer_AC_amr_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_alt_allele_count_eas = str_match(INFO, "non_cancer_AC_eas=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_eas_female = str_match(INFO, "non_cancer_AC_eas_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_eas_male = str_match(INFO, "non_cancer_AC_eas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_alt_allele_count_sas = str_match(INFO, "non_cancer_AC_sas=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_sas_female = str_match(INFO, "non_cancer_AC_sas_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_sas_male = str_match(INFO, "non_cancer_AC_sas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_alt_allele_count_nfe = str_match(INFO, "non_cancer_AC_nfe=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_nfe_female = str_match(INFO, "non_cancer_AC_nfe_female=(.*?);")[,2]) %>%
  #   mutate(nc_alt_allele_count_nfe_male = str_match(INFO, "non_cancer_AC_nfe_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_alt_allele_count_nfe_nwe = str_match(INFO, "non_cancer_AC_nfe_nwe=(.*?);")[,2]) %>%
  # 
  #   # Generate allele count variables from the INFO var
  #   mutate(nc_freq_allele_count = str_match(INFO, "non_cancer_AF=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_female = str_match(INFO, "non_cancer_AF_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_male = str_match(INFO, "non_cancer_AF_male=(.*?);")[,2]) %>%
  # 
  #   mutate(nc_freq_allele_count_afr = str_match(INFO, "non_cancer_AF_afr=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_afr_female = str_match(INFO, "non_cancer_AF_afr_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_afr_male = str_match(INFO, "non_cancer_AF_afr_male=(.*?);")[,2]) %>%
  # 
  #   mutate(nc_freq_allele_count_amr = str_match(INFO, "non_cancer_AF_amr=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_amr_female = str_match(INFO, "non_cancer_AF_amr_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_amr_male = str_match(INFO, "non_cancer_AF_amr_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_freq_allele_count_eas = str_match(INFO, "non_cancer_AF_eas=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_eas_female = str_match(INFO, "non_cancer_AF_eas_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_eas_male = str_match(INFO, "non_cancer_AF_eas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_freq_allele_count_sas = str_match(INFO, "non_cancer_AF_sas=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_sas_female = str_match(INFO, "non_cancer_AF_sas_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_sas_male = str_match(INFO, "non_cancer_AF_sas_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_freq_allele_count_nfe = str_match(INFO, "non_cancer_AF_nfe=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_nfe_female = str_match(INFO, "non_cancer_AF_nfe_female=(.*?);")[,2]) %>%
  #   mutate(nc_freq_allele_count_nfe_male = str_match(INFO, "non_cancer_AF_nfe_male=(.*?);")[,2]) %>%
  #   
  #   mutate(nc_freq_allele_count_nfe_nwe = str_match(INFO, "non_cancer_AF_nfe_nwe=(.*?);")[,2]) %>%
  # 
  #   # "Total number of alleles in samples in the non_cancer subset, before removing low-confidence genotypes">
  #   mutate(nc_freq_allele_count_raw = str_match(INFO, "non_cancer_AF_raw=(.*?);")[,2]) %>%
  # 
  #   # Depth of informative coverage for each sample
  #   mutate(depth_allele = str_match(INFO, "DP=(.*?);")[,2]) %>%
  # 
  #   # allele type
  #   mutate(allele_type = str_match(INFO, "allele_type=(.*?);")[,2]) %>%
  #   # variant type
  #   mutate(variant_type = str_match(INFO, "variant_type=(.*?);")[,2]) %>%
  #   # Number=A,Type=String,Description="Population with maximum AF in the controls subset">
  #   mutate(controls_popmax = str_match(INFO, "controls_popmax=(.*?);")[,2]) %>% 
  #   mutate(popmax = str_match(INFO, ";popmax=(.*?);")[,2]) %>% 
    
  gnomad <- gnomad %>% 
    # Create filter for removing noise in AB
    mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
    separate(col = ab_hist_alt_bin_freq,
             into = c("0.00-0.05"),
             sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
    mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
    mutate(contaminated_AB = case_when(
      `0.00-0.05` == total_allele_balance               ~ "sequencing noise",
      TRUE                                              ~ "clean sequencing"
    )) %>% 
    select(-c("0.00-0.05"))

  write_delim(gnomad, paste0(i, "noise_filter.vcf.gz"))
  
}


###################################################################### II ### Add var necessary for filtering noise variants
chr_vec <- paste0("chr", c("X", "Y", 
                           seq(1, 22, 1)
                           ), "_")

for (i in chr_vec){
  
  ### I ### List and load files
  file_list_gnomad <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants/allele_balance/age_distribution/variant_depth",
                                 pattern = i,
                                 recursive=FALSE,
                                 full.names = TRUE)
  
  gnomad <- do.call("rbind",lapply(Sys.glob(file_list_gnomad), read.delim,
                                   header = TRUE, sep = " ",
                                   colClasses=c("X.CHROM"="character"))) %>% 
    mutate(X.CHROM1 = X.CHROM,
           X.CHROM = str_match(X.CHROM1, "^(\\d+|X|Y)")[,1]) 
  
  # file_list_cosmic <- list.files(path = "/share/lab_gillis/Christelle/cosmic_raw_data/splitted_data/filter_cosmic_withIDs",
  #                                pattern = i,
  #                                recursive=FALSE,
  #                                full.names = TRUE)
  # 
  # cosmic <- read.delim(file_list_cosmic,
  #                      header = TRUE, sep = " ") %>%
  #   mutate(X.CHROM = as.character(chr))

  ### II ### CALL FUNCTION
  noise_filter(gnomad = gnomad#, cosmic = cosmic
              )

  # n = n + 1

}


