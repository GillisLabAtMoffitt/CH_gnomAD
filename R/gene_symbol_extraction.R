######## Create function to extract gene symbol from best consequence rank in gnomAD ######## 

# Import library
library(tidyverse)

# Load data
path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
gnomad <-
  read.delim(paste0(path, "/gnomAD_raw_data/splitted_data/age_comparison_pval/padjust_age_selected_variants/chr2_padjusted.vcf.gz"),
             sep = " ")
gencode <- read.delim(paste0(path, "/gnomAD_skweness/other_data/gencode_filtered.vcf.gz"),
                      header = TRUE, sep = " ")
consequence_ranking <- readxl::read_xlsx(
  paste0(path, "/gnomAD_skweness/other_data/gnomAD_mut_consequence_grps_12.12.22_NG_MT.xlsx"))

# Clean ranking
consequence_ranking <- consequence_ranking %>% 
  select(mut_consequence = `Mutation Consequence in gnomAD`, ConsequenceGroup, consequence_group = ...4) %>% 
  mutate(consequence_group = coalesce(consequence_group, ConsequenceGroup)) %>% 
  select(-ConsequenceGroup) %>% 
  full_join(., consequence_ranking %>% 
              select(`Row Labels`, Rank),
            by = c("consequence_group" = "Row Labels"))

write_csv(consequence_ranking, 
          paste0(path, "/gnomAD_skweness/other_data/gnomAD_mut_consequence_ranking.csv"))


# Extract data
gnomad_vep <- gnomad %>% 
  # Extract SYMBOL from gnomAD VEP
  mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2]) %>% 
  # Evaluate the number of vep in each vep string
  mutate(number_of_vep = sapply(strsplit(ens_vep, ","), length)) %>% 
  # Separate each vep
  separate(col = ens_vep, paste("ens_vep", 1:max(.$number_of_vep), sep="_"),
           sep = ",", remove = T, extra = "warn", fill = "right") %>%
  # pivot longer VEPs
  pivot_longer(cols = starts_with("ens_vep_"), names_to = NULL, values_to = "ens_vep") %>%
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

gnomad_vep <- gnomad_vep %>%
  # Extract Consequence and SYMBOL
  separate(col = ens_vep, into = ens_vep_var_names,
           sep = "\\|", remove = F, extra = "warn", fill = "right") %>% 
  select(IDs, Consequence, SYMBOL)

gnomad_vep1 <- gnomad_vep %>% 
  # Add rank to variants Ids
  left_join(., consequence_ranking,
            by = c("Consequence" = "mut_consequence")) %>% 
  arrange(IDs, Rank) %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  select(IDs, SYMBOL, Consequence, consequence_rank = Rank) %>% 
  # Add protein coding gene var 
  left_join(., gencode %>% 
              rename(gencode_gene_type = gene_type), 
            by = c("SYMBOL" = "gene_name")) %>% 
  mutate(is_gencode_protein_coding_gene = case_when(
    is.na(gencode_gene_type)     ~ "No",
    !is.na(gencode_gene_type)    ~ "Yes"
  ))

# Merge with initial gnomAD
gnomad <- left_join(gnomad, 
                    gnomad_vep, 
                    by = "IDs")











