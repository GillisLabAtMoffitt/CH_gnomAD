############## Create list genes that have at least 1 variant with gene_type | transcript_type == "protein_coding"

# Import library
library(tidyverse)

# Load data

path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")

gencode <-
  read.delim(#"~/Downloads/gencode.v42.annotation.gtf.gz", 
    paste0(path, "/gnomAD_skweness/other_data/gencode.v42.annotation.gtf.gz"),
    comment.char="#",
    # sep = " "
    col.names = c("chromosome", "annotation source", "feature type",
                  "genomic start location", "genomic end location", 
                  "score(not used)", "genomic strand", "genomic phase (for CDS features)",
                  "additional information as key-value pairs")
  ) %>% 
  janitor::clean_names()

gencode <- gencode %>% 
  mutate(chrom = str_match(chromosome, "chr(.*)")[,2]) %>% 
  mutate(gene_id = str_match(additional_information_as_key_value_pairs, 
                             "gene_id (.*?);")[,2]) %>% 
  mutate(transcript_id = str_match(additional_information_as_key_value_pairs, 
                                   "transcript_id (.*?);")[,2]) %>% 
  mutate(gene_type = str_match(additional_information_as_key_value_pairs, 
                               "gene_type (.*?);")[,2]) %>% 
  mutate(gene_name = str_match(additional_information_as_key_value_pairs, 
                               "gene_name (.*?);")[,2]) %>% 
  mutate(transcript_type = str_match(additional_information_as_key_value_pairs, 
                                     "transcript_type (.*?);")[,2])

sum(is.na(gencode$gene_name))
table(gencode$gene_type)
table(gencode$transcript_type)

gencode1 <- gencode %>% 
  filter(gene_type == "protein_coding") %>% 
  select(gene_name, gene_type) %>% 
  distinct()

write_delim(gencode1, paste0(path, "/gnomAD_skweness/other_data/gencode_filtered.vcf.gz"))

# Effective exon size
library(GenomicFeatures)

gtf <- paste0(path, "/gnomAD_skweness/other_data/gencode.v42.annotation.gtf.gz")

txdb <- makeTxDbFromGFF(gtf)
write_rds(txdb, "txdb.rds")
y= exonsBy(txdb, "gene")
Z = sapply(y, function(y) sum(width(reduce(y))))
z <- Z %>% as_tibble(rownames = "gene_id")

exon_effective_length <- gencode %>% 
  filter(gene_type == "protein_coding") %>% 
  full_join(., 
            z, 
            by = "gene_id") %>% 
  dplyr::select(gene_name, gene_type, exon_effective_length = value)

write_delim(exon_effective_length, paste0(path, "/gnomAD_skweness/other_data/gencode_filtered.vcf.gz"))




