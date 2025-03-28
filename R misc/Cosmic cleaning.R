# Import library
library(tidyverse)

# Load data
path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
# From COSMIC ANNOVAR
annovar_cosmic_data <-
  read.delim(paste0(path, "/cosmic_raw_data/hg38_cosmic70.txt")) %>% 
  `colnames<-`(c("#Chr","Start","End","Ref","Alt","INFO"))
# From COSMIC website
cosmic_data <-
  read.delim(paste0(path, "/cosmic_raw_data/CosmicMutantExport.tsv.gz")) #%>% 
  # `colnames<-`(c("#Chr","Start","End","Ref","Alt","INFO"))

head(cosmic_data, n = 30)

# cosmic_data <-
#   read.delim(paste0(path, "/cosmic_raw_data/CosmicMutantExport.tsv.gz"))
classification <-
  read_csv(paste0(path, "/cosmic_raw_data/classification.csv"))




# Extract cosmic id, occurrence and occurrence site
du1 <- annovar_cosmic_data %>% 
  mutate(cosmic_id = str_match(INFO, "ID=(.*?);")[,2]) %>% 
  mutate(occurrence = str_match(INFO, ";OCCURENCE=(.*?)$")[,2]) %>% 
  mutate(occ_sum = sapply(str_extract_all(occurrence, "(\\d)"),
                          function(x) sum(as.numeric(x)))) %>% 
  mutate(cosmic_presence = "In Cosmic")
head(cosmic_data)

write_rds(cosmic_data, "cosmic_data.rds")
