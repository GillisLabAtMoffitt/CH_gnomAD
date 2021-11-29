# Import library
library(tidyverse)

# Load data
path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
cosmic_data <-
  read.delim(paste0(path, "/cosmic_raw_data/hg38_cosmic70.txt")) %>% 
  `colnames<-`(c("#Chr","Start","End","Ref","Alt","INFO"))
head(cosmic_data)


# Extract cosmic id, occurrence and occurrence site
cosmic_data1 <- cosmic_data %>% 
  mutate(cosmic_id = str_match(INFO, "ID=(.*?);")[,2]) %>% 
  mutate(occurrence = str_match(INFO, ";OCCURENCE=(.*?)$")[,2]) %>% 
  mutate(occ_sum = sapply(str_extract_all(cosmic_data1$occurrence, "(\\d)"),
                          function(x) sum(as.numeric(x))))
head(cosmic_data1)
