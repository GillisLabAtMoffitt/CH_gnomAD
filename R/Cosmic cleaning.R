# Cosmic

cosmic_data <-
  read.delim(paste0(path, "/cosmic_raw_data/hg38_cosmic70.txt")) %>% 
  `colnames<-`(c("#Chr","Start","End","Ref","Alt","INFO"))
head(cosmic_data)

