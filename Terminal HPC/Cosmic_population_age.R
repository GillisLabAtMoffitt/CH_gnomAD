## Import library
library(tidyverse)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants",
                        pattern = "*.vcf.gz",
                        recursive=FALSE,
                        full.names = TRUE)


###################################################################### II ### Plots
gnomad <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " "))

pdf("Cosmic population age.pdf")
gnomad %>% 
  ggplot(aes(x= Age))+
  geom_histogram(binwidth = 1, fill= "purple")+
  xlim(0, 100)
dev.off()
