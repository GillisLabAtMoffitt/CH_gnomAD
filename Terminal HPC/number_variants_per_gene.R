## Import library
library(tidyverse)


###################################################################### I ### List files
file_list <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_skweness/potential_CH_variants",
                        pattern = "*.vcf.gz",
                        recursive=FALSE,
                        full.names = TRUE)

known_CH_genes <-
  read_table("/share/lab_gillis/Christelle/gnomAD_skweness/data/3390851_Covered.bed", 
             col_names = FALSE, skip = 2)


###################################################################### II ### Plots
gnomad <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                 header = TRUE, sep = " "))

known_CH_gene_list <- paste0(c(unique(known_CH_genes$X4)), collapse = "|")

pdf("Number of variants per gene.pdf")
gnomad %>% 
  filter(!is.na(Gene.name)) %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  
  mutate(known_CH_gene = case_when(
    str_detect(Gene.name, known_CH_gene_list)       ~ "known CH gene",
    TRUE                                            ~ Gene.name
  )) %>% select(IDs, Gene.name, known_CH_gene) %>% 
  arrange(known_CH_gene) %>% 
  group_by(known_CH_gene) %>% 
  mutate(number_of_variants_per_gene = n()) %>% 
  distinct(Gene.name, .keep_all = TRUE) %>% 
  select(-IDs) %>% 
  
  group_by(known_CH_gene) %>% 
  mutate(number_of_gene = n()) %>%
  ungroup() %>% 
  mutate(average_of_variants_per_gene = number_of_variants_per_gene / number_of_gene) %>%
  distinct(known_CH_gene, .keep_all = TRUE) %>% 
  
  mutate(known_CH_gene = fct_relevel(known_CH_gene, "known CH gene")) %>%
  ggplot(aes(x= known_CH_gene, y= average_of_variants_per_gene,
             fill= known_CH_gene == "known CH gene"
  ))+
  geom_bar(stat = "identity", show.legend = FALSE)+
  ggtitle("Number of overall variants per gene")+
  labs(x= NULL)+
  coord_flip()
dev.off()
