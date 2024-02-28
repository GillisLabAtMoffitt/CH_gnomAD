## Import library
library(tidyverse)


###################################################################### I ### List and read files
file_list <- list.files(
  path =
    "/share/lab_gillis/Christelle/cosmic_raw_data/splitted_data/filter_cosmic_withIDs",
  pattern = "*.vcf.gz",
  recursive=FALSE,
  full.names = TRUE)

cosmic_data <- do.call("rbind",lapply(Sys.glob(file_list), read.delim,
                                          header = TRUE, sep = " "))

# path <- fs::path("", "Volumes", "Lab_Gillis", "Christelle")
# cosmic_data <-
#   read.delim(paste0(path, "/cosmic_raw_data/splitted_data/filter_cosmic_withIDs/chr2_cosmic.vcf.gz"),
#               sep = " ")


###################################################################### II ### Analysis
print(paste0("The number row in cosmic is ", 
             nrow(cosmic_data)
             ))

print(paste0("The number single Sample.name in cosmic is ", 
             nrow(cosmic_data %>% distinct(Sample.name))
))

print(paste0("The number single ID_sample in cosmic is ", 
             nrow(cosmic_data %>% distinct(ID_sample))
))

print(paste0("The number single ID_tumour in cosmic is ", 
             nrow(cosmic_data %>% distinct(ID_tumour))
))

print("Cancers represented in Cosmic")
table(cosmic_data$Primary.site)

print("Somatic mutations in cosmic")
table(cosmic_data$Mutation.somatic.status)

print("Sample.Type in cosmic")
table(cosmic_data$Sample.Type)

print("Tumour.origin in cosmic")
table(cosmic_data$Tumour.origin)

a <- cosmic_data %>% distinct(gene_name_cosmic)
write_delim(a, "genes in cosmic data.txt")

