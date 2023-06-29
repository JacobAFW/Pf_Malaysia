# Packages
library(tidyverse)
library(janitor)
library(data.table)

# Metadata

Pf7_meta <- read_tsv("/g/data/pq84/malaria/Pf_Malaysia/data/malariaGEN/Pf7/Pf7_analysis_metadata.tsv") %>% 
  rename(Region = "Admin level 1") %>% 
  select(1:4, Population, Year) %>% 
  add_column(partner_sample_id = NA)

anstey_meta <- read.table("/g/data/pq84/malaria/Pf_Malaysia/data/malariaGEN/Pf.7.1/Pf_71__1202-PF-MY-ANSTEY__sample_metadata.txt", sep ="\t") %>% 
  row_to_names(1) %>% 
  rbind(
    read.table("/g/data/pq84/malaria/Pf_Malaysia/data/malariaGEN/Pf.7.2/Pf_72__1202-PF-MY-ANSTEY__sample_metadata.txt", sep ="\t") %>%
      row_to_names(1) %>%
      mutate(partner_sample_id = Sample_External_ID)
  ) %>% 
  select(1,2,4) %>% 
  add_column(Country = "Malaysia") %>% 
  add_column(Region = "Sabah") %>% 
  add_column(Population = NA) %>% 
  add_column(Year = "2018") %>%
  rename(Sample = sample) %>% 
  rename(Study = study) 

combined_meta <- Pf7_meta %>%
  rbind(anstey_meta) 

Plink_fam <- read.table("Pf.fam", sep =" ") %>% 
  rename(Sample = V1) %>% 
  left_join(combined_meta) 

write.table(Plink_fam, "Pf.fam", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
Plink_fam %>% 
  select(-c(2:6)) %>% 
  write_csv(Plink_fam, "Pf.csv")
