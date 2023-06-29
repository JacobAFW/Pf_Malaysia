# Converting VCF to GDS for estMOI
library(tidyverse)
library(SeqArray)
library(moimix)

seqVCF2GDS("moimix/Pf3D7_moimix_subset_2.vcf.gz", "moimix/Pf3D7_moimix_subset_2.gds", storage.option = "LZ4_RA")
isolates <- seqOpen("moimix/Pf3D7_moimix_subset_2.gds") # read into R
seqSummary(isolates)
sample.id <- seqGetData(isolates, "sample.id")

coords <- getCoordinates(isolates) # get genomic coordinates of all variants

## estimating BAF matrix 
isolate_baf <- bafMatrix(isolates)

baf_df <- isolate_baf$coords %>%
  cbind(
    as.data.frame(isolate_baf$baf_site) %>% 
      rename(baf = "isolate_baf$baf_site")
    ) %>%
  left_join(
    isolate_baf$baf_matrix %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("variant.id") %>%
      mutate(variant.id = as.numeric(variant.id)) 
    ) 

write_tsv(baf_df, "moimix/BAF_dataframe_2.tsv")


# Estimate Fws
set.seed(2023)

fws_all <- getFws(isolates) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  rename("Proportion" = ".")

## Plot MOI distribution and export ######################

## output MOI data
fws_all %>%  
  write_tsv("moimix/fws_MOI_2.tsv")

## output samples with HIGH MOI data
fws_all %>%  
  filter(Proportion < 0.85) %>% 
  select(sample) %>%
  write_tsv("moimix/fws_high_MOI_2.tsv")