# Load packages
library(tidyverse)
library(data.table)

# Sample-based missing data
## Read in data
sample_mis <- read_table("Pf.imiss", col_names=T) 

# Density of miss by sample
sample_miss_density <- sample_mis %>% 
    rename(Sample = IID,
        Geno_Miss = F_MISS) %>%
    select(2:ncol(.)) %>%
    ggplot(aes(x = Geno_Miss)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept = 0.25, colour = "#1F968BFF") 

ggsave("filtering/sample_miss_density_plot.png", width = 12, dpi = 600, sample_miss_density)

# Malaysia mean missingness
sample_mis %>%
    left_join(
        read_csv("Pf.csv") %>% 
            rename(FID = Sample)
    ) %>%
    filter(Country == 'Malaysia') %>%
    summarise(mean = mean(F_MISS), sd = sd(F_MISS), min = min(F_MISS), max = max(F_MISS)) %>% 
    write_csv("filtering/malaysia_missingness.csv")


# Variant-based missing data
## Read in data
var_miss <- read_table("Pf.lmiss", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*")) %>%
    rename(Geno_Miss = F_MISS,
            Chr = CHR)

# Density of miss by variant
var_miss_density <- var_miss %>% 
    ggplot(aes(x = Geno_Miss)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept = 0.25, colour = "#1F968BFF") +
    geom_vline(xintercept = 0.1, colour = "#1F968BFF") 

ggsave("filtering/variant_miss_density_plot.png", width = 12, dpi = 600, var_miss_density)


# Summary stats

miss_summary <- var_miss %>%
    summarise(Miss_mean = mean(Geno_Miss),
                Miss_SD = sd(Geno_Miss),
                Miss_SE = (sd(Geno_Miss))/sqrt(n()))

write_csv(miss_summary, "filtering/pf_miss_summary.csv")