# Load packages
library(tidyverse)

# Read in data
pf_maf <- read_table("Pf.frq", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*"))

# Density of MAF by SNP
pk_maf_density <- pk_maf  %>%
    ggplot(aes(x = MAF)) +
    geom_density(alpha=.3) 

ggsave("filtering/Pf_maf_density.png", dpi = 600, pk_maf_density)

# Plot MAF by chr
pf_maf_chr <- pf_maf %>% 
    mutate(SNP = str_remove(SNP, "Pf3D7_"),
        SNP = str_remove(SNP,"_v3.*")) %>% 
    group_by(SNP) %>%
    summarise(MAF = mean(MAF)) %>%
    ggplot(aes(x = SNP, y = MAF)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
    scale_fill_viridis_d() 

ggsave("filtering/Pf_maf_plot_chr.png", dpi = 600, pf_maf_chr)

# Summary stats

pf_maf %>% 
    mutate(SNP = str_remove(SNP, "Pf3D7_"),
        SNP = str_remove(SNP,"_v3.*")) %>% 
    group_by(SNP) %>%
    summarise(MAF_mean = mean(MAF),
            MAF_SD = sd(MAF),
            MAF_SE = (sd(MAF))/sqrt(n())) %>%
rbind(
    pk_maf %>%
        summarise(MAF_mean = mean(MAF),
                MAF_SD = sd(MAF),
                MAF_SE = (sd(MAF))/sqrt(n())) %>%
        add_column(SNP = "Total")
    ) %>%
    write_csv(col_names = T, "filtering/pf_maf_summary.csv")
