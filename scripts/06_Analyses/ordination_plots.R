# Load packages
library(tidyverse)
library(ape)
library(ggtree)
library(viridis)

# PCA

## Read in data
pca <- read_table("Pf.eigenvec", col_names=F) %>%
    select(-1) %>%
    rename("sampleid" = "X2") %>%
    rename_at(vars(starts_with("X")), ~str_replace(., "X.*", paste0("PC", seq_along(.))))

eigenval <- scan("Pf.eigenval")

## Add metadata 
pca <- pca %>%
    rename(Sample = sampleid) %>% 
    left_join(read_csv("Pf.csv")) %>% 
    mutate(Population = ifelse(Study == "1202-PF-MY-ANSTEY", "Malaysian_outbreak", .$Population)) %>% 
    mutate(Malaysian_outbreak = ifelse(Population == "Malaysian_outbreak", "Yes", "No")) %>% 
    mutate(Malaysian_outbreak = fct_relevel(Malaysian_outbreak, "Yes"))


## Percentage variance explained by each PC
pve_plot <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100) %>%
    ggplot(aes(PC, pve)) + 
    geom_bar(stat = "identity") + 
    ylab("Percentage variance explained") + 
    theme_light()

ggsave("clustering/Percentage_variance_explained.png", dpi=600, pve_plot)


## PCA - clusters
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Population, shape = Malaysian_outbreak)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("clustering/PCA_population.png", dpi = 600, pca_plot)


############################################################################################################################################################

# MDS - clusters

## Read in data
mds <- read_table("Pf.mds", col_names=T) %>%
    select(-c("FID", "X14")) %>%
    rename("sampleid" = "IID") %>%
    rename_at(vars(starts_with("C")), ~str_replace(., "C", "MDS"))

## Add metadata 
mds <- mds %>%
    rename(Sample = sampleid) %>% 
    left_join(read_csv("Pf.csv")) %>% 
    mutate(Population = ifelse(Study == "1202-PF-MY-ANSTEY", "Malaysian_outbreak", .$Population)) %>% 
    mutate(Malaysian_outbreak = ifelse(Population == "Malaysian_outbreak", "Yes", "No")) %>% 
    mutate(Malaysian_outbreak = fct_relevel(Malaysian_outbreak, "Yes"))

## Plot MDS
mds_plot <- ggplot(mds, aes(MDS1, MDS2, colour = Population, shape = Malaysian_outbreak)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light()

ggsave("clustering/MDS_population.png", dpi=600, mds_plot)

############################################################################################################################################################


# Neighbour-joining tree

## Create distance matrix from PLINK data and build NJT
NJT_ID <- read_table("Pf.dist.id", col_names=F) %>%
    as.data.frame()

NJT_matrix <- read_table("Pf.dist", col_names=NJT_ID$X1) %>%
    as.data.frame() %>%
    add_column(Row_Names = NJT_ID$X1) %>%
    column_to_rownames("Row_Names") %>%
    as.matrix()

NJT_tree <- nj(NJT_matrix)

## Plot tree

# Read in metadata
NJT_metadata <- read_csv("Pf.csv") %>% 
    mutate(Population = ifelse(Study == "1202-PF-MY-ANSTEY", "Malaysian_outbreak", .$Population)) %>% 
    mutate(Malaysian_outbreak = ifelse(Population == "Malaysian_outbreak", "Yes", "No")) %>% 
    mutate(Malaysian_outbreak = fct_relevel(Malaysian_outbreak, "Yes"))

#Plot 
options(ignore.negative.edge=TRUE)

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = Population)) %<+% NJT_metadata +
    theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.key = element_blank()) +
    scale_color_viridis_d()
    
ggsave("clustering/NJT_tree_population_unrooted.png", dpi = 300, NJT_tree_plot)
