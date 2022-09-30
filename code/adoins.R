library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)
library(ggtext)



# Ensure analysis is reproducible
set.seed(19861015)

# Set permutation parameter 
permutation <- 1000 

# Change working directory 
setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/") 




# Read in Bray-Curtis distance matrix with rarefaction 
bc_dist_rare <- read_csv("processed_data/beta_dist_rarefied_long.csv") %>%
        select(-1) # drop the first column (index)

# Read in the meta data 
meta <- read_csv("processed_data/metadata2.csv")  %>% 
        rename(sample_id = sample.id)

# join distance matrix and metadata 
dist_meta <- inner_join(bc_dist_rare, 
                        meta, 
                        by="sample_id")   

# convert to dist data 
all_dist <- dist_meta %>%
        select(all_of(.[["sample_id"]])) %>% 
        as.dist() 



# Hypothesis testing to compare centers by two diet groups using vegan::adonis   
adonis_test <- adonis(all_dist~DIET,
                      data = dist_meta,
                      permutations = permutation)

p.val_diet <- adonis_test[["aov.tab"]][["Pr(>F)"]][1]  # significant different at 0.05 level  
















