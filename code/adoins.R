library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)
library(ggtext)


set.seed(1015)
permutation <- 1000 


# read in bray-curtis distance matrix 

dist <- read_csv("./processed_data/bc_dist_mat.csv",
                 skip=1,
                 col_names = FALSE)  

colnames(dist) <- c("sample.id", dist$X1)

# read in the meta-data 
meta <- read_csv("./processed_data/metadata2.csv")

dist_meta <- inner_join(dist, meta, by="sample.id")  

all_dist <- dist_meta %>%
        select(all_of(.[["sample.id"]])) %>% 
        as.dist() 


# test center variation between two groups with adonis 
adonis_test <- adonis(all_dist~DIET,
                      data = dist_meta,
                      permutations = permutation)


p.val_diet <- adonis_test[["aov.tab"]][["Pr(>F)"]][1]  # significant different at 0.05 level  
















