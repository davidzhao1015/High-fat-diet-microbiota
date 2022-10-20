library(tidyverse)
library(vegan)
library(GUniFrac)
library(ape)
library(ade4)



# Compare microbiome communities using the GUniFrac Package 

# UniFrac (Lozupone and Knight 2005), Weighted UniFrac (Lozupone et al. 2007) 
# and General UniFrac Distance Metrics (Chen et al. 2012)

# Pros and cons: 
## 1. UniFrac - considers only taxa presence and absence information; it most efficient in detecting abundance change 
## in rare lineages; completely ignore the taxa abundance information 
## 2. Weighted UniFrac - add a proportional weighting to UniFrac method; detect both changes in how many 
## sequences from each lineage are present, as well as in which taxa are present; it most sensitive to detect change 
## in abundant lineages; either UniFrac or weighted UniFrac is not powerful to detect change in moderately abundant 
## lineages 
## 3. Generalized UniFrac: contains extra parameter alpha controlling the weight on abundant lineage 
## so that the distance is not dominated by highly abundant lineages. 


# GUniFrac package is dependent on vegan and ape packages and the author also suggests using package ade4. 





# Example code - 
## Ref. 1) Chap 9.5 - Book: Statistical analysis of microbiome data with R 
## 2) Help document of GUniFrac package 

data("throat.otu.tab")
data("throat.tree")
data("throat.meta")

groups <- throat.meta$SmokingStatus  


# Rarefaction - GUniFrac is sensitive to different sequencing depth, to compare microbiomes 
# on an equal basis, rarefaction might be used. 
otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff 

# # alternative script to implement rarefaction using vegan package 
# otu.tab.rff_vegan <- rrarefy(throat.otu.tab, min(apply(throat.otu.tab, 1, sum)))    




# Calculate the UniFracs - GUniFrac requires a rooted phylogeny tree 
unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs   

str(unifracs)

dw <- unifracs[, , "d_1"] # weighted UniFrac 
du <- unifracs[, , "d_UW"] # unweighted UniFrac 
dv <- unifracs[, , "d_VAW"] # variance adjusted weighted UniFrac 
d0 <- unifracs[, , "d_0"] # GUniFrac with alpha 0
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5 


# Permanova - Distance based multivariate analysis of variance 
set.seed(19861015)

# merge d5 matrix and metadata to make sample id consistent 
d5_df <- d5 %>% 
        as_tibble(rownames = NA) %>% 
        rownames_to_column("sample") 


meta_df <- throat.meta %>% 
        as_tibble(rownames = NA) %>% 
        rownames_to_column("sample") 

d5_meta <- inner_join(d5_df, meta_df, by = "sample")

d5_dist <- d5_meta %>%
        select(2:61) %>% 
        as.dist()

ado <- adonis2(d5_dist ~ SmokingStatus,
               data = d5_meta, 
               permutations = 999)

ado # PERMANOVA test shows that smoking status explain 4% variance of generalized UniFrac (p=0.002) 



































