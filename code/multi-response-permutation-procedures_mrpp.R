library(tidyverse)
library(vegan)



# MRPP = multi-response permutation procedures 

# MRPP tests whether there is a significant difference between two or more groups of sampling units 
# MRPP is nonparametric method. The testing difference may be differences in mean (location) or differences 
# in within-group distance (spread) 

# A = 1 - obs-delta/exp-delta 

# effect size A = 0 : within-group heterogeneity equals expectation by chance; effect size A =1: all items 
# are identical within groups. A > 0.3 is fairly high in ecology. However, based on our knowledge, there is 
# no reporting criterion of effect size in microbiome literature. 





# implement MRPP using mrpp() of vegan
## eligible input data to mrpp() function include  1)raw OTU matrix, 2)a dissimilarity matrix, 
# and 3)a symmetric square matrix  
setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")
set.seed(19861015)


# read rarefied bray-curtis distance matrix 
raried_beta <- read_csv("processed_data/beta_dist_rarefied_wide.csv") %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>% 
        as.matrix()

# read metadata 
meta <- read_csv("processed_data/metadata2.csv") %>% 
        select(sample_id = sample.id, diet = DIET) %>% 
        mutate(diet = factor(diet))  

# raw otu count data 
otu_count <- read_csv("processed_data/otutable_wide.csv") %>% 
        select(-1) %>% 
        as_tibble() %>% 
        inner_join(meta, by="sample_id") %>% 
        column_to_rownames("sample_id")

diet <- otu_count$diet

otu_count2 <- otu_count %>%
        select(-diet)  


# mrpp with rarefied bray-curtis matrix 
mrpp(raried_beta, diet, permutations=1000)  # A=0.0496 with p=0.000999  

# alternative input: raw otu count matrix 
mrpp(otu_count2, diet, permutations = 1000, distance = "bray") # A=0.04937 with p=0.000999  




# obtain mean distance matrix using meandist() 
## the function meandist() calculates a matrix of mean within-cluster (block) dissimilarities (diagonal)
# and between-cluster dissimilarities (off-diagonal elements), and an attribute n of grouping counts 

raried_beta2 <- raried_beta %>% as.dist()  # convert to a dist object 
        
bray_mrpp <- meandist(dist =raried_beta2, diet, permutations=1000)   
summary(bray_mrpp)


















