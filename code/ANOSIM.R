library(tidyverse)
library(vegan)


# ANOSIM = ANalysis Of SIMilarity test - developed by Clarke in 1993  


# Ref. 1) GUSTA ME https://sites.google.com/site/mb3gustame/hypothesis-tests/anosim 
# 2) chap9.3 - Book: statistical analysis of microbiome data with R 

## principles: 
# Input to ANOSIM test is often ranked (dis)similarities. 
# The dimension reduction and visualization capacities of NMDS and hypothesis testing 
# offered by ANOSIM are complementary approaches in evaluating nonparametric multivariate data 

# one- and two- way ANOSIM 

# The ANOSIM R compares the mean of ranked dissimilarities between groups to the mean of ranked 
# dissimilarities within groups. R value close to 1 suggests dissimilarity between groups while R 
# close to 0 suggests an even distribution of high and low ranks within and between groups 

# permutation method to test significance 

# Key assumption: the ranges (ranked) dissimilarities within groups are equal, or at least very similar 





setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")
set.seed(19861015)


# 5 steps to implement ANOSIM 
# step 1. calculate dissimilarity matrix 

# step 2. calculate rank dissimilarities and assign a rank of 1 to be the smallest dissimilarity 

# step 3. calculate the mean among- and within-group rank dissimilarities 

# step 4. calculate test statistic R using the above formula R ... 

# step 5. test for significance 


raried_beta <- read_csv("processed_data/beta_dist_rarefied_wide.csv") %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>% 
        as.matrix()


meta <- read_csv("processed_data/metadata2.csv") %>% 
        select(sample_id = sample.id, diet = DIET) %>% 
        mutate(diet = factor(diet))  # metadata 



require(vegan)
set.seed(19861015)
anosim_bc <- anosim(raried_beta, grouping = meta$diet, permutations = 1000)  # r= 0.007, p=0.24 
summary(anosim_bc)

plot(anosim_bc) # visualize mean ranked dissimilarity between-group vs within-group 

























































