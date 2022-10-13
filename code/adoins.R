library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)
library(ggtext)
library(broom)



# Permutation Multivariate Analysis of Variance (PERMANOVA) with adonis() [adonis2() may be the upgraded version?] 
# in vegan package 

# Ref. 
# 1) minimalR by Dr. Patrick Schloss
# 2) Paper: Application of multivariate statistical techniques in microbial ecology 
# 3) Chap 9 - Book: Statistical analysis of microbiome data with R  



# Objective of codes: test significance of difference in Bray-Curtis distance between two diet groups; 
# in addition, test significance of difference in Bray-Curtis distance by diet and sex. 





# Ensure analysis is reproducible
set.seed(19861015)

# Set permutation parameter 
permutation <- 1000 

# Change working directory 
setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/") 




# Read in rarefied Bray-Curtis distance matrix 
bc_dist_rare <- read_csv("processed_data/beta_dist_rarefied_long.csv") %>%
        select(-1) # drop the first column (index)

# Read in the meta data 
meta <- read_csv("processed_data/metadata2.csv")  %>% 
        rename(sample_id = sample.id) %>% 
        select(sample_id, diet = DIET, sex = SEX) 



# join distance matrix and metadata 
dist_meta <- inner_join(bc_dist_rare, 
                        meta, 
                        by="sample_id") %>% 
        mutate(diet_sex = interaction(diet, sex))  # add a new interaction grouping diet * sex 
        
        

# convert to dist data 
all_dist <- dist_meta %>%
        select(all_of(.[["sample_id"]])) %>% 
        as.dist() 




# test significance for diet alone 
set.seed(19861015) 
adonis_diet <- adonis(all_dist~diet,
                      data = dist_meta,
                      permutations = permutation)

p.val_diet <- adonis_diet[["aov.tab"]][["Pr(>F)"]][1]  # significant different at 0.05 level  


# sequential test for [diet x sex] OR [sex X diet]  
set.seed(19861015)
adonis_diet_sex <- adonis(all_dist~diet*sex, data = dist_meta, permutations = 1000) 
adonis_diet_sex$aov.tab  # view output 
# diet explains 8.27% variance with p=0.001 while sex explains 8.54% variance 
# diet and sex have interaction effect on community structure with p=0.001 

set.seed(190861015)
adonis_sex_diet <- adonis(all_dist~sex*diet, data = dist_meta, permutations = 1000) 
adonis_sex_diet$aov.tab # result seems similar to the above although orders of two explanatory variables matters 



# test for newly-created four groups diet_sex 
set.seed(19861015)
adonis_group4 <- adonis(all_dist~diet_sex, data = dist_meta, permutations = 1000) 
adonis_group4$aov.tab
# based on Bray-Curtis distance, the group explains about 26.6% variance with significance at p=0.001 


# change the default dummy contrasts in R to the default "Sum" contrasts in vegan package 
set.seed(19861015) 
adonis_group4_contr.sum <- adonis(all_dist~diet_sex, 
                                  data = dist_meta, 
                                  permutations = 1000,
                                  contr.unordered = "contr.sum")

adonis_group4_contr.sum$aov.tab  # groups explain about 26.6% with p=0.001 


# add the contrasts for ordered factors 
set.seed(19861015)
adonis_group4_contr.poly <- adonis(all_dist~diet_sex, 
                                   data = dist_meta, 
                                   permutations = 1000,
                                   contr.unordered = "contr.sum",  # not quite sure this argument ... 
                                   contr.ordered = "contr.poly") # not quite sure this argument ... 

adonis_group4_contr.poly$aov.tab # groups explain about 26.6% with p=0.001  








