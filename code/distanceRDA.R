library(tidyverse) 
library(vegan) 



# distance-based RDA 


setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota")



# load example data for the following analysis 
data("throat.otu.tab") # load otu count data from GUniFrac package 
data("throat.meta") # load meta data 

throat_meta <- throat.meta %>% 
        select(SmokingStatus, Age, Sex, PackYears) 


db_rda <- capscale(throat.otu.tab ~ SmokingStatus + Age + Sex + PackYears, throat_meta, dist = "bray")



plot(db_rda)

anova(db_rda) # permutation test on significance of global model 

anova(db_rda, by="terms", permutations = 200) # permutation test on significance of env variables 


# extract scores for customized plots  
scores_db_rda <- scores(db_rda)  
scores_sample <- scores_db_rda$sites 
scores_otu <- scores_db_rda$species 













