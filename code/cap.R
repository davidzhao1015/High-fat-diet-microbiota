library(tidyverse) 
library(vegan)
library(GUniFrac) 


# Constrained analysis of principal coordinates (CAP)

# CAP = constrained analysis of proximity (in vegan package) 
# CAP allows non-Euclidean dissimilarity distance, including Manhattan or Bray-Curtis distance 

# Ref: Chap 7 - Book: Statistical analysis of microbiome data with R 




setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")  

set.seed(19861015)  



# step 1. load example data for the following analysis 
data("throat.otu.tab") # load otu count data from GUniFrac package 
data("throat.meta") # load meta data 

throat_meta <- throat.meta %>% 
        select(SmokingStatus, Age, Sex, PackYears) 




# step 2. implement CAP with controlling for Age on Bray-Curtis distance 
throt_cap <- capscale(throat.otu.tab~SmokingStatus + Sex + PackYears + Condition(Age), 
                      throat_meta,
                      dist= "bray") 

throt_cap # smoking status, sex and pack years explain 8.83% of total variation  




# step 3.  test on significance of the global model 
set.seed(19861015)
anova(throt_cap, step = 1000) # permutation test on significance of the global model p=0.006 

anova(throt_cap, by = "axis", step = 1000) # sig. for axis p=0.007 

anova(throt_cap, by="term", step = 1000) # sig. for env variables smoking status and sex stat sig! 




# step 4. forward selection to reduce environmental variables 
step_forward <- ordistep(capscale(throat.otu.tab~1,
                                  data = throat_meta),
                         scope = formula(throt_cap),
                         direction = "forward", 
                         pstep = 1000) # SmokingStatus is selected 


cap_final <- capscale(throat.otu.tab~SmokingStatus, throat_meta, dist = "bray")

cap_final 

anova(cap_final, step = 1000) # the final model is statistical significant p=0.003 



# step 5. quick plot for CAP 
plot(throt_cap) 


## customize CAP plot 
groups <- throat.meta$SmokingStatus 
groups 

plot(cap_final, type ="n") # create an empty CAP ordination diagram 
points(cap_final, col = as.numeric(as.factor(groups)),
       pch = as.numeric(as.factor(groups))) # add points 
ordispider(cap_final, groups, lty=2, col="grey", label = T) # spider plot to connect individual members of the group
# with the group centroid 
ordiellipse(cap_final, groups, lty=2, col="grey", label =F) # encircles clouds of points within the group 
# by ellipse like the envelopes 




        
        
        
        
        
        
        
        
        
        
        








