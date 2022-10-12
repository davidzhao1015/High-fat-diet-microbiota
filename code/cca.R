library(tidyverse) 
library(vegan)
library(GUniFrac)



# Canonical correspondence analysis (CCA) - constrained ordination 

# Ref: 
# 1) GUSTA ME https://sites.google.com/site/mb3gustame/constrained-analyses/cca 
# 2) Chap7.4.6 Book: Statistical analysis of microbiome data with R 

# Assumptions: 
# 1) response variable show unimodal distribution across objects. A sampling gradient must be long enough to allow 
# the increase and decrease of a given species or OTU across the sites sampled. Gradients that are too short may 
# manifest linear responses and may be netter handled by RDA, although CCA may also handle linear relationships. 
# 2) Explanatory variables show linear, causal relationships to the response data. If one is unsure if their is 
# causal relationship between an explanatory variable and the response data, interpreting should performed with care. 



setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")  

set.seed(19861015)  


# step 1. load example data for the following analysis 
data("throat.otu.tab") # load otu count data from GUniFrac package 
data("throat.meta") # load meta data 

throat_meta <- throat.meta %>% 
        select(SmokingStatus, Age, Sex, PackYears)  




# step 2. implement CCA. please do not transform OTU data 
smoker_cca <- cca(throat.otu.tab ~ ., throat_meta) 

smoker_cca 

summary(smoker_cca)
# CCA axis (constrained proportions) explains 7.7% of variance in bacterial community across samples 
# while 92.3% variance can not be explained. 




# step 3. test significance of the global model 
anova.cca(smoker_cca, step = 1000)  # p=0.099 not statistically significant! 




# step 4. select variables with forward selection method 
ordistep(cca(throat.otu.tab ~ 1,
             data = throat_meta),
         scope = formula(smoker_cca),
         direction = "forward",
         pstep = 1000) # SmokingStatus selected 



# step 5 . build final model 
final_cca <-  cca(throat.otu.tab ~ SmokingStatus, throat_meta) 

anova.cca(final_cca,  step =1000)  # p=0.003, statistically significant! however, the final model paid off with 
# larger residual compared to the full model 


# step 6. quick plot final model of CCA 

plot(smoker_cca,
     scaling = 1, 
     display = c("lc", "cn"),
     main = "Biplot CCA-scaling 1")    























