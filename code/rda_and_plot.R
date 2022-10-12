library(tidyverse) 
library(vegan)
library(GUniFrac)


# Redundancy Analysis (RDA) + partial RDA  

# RDA extract and summarize variation in a set of response variables that can be explained by a set of explanatory
# variables, namely a direct gradient analysis. 
# RDA is considered as a constrained version of PCA 

# Ref: 
# 1) Book: Statistical analysis of microbiome data with R; 
# 2) https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html 
# 3) partial RDA https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html 



setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")

set.seed(19861015)  


# step 1. load example data for the following analysis 
data("throat.otu.tab") # load otu count data from GUniFrac package 
data("throat.meta") # load meta data 
 
throat_meta <- throat.meta %>% 
        select(SmokingStatus, Age, Sex, PackYears) 



 
# step 2. transform otu count data with Hellinger method 
# standardize response variables and explanatory variables if they are not dimentionally homogeneous 
# turn qualitative variable to dummy variables 
# examine distribution of response and explanatory variables, and apply transformation when needed 
# non-Euclidean distance matrix (eg. Hellinger distance) is applicable to RDA  

otu_hell <- decostand(throat.otu.tab, "hell")   

# standardize quantitative explanatory data 
throat_meta$Age <- decostand(throat_meta$Age, method = "standardize")

throat_meta$PackYears <- decostand(throat_meta$PackYears, method = "standardize") 




# step 3. implement RDA between otu and the explanatory variable, diet 
rda_out_hell <- rda(otu_hell ~ ., throat_meta)  

summary(rda_out_hell)
rda_out_hell  # print output  - constrained axis explain about 4.8% variation of community data 




# step 4. forward selection to reduce the number of variables entering the analysis 
# two alternative functions of vegan can do this job, including ordistep() and ordiR2step 

step_forward <- ordistep(rda(otu_hell ~ 1, data=throat_meta),
                         scope = formula(rda_out_hell),
                         direction = "forward",
                         pstep=1000)  # SmokingStatus, Sex 

# alternative method 
step_forwad2 <- ordiR2step(rda(otu_hell ~1, data=throat_meta),
                           scope = formula(rda_out_hell),
                           direction = "forward",
                           pstep = 1000) # SmokingStatus





# step 5. final model including PackYear and sex 
rda_final <- rda(otu_hell ~ SmokingStatus + Sex, data = throat_meta) 

# what is the model's explanatory power?  
RsquareAdj(rda_final)$adj.r.squared  # adjusted R^2 = 4.4% 
# adjusted R^2 measures the strength of the relationship between response and explanatory variables, but applies 
# a correction of the R^2 to take into account the number of explanatory variables. This is the statistic that 
# should be reported. 

# is the model statistically significant? Yes, p=0.001 
anova.cca(rda_final, step =1000)  # significance of global model 

anova.cca(rda_final, by="axis", step = 1000)  # sig. for axis 

anova.cca(rda_final, by = "term", step = 1000) # sig. for explanatory variables  





# partial RDA (pRDA): remove the effect of one or more explanatory variables on a set of response variables 
# prior to a standard RDA. This may be useful when well-characterized variables with strong effects obscure 
# the effects of more interesting explanatory variables 

## focus on effect of SmokingStatus controlling for rest covariables including Sex, Age and PackYears 

# method 1: subset explanatory (environmental) variables into SmokingStatus and covariables 
env.smoke <- subset(throat_meta, select=c(SmokingStatus)) 
env.covariables <- subset(throat_meta, select = c(-SmokingStatus)) 

pRDA <- rda(X=otu_hell, Y=env.smoke, Z=env.covariables)
summary(pRDA) 
# interpretation: smoking status (constrained proportion) explains 2.6% of the variation in bacterial community 
# composition across fecal samples, while covariables (conditioned proportion) explains 7.8% of this variation. In
# addition, unexplained variation in bacterial community across fecal samples reaches 89.6%. 

RsquareAdj(pRDA)$adj.r.squared  # explanatory power is 1.1% 

anova.cca(pRDA, step = 1000)  # the model is statistically significant with p = 0.05  

# method 2: alternative syntax 
pRDA2 <- rda(otu_hell~ SmokingStatus +  # effect of interest 
                     Condition(Sex + Age + PackYears),  # effect to control 
             data = throat_meta) 

summary(pRDA2) 





# step 6. look up canonical coefficients (equivalent of regression coefficients) of 
# each explanatory variable on each canonical axis 
coef(rda_final) 



# step 7. apply Kaiser-Guttman criterion to residual axes 
rda_final$CA$eig >= mean(rda_final$CA$eig)



# step 8. plot the results  of RDA 
## default/ quick plot 
ordiplot(rda_final, scaling = 1, type="text")  # distances among objects reflect their similarity 
ordiplot(rda_final, scaling = 2, type = "points") # angles between variables reflect their correlation 


# customize RDA plot with ggplot2 
## extract % explained by the first 2 axes 
perc <- round(100*(summary(rda_final)$cont$importance[2, 1:2]),2)   

## extract scores - these are coordination in the RDA space 
sc_samples <- scores(rda_final, display = "sites", choices = c(1,2), scaling =1) 

sc_otu <- scores(rda_final, display = "species", choices = c(1,2), scaling =1) 

sc_exp <- scores(rda_final, display = "bp", choices = c(1,2), scaling = 1)  


plot(rda_final,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_samples, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_otu, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_otu + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_otu), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the explanatory variables
arrows(0,0, # start them from (0,0)
       sc_exp[,1], sc_exp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_exp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_exp[,2] - 0.03, 
     labels = rownames(sc_exp), 
     col = "red", 
     cex = 1, 
     font = 2)





















