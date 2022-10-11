library(tidyverse) 
library(vegan)
library(GUniFrac)


# redundancy analysis (RDA) 
# RDA extract and summarize variation in a set of response variables that can be explained by a set of explanatory
# variables, namely a direct gradient analysis. 
# RDA is considered as a constrained version of PCA 

setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")

set.seed(19861015)



# pre-analysis 
# standardize response variables and explanatory variables if they are not dimentionally homogeneous 
# turn qualitative variable to dummy variables 
# examine distribution of response and explanatory variables, and apply transformation when needed 
# non-Euclidean distance matrix (eg. Hellinger distance) is applicable to RDA  


# total variance = constrained and unconstrained variance  
# scores of response and explanatory variables 
# examine non-canonical (unconstrained) vectors of RDA solution by ordination and correlation, resulting in 
# insights into the behavior of these residuals 
# scaling type I and II plots 


$# Ref 1) Book, Statistical analysis of microbiome data with R; 
# 2) https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html 





# step 0. load example data for the following analysis 
data("throat.otu.tab") # load otu count data from GUniFrac package 
data("throat.meta") # load meta data 
 
throat_meta <- throat.meta %>% 
        select(SmokingStatus, Age, Sex, PackYears) 

 
# step 1. transform otu count data with Hellinger method 
otu_hell <- decostand(throat.otu.tab, "hell")   



# step 2. implement RDA between otu and the explanatory variable, diet 
rda_out_hell <- rda(otu_hell ~ ., throat_meta)  

summary(rda_out_hell)
rda_out_hell  # print output  - constrained axis explain about 4.8% variation of community data 



# step 3. look up canonical coefficients (equivalent of regression coefficients) of 
# each explanatory variable on each canonical axis 
coef(rda_out_hell) 



# step 4. look up R^2 and adj R^2 
RsquareAdj(rda_out_hell)  
# adj R^2 is negative meaning that explanatory variables explain less variation than 
# the same number of randomly generate variables 



# step 5. apply Kaiser-Guttman criterion to residual axes 
rda_out_hell$CA$eig >= mean(rda_out_hell$CA$eig)


# step 6. plot the results  of RDA 
ordiplot(rda_out_hell, scaling = 1, type="text")  # distances among objects reflect their similarity 
ordiplot(rda_out_hell, scaling = 2, type = "text") # angles between variables reflect their correlation 

# customize RDA plot with ggplot2 
## extract % explained by the first 2 axes 
perc <- round(100*(summary(rda_out_hell)$cont$importance[2, 1:2]),2)   

## extract scores - these are coordination in the RDA space 
sc_samples <- scores(rda_out_hell, display = "sites", choices = c(1,2), scaling =1) 

sc_otu <- scores(rda_out_hell, display = "species", choices = c(1,2), scaling =1) 

sc_exp <- scores(rda_out_hell, display = "bp", choices = c(1,2), scaling = 1)  


plot(rda_out_hell,
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




# step 7. permutation to test significance of variations in response variables explained by exploratory variables 
anova.cca(rda_out_hell, step = 1000)  # permutation test on the global model 
anova.cca(rda_out_hell, by="axis", step=1000) # test on axis 
anova.cca(rda_out_hell, by = "terms", step =1000) # test on explanatory variable  


# step 8. forward selection to reduce the number of variables entering the analysis 
# two alternative functions of vegan can do this job, including ordistep() and ordiR2step 

step_forward <- ordistep(rda(otu_hell ~ 1, data=throat_meta),
                         scope = formula(rda_out_hell),
                         direction = "forward",
                         pstep=1000)   # select Sex and PackYears  



# alternative method 
step_forwad2 <- ordiR2step(rda(otu_hell ~1, data=throat_meta),
                           scope = formula(rda_out_hell),
                           direction = "forward",
                           pstep = 1000) 



# step 9. final model including PackYear and sex 
rda_final <- rda(otu_hell ~ SmokingStatus + Sex, data = throat_meta) 
RsquareAdj(rda_final)


anova(rda_final, step =1000)  # significance of global model 

anova(rda_final, by="axis", step = 1000)  # sig. for axis 

anova(rda_final, by = "term", step = 1000) # sig. for explanatory variables  






















