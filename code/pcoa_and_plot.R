library(tidyverse)
library(vegan)
library("BiodiversityR")  
require(Rcmdr)
require(rgl) 
require(multcomp)   

search()


# unconstrained PCoA (principal coordinates analysis) and bi-plot 

setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota")

set.seed(19861015)



# read in rarefied Bray-Curtis distance matrix 
bc_rare_matrix <- read_csv("processed_data/beta_dist_rarefied_long.csv") %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>%
        as.matrix()


# read in absolute abundance data of OTU 
otu_abs_abund <- read_csv("processed_data/otutable_wide.csv") %>% 
        select(-1) %>% 
        as_tibble() %>% 
        column_to_rownames("sample_id") 


# relative abundance data of OTU  
otu_rel_abund <- otu_abs_abund 
for(i in 1:nrow(otu_abs_abund)){
        otu_rel_abund[i, ] <- otu_abs_abund[i, ]/ sum(otu_abs_abund[i, ])
}



# write in metadata 
metadata <- read_csv("processed_data/metadata2.csv")  

metadata_diet <- metadata %>% 
        select(sample_id = sample.id, diet =DIET)



# implement PCoA on BC distance matrix with vegan package 
PCoA <- cmdscale(bc_rare_matrix, # input is bc distance matrix 
                 eig = TRUE, # which saves the eigenvalues 
                 k=2) # default values for the number of dimensions to return 

PCoA


# extract eigenvalues of first two axis 
# evaluate PCoA output using 1) Kaiser-Gutman criterion and 2) broken-stick model 
explainedvar1 <- round(PCoA$eig[1] / sum(PCoA$eig), 2) * 100  # % of total variation explained by PC1 
explainedvar2 <- round(PCoA$eig[2] / sum(PCoA$eig), 2) * 100 # % of total variation explained by PC2 

sum_eig <- sum(explainedvar1, explainedvar2) 
sum_eig # first two axis explain 38% variation of data set 

# scree plot 
barplot(PCoA$eig, 
        names = paste('PCoA', 1:length(PCoA$eig)),
        las= 3,
        ylab = "eigenvalues")



# Investigate performance of PCoA with Kaiser-Guttman criterion and broken-stick model criterion 

# Kaiser-Guttman criterion: the eigenvalue associated with first few axes should be larger 
# than the average of all the eigenvalues 
# the criterion of broken-stick model compares the eigenvalues associated with first few axis to 
# the expectations of the broken-stick model 
par(mar = c(5,5,1,2) + 0.1)
plot(PCoA$eig[1:20], xlab="PCoA", ylab = "Eigenvalue",
     las =1, cex.lab=1.5, pch = 16)  # plot eigenvalues 
abline(h=mean(PCoA$eig), lty=2, lwd=2, col = "blue")
# add expectation based on Kaiser-Guttman criterion and broken stick model  
b_stick <- bstick(length(PCoA$eig), sum(PCoA$eig)) 
lines(1:length(PCoA$eig), b_stick, type="l", lty=4, lwd=2, col="red") 
legend("topright", 
       legend = c("Average Eigenvalue", "Broken-Stick"),
       lty = c(2,4), bty="n", col= c("blue", "red"))  # add legend  



# ordination plot PCoA with ggplot2 
sample_score <- PCoA$points %>% 
        as_tibble(rownames = "sample_id") %>%  # convert to tittble 
        select(PC1 = V1, PC2 = V2, sample_id) %>% 
        inner_join(., metadata_diet, by= "sample_id") %>% 
        mutate(diet = factor(diet))



ggplot(sample_score, aes(x= PC1, y= PC2)) +
        geom_point()+
        labs(x = "PCoA1 (25%)", y = "PCoA2 (13%)") 



# which genera drive the observed divergence among points? 
require(BiodiversityR) # some dependent packages are missing so that some functions may not be available 

PCoA <- add.spec.scores(PCoA, otu_rel_abund, method="pcoa.scores") 



# correlation of each toxon along the PCoA axes 
taxa_corr <- add.spec.scores(PCoA, otu_rel_abund, method = "cor.scores")$cproj %>% 
        as_tibble(rownames = "taxa") %>% 
        select(PC1 = Dim1, PC2 = Dim2, taxa)
        
corrcut <- 0.55 # set threshold 

import_taxa <- taxa_corr %>% 
        filter(abs(PC1) >= corrcut | abs(PC2) >= corrcut)  %>% # import taxa influencing first two axis  
        separate(taxa, 
                 into=c("kingdom", "phylum", "class", "order", "family","genus","species","strain"), 
                 sep=";") %>% 
        select(PC1, PC2, species) 



# envfit() function from vegan package to conduct a permutation test for general abundances across axes 
# on these correlations 

fit <- envfit(PCoA, otu_rel_abund, perm = 999)       
 
otu_score <- fit$vectors$arrows %>% 
         as_tibble(rownames = "taxa") %>%
         select(taxa,  PC1 = Dim1, PC2 = Dim2)

r <- fit$vectors$r %>%
        as_tibble(rownames = "taxa") %>% 
        select(taxa, r= value)

pval <- fit$vectors$pvals %>% 
        as_tibble(rownames = "taxa") %>% 
        select(taxa, pvalue = value)

otu_score2 <- otu_score %>% 
        inner_join(r, by = "taxa") %>% 
        inner_join(pval, by = "taxa") %>% 
        filter(pvalue < 0.05 & r > 0.5) %>% 
        mutate(x=0, y=0, xend = PC1, yend = PC2) 


ggplot(sample_score, aes(x= PC1, y= PC2)) +
        geom_point(aes(color = diet), size = 3)+
        labs(x = "PCoA1 (25%)", y = "PCoA2 (13%)") +
        geom_segment(data = otu_score2, 
                     aes(x=x, xend=xend*0.5, y=y, yend= yend*0.5),
                     arrow = arrow(length = unit(0.2, "cm"))) + 
        geom_text(data= otu_score2, 
                  aes(x= xend*0.5, y = yend*0.5, label = taxa),
                  size=4,
                  color = "red")+
        theme_classic()

ggsave("documentation/pcoa_plot.tiff", width = 6, height = 7)
















