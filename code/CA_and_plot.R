library(tidyverse)
library(vegan)



# Correspondence analysis (CA)
# CA uses chi-square distances and gives high weight to rare species (low occurence species with many zero)
# answer "which sites do my species prefer?" or "which sites to my species correspond to?" 

# key assumption: 
# 1) variables have a uni-modal distribution across objects 
# 2) variables are dimensionally homogeneous (i.e., have the same units)
# 3) all variables are either zero or positive values 



setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")

set.seed(19861015)



# read in otu count (abs abundance) data 
otu_abs_abund <- read_csv("processed_data/otutable_wide.csv")  %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>% 
        as.matrix() 


otu_cca <- cca(otu_abs_abund) # compute CA 

otu_cca # print output 
# total heterogeneity of data (i.e., inertia) is 2.6, and first axis captures 14.3% = 0.3751/2.605 where
# 0.3751 is eigenvalue of the first axis CA1, and 2.605 is the total heterogeneity of the data 


plot(otu_cca, display = "sites")  # display samples only 
ordiplot(otu_cca, type= "points", display = "all") # add taxa to plot 

# customize CA plot 
metadata <- read_csv("processed_data/metadata2.csv") %>% 
        select(sample_id = sample.id, diet = DIET) %>% 
        mutate(diet = factor(diet))  # metadata 



sample_score <- otu_cca$CA$u %>% 
        as_tibble(rownames = "sample_id") %>% 
        select(sample_id, CA1, CA2) %>%
        inner_join(., metadata, by="sample_id")


taxa_score <- otu_cca$CA$v %>% 
        as_tibble(rownames = "taxa") %>% 
        select(taxa, CA1, CA2) 

ggplot(sample_score, aes(x=CA1, y=CA2,color=diet))+
        geom_point(size = 3)+
        geom_point(data = taxa_score, aes(x=CA1, y=CA2), 
                   inherit.aes = F,
                   color = "red",
                   shape= 2, 
                   alpha = 0.6)+
        theme_bw()




# post- CA analysis 
# apply Keiser-Guttman criterion and Broken-stick model to evaluate performance of CA axis 

source("code/evplot.R")  # import the function, evplot() 
ev <- otu_cca$CA$eig # extract eigenvalues 

evplot(ev)  # plot 





























