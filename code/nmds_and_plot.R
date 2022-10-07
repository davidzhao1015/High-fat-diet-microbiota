library(tidyverse)
library(vegan) 



# NMDS, non-metric multidimensional scaling 
# NMDS is a rank-based approach which means that the original distance data is substituted with ranks 
# While information about the magnitude of distances is lost, rank-based methods are generally more robust to 
# data which do not have an identifiable distribution 

# NMDS can 1) tolerate missing pairwise distance, 
# 2) be applicable to a dissimilarity matrix with any dissimilarity measure
# 3) use quantitative, semi-quantitative, qualitative, or mixed variables 


setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/") 

set.seed(19861015)


# read in rarefied Bray-Curtis distance matrix 

bc_rare <- read_csv("processed_data/beta_dist_rarefied_long.csv")  %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>% 
        as.matrix() # convert to a matrix 




# calculate NMDS 
bc_rare_nmds <- metaMDS(bc_rare) 

bc_rare_nmds # stress = 0.183   
# examine stress: >= 0.2 suspect; = 0.3 arbitrary; <= 0.1 fair; =<0.05 good 




# NMDS ordination plot with vegan package 
ordiplot(bc_rare_nmds, type="point") # plot samples 

# use ggplot2 to plot NMDS ordination 
metadata <- read_csv("processed_data/metadata2.csv") %>% 
        select(sample_id = sample.id, diet = DIET) %>% 
        mutate(diet = factor(diet))  # metadata 

sample_score <- scores(bc_rare_nmds, display = "sites") %>% 
        as_tibble(rownames = "sample_id") %>% 
        inner_join(., metadata, by = "sample_id")  # add metadata to sample scores 

ggplot(sample_score, aes(x= NMDS1, y=NMDS2, color=diet)) +
        geom_point()+
        stat_ellipse(level = 0.8)+
        annotate("text", x=-0.4, y= -0.3, label = "Stress=0.183", color = "blue")+
        theme_bw()
        



# stress plot 
par(mfrow = c(1,2))
stressplot(bc_rare_nmds) 
plot(bc_rare_nmds, display = "sites", type="point", main = "Goodness of fit")
points(bc_rare_nmds, display = "sites", cex=goodness(bc_rare_nmds)*300) # points bigger = worse fit 



























