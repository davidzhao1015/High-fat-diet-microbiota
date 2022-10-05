library(tidyverse) 
library(vegan)



# Ordination of interest includes: PCA, PCoA, NMDS, CA, RDA, CCA and CAP 
# implement those ordination and then conduct hypothesis testing over groups 

# References: 1) Statistical analysis of microbiome data with R by Xia et al 
# 2) GUSTA ME (explains illustration for principle and ordination interpretation) 
# 3) minimalR (by Pat Schloss - doing almost everything with only tidyverse and vegan) 
# 4) QCBS R workshop series 



setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/")


## Principal Component Analysis (PCA)
# read in absolute abundance data 
abs_abund <- read_csv("processed_data/otutable_wide.csv") %>% 
        select(-1) %>% 
        column_to_rownames("sample_id") %>% 
        as.matrix()  # convert to matrix 
        

# read in meta data 
metadata <- read_csv("processed_data/metadata2.csv")  

metadata_diet <- metadata %>% 
        select(sample_id = sample.id, diet =DIET)


# absolute count data likely violates assumption of PCA so that data transformation is required 
# transformation method includes "total" or "hellinger" 

stand_abs_abund_hellinger <- decostand(abs_abund, 
                                       method = "hellinger") # apply hellinger method to transform count data 


stand_abs_abund_total <- decostand(abs_abund, 
                                   method = "total") # apply total method to transform count data  


# implement PCA with rda() of vegan package 
PCA_hellinger <- rda(stand_abs_abund_hellinger) # run PCA on hellinger-transformed count data 
PCA_hellinger # print summary 


PCA_total <- rda(stand_abs_abund_total) # run PCA on total-transformed count data 
PCA_total # print summary 



# quick plots for PCA of hellinger-transformed count data 
biplot(PCA_hellinger, display = "species")


ordiplot(PCA_hellinger, 
         display = "sites", 
         type="text")  # show sample names 

ordiplot(PCA_hellinger, 
         display = "sites", 
         type="points",
         scaling= 1)  # show points for samples  


# customize PCA plots 
# extract % explained by the first 2 axes 

# extract scores 

score_sample <- scores(PCA_hellinger,
       display = "sites",
       choices = c(1,2),
       scaling = 1) %>% 
        as_tibble(rownames = NA)  %>% # keep row names when converting to tittble 
        mutate(sample_id = rownames(.)) %>% 
        inner_join(metadata_diet, by="sample_id") %>%  # add diet 
        mutate(diet = factor(diet))



score_otu <- scores(PCA_hellinger,
                    display = "species",
                    choices = c(1,2),
                    scaling = 1) %>% 
        as_tibble(rownames = NA) 

# top 3 taxa contributing to PC1 and PC2 
pc1_taxa <- score_otu %>% 
        arrange(desc(abs(PC1))) %>% # arrange otu in descending order 
        mutate(taxa = rownames(.)) %>% 
        head(3)

pc2_taxa <- score_otu %>% 
        arrange(desc(abs(PC2))) %>% # arrange otu in descending order 
        mutate(taxa = rownames(.)) %>% 
        head(3)

pc_taxa <- bind_rows(pc1_taxa, pc2_taxa) %>% 
        mutate(x = 0, y=0) %>% 
        mutate(xend = PC1*0.3, yend=PC2*0.3) %>%  # shorten original arrow of each top variable 
        select(x, xend, y, yend, taxa) %>% 
        separate(taxa, 
                 into=c("kingdom", "phylum", "class", "order", "family","genus","species","strain"), 
                 sep=";")


ggplot(data = score_sample, aes(x=PC1, y= PC2)) +
        geom_point(aes(color = diet),
                   alpha=0.9,
                   size = 4,
                   fill = "black") +
        geom_segment(data = pc_taxa, 
                     aes(x=x, xend=xend, y=y, yend= yend),
                     arrow = arrow(length = unit(0.2, "cm"))) +
        geom_text(data= pc_taxa, 
                  aes(x= xend, y = yend, label = species),
                  size=4,
                  color = "red")+
        theme_classic()

ggsave("documentation/pca_plot.tiff",
       width = 8, height =5)        



































