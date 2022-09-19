library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)


# read in meta.tidied_otu file 
otu_trimmed <- read_csv("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/otu_trimmed.csv")

head(otu_trimmed)
otu_trimmed <- otu_trimmed[,-1] # remove repeated id column  
dim(otu_trimmed)

# read in metadata 
metadata <- read_delim(file = "metadata_wu2011.txt")   # read in the metadata 
str(metadata)
colnames(metadata)

metadata2 <- metadata %>%  rename(sample.id = "#SampleID")


# inner join meta-data and otu table 

otu_meta <- otu_trimmed %>% inner_join(metadata2, by="sample.id")
head(otu_meta)
dim(otu_meta)

# extract only the abundance information 
com <- otu_meta[,2:179]  


# turn data frame into matrix 
m_com <- as.matrix(com)


# make community matrix  
set.seed(1015)
nmds <- metaMDS(m_com, distance = "bray") 
nmds 

# extract NMDS scores (x and y coordinates) 
data.scores <- as.data.frame(scores(nmds)) 


# add meta-data to the scores 
data.scores$diet <- otu_meta$DIET 

head(data.scores)


# plot ordination with ggplot2 
ggplot(data.scores, aes(x= NMDS1, y=NMDS2, color=diet))+
        geom_point(size = 3)+
        theme_classic()

 





