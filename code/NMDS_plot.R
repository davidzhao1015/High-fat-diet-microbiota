library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)
library(ggtext)


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
data.scores$diet <- otu_meta$DIET  # add diet 
data.scores$sample.id <- otu_meta$sample.id # add sample id 

head(data.scores)

write.csv(data.scores, "./processed_data/data.scores.csv")  # write out nmds score data frame 

# plot ordination with ggplot2 
nmds_plot <- ggplot(data.scores, aes(x= NMDS1, y=NMDS2, color=diet))+
        geom_point(size = 3)+
        theme_classic()

ggsave("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/documentation/nmds_plot.tiff",
       width = 5, height = 4)
 

# add ellipse to nmds plot - refer to Pat Schloss instructions 
nmds_ellipse_plot <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, color = diet, fill= diet))+
        geom_point(size =3, show.legend = TRUE)+
        stat_ellipse(geom = "polygon", type="norm", level=0.8, alpha=0.2, show.legend = FALSE)+
        #geom_richtext(...) 
        coord_cartesian(xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))+
        labs(title = "Bray-Curtis dissimilarity \n between high- and low-fat groups",
             x= "NMDS Axis 1",
             y= "NMDS Axis 2")+
        scale_color_manual(name = NULL, 
                           breaks =c("HighFat", "LowFat"),
                           values= c("blue", "red"),
                           labels = c("High fat",  "Low fat"))+
        scale_fill_manual(name = NULL, 
                          breaks =c("HighFat", "LowFat"),
                          values= c("dodgerblue", "pink"),
                          labels = c("High fat",  "Low fat"))+ 
        theme_classic()

ggsave("./documentation/nmds_ellipse_plot.tiff", width = 5, height=5)
        
        


















