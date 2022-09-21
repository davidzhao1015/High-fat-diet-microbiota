library(ape) # load library for reading fasta files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)


set.seed(1015)
setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/raw_data")


# READ IN RAW DATA 
metadata <- read_delim(file = "./metadata_wu2011.txt")   # read in the metadata 
str(metadata)
colnames(metadata)

metadata2 <- metadata %>%  
        rename(sample.id = "#SampleID")   # rename the sample id 


otu <- read_delim(file= "./refseq-based_otutable_wu2011.txt") # read in refseq-based otu table 
dim(otu)
otu_id <- otu[[1]]   # extract taxa 

# PREPROCESS INPUT DATAFRAME 
# transpose the otu table 
otu_t <- data.frame(t(otu[,-1])) # flip row to column with base R function 

colnames(otu_t) <- otu_id 

otu_t$sample.id <- rownames(otu_t) 

otu_t <- otu_t %>% select(sample.id, 1:361)  

otu_t_long <- otu_t %>% gather("taxa", "count", 2:362) 

otu_t_long2 <- otu_t_long %>% 
        separate(taxa, 
                into=c("kingdom", "phylum", "class", "order", "family","genus","species","strain"), 
                sep=";")

# re-calculate relative abundance at phylum level 
# check sequence counts in each sample - seems not equal 
otu_t_long2 %>% 
        group_by(sample.id) %>%
        summarise(n=sum(count)) %>%
        summary()

# add relative abundance for every otu 
otu_t_long3 <- otu_t_long2 %>% 
        group_by(sample.id) %>% 
        mutate(total_count = sum(count)) %>% 
        ungroup()

otu_t_long4 <- otu_t_long3 %>% 
        mutate(rel_abundance = count/total_count) 

otu_t_long4 %>% group_by(sample.id) %>% summarise(n=sum(rel_abundance)) %>% filter(n != 1)

# aggregate phylum count 
agg_phylum_data <- otu_t_long4 %>%
        group_by(sample.id, phylum) %>%
        summarise(prop = sum(rel_abundance)) 

agg_phylum_data %>% 
        group_by(phylum) %>% 
        summarise(median = median(prop)) %>% 
        arrange(desc(median))

# aggregate genus count 
agg_genus_data <- otu_t_long4 %>%
        group_by(sample.id, genus) %>%
        summarise(prop = sum(rel_abundance)) 

agg_genus_data %>% 
        group_by(genus) %>% 
        summarise(median = median(prop)) %>% 
        arrange(desc(median))


write.csv(agg_genus_data, 
          "C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/agg_genus_data.csv")

