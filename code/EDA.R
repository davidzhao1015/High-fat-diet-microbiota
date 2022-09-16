library(ape) # load library for reading fasta files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)


# READ IN RAW DATA 
metadata <- read_delim(file = "metadata_wu2011.txt")   # read in the metadata 
str(metadata)
colnames(metadata)


otu <- read_delim(file= "refseq-based_otutable_wu2011.txt") # read in refseq-based otu table 
str(otu)
head(otu)

otu_id <- otu[[1]]   # extract taxa 


# PREPROCESS INPUT DATAFRAME 
# transpose the otu table 
otu_t <- data.frame(t(otu[,-1])) # flip row to column with base R function 
head(otu_t)

colnames(otu_t) <- otu_id # assign new column names as the taxa 

colnames(otu_t) 


# full taxa seems too long as column labels. 
# how to show only the name at the highest resolution?
require(stringr)

otu_taxa_matrix <- str_split(otu_id, pattern = ";", simplify = TRUE)
head(otu_taxa_matrix) 
dim(otu_taxa_matrix)


dim(otu_t)  # 95 observations, 360 taxa 

# missing values ?
sum(is.na(otu_t)) # no missingness 


# calculate relative abundance across the row based on count numbers 
sample.id <- rownames(otu_t) 

otu_t_id <- cbind(sample.id, otu_t)  

otu_t2 <- as.matrix(otu_t) 

otu_t_relative <- proportions(otu_t2, margin=1)   

rowSums(otu_t_relative)


otu_t_relative2 <- as.data.frame(otu_t_relative) 
rownames(otu_t_relative2)


# trim rare taxa with presence less than 10% or proportion lower than 0.1% in all samples 
presence <- apply(otu_t_relative2, 2, function(x) mean(x !=0))
sum(presence < 0.1)  # 181 rare taxa whereas 180 common taxa 

rara_taxa <- names(presence)[c(which(presence <0.1))]   # store names of rare taxa 

common_taxa <- names(presence)[c(which(presence >=0.1))] # store 180 common taxa 

otu_t_relative2_trim  <- otu_t_relative2 %>% select(all_of(common_taxa)) 
dim(otu_t_relative2_trim)  # done! 180 common taxa with 95 samples 



low.prop <- apply(otu_t_relative2_trim, 2, function(x) sum(x <0.001)) # prop less than 0.1% in each sample 
non.low.prop_taxa <- names(low.prop)[c(which(low.prop !=95))] # drop additional 2 rare taxa 

otu_t_relative2_trim2 <- otu_t_relative2_trim %>% select(all_of(non.low.prop_taxa))  
dim(otu_t_relative2_trim2) # done! 178 common taxa with 95 samples 

head(otu_t_relative2_trim2) # df used for the following analysis 


# reshape data frame 

otu_t_relative2_trim3 <- cbind(sample.id, otu_t_relative2_trim2) # add sample id to the first column 
dim(otu_t_relative2_trim3)

otu_t_relative2_trim3_long <- otu_t_relative2_trim3 %>% gather(taxa, prop, 2:179) # make a long table 
head(otu_t_relative2_trim3_long)

# bar plot without showing y-axis label ... 














