library(ape) # load library for reading fasta files 
library(tidyverse)
library(readr) # package to read in any tabular data  



metadata <- read_delim(file = "metadata_wu2011.txt")   # read in the metadata 
str(metadata)
colnames(metadata)


otu <- read_delim(file= "refseq-based_otutable_wu2011.txt") # read in refseq-based otu table 
str(otu)
