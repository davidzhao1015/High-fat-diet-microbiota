library(vegan)
library(tidyverse)

set.seed(19861015)  # ensure data analysis reproducibility  

setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/") 


# Read in OTU table 
add_genus_data <- read_csv("processed_data/agg_genus_data.csv") 
otu_raw <- read_delim("raw_data/refseq-based_otutable_wu2011.txt")

otu_raw_t <- otu_raw %>% 
        t() %>%
        as.data.frame() 

colnames(otu_raw_t) <- otu_raw_t[1, ] 

otu_raw_t2 <- otu_raw_t[-1, ]  # 95 samples * 361 taxa 
dim(otu_raw_t2)

otu_raw_t3 <- otu_raw_t2 %>% 
        mutate(sample_id = rownames(otu_raw_t2)) %>% # add sample id column  
        gather(key=taxa, value=read_count, 1:361) %>% # gather data frame  
        mutate(read_count = as.numeric(read_count))



# Library size per sample 
# determine the minimal library size for rarefaction   

min_librarysize <- otu_raw_t3 %>%
        mutate(sample_id =  factor(sample_id)) %>% # convert to factor variable 
        mutate(read_count = as.numeric(read_count)) %>% # convert to numeric variable 
        group_by(sample_id) %>% 
        summarize(n= sum(read_count))%>%
        summarize(min = min(n),
                  max = max(n)) %>% # max and min library size, min = 1892 

        pull(min)  # obtain the lowest library size for rarefaction  
        



# Apply vegan functions to calculate alpha diversity with rarefaction  
otu_raw_t3_df <- otu_raw_t3 %>% 
        spread(key= taxa, value = read_count, fill=0) %>% 
        as.data.frame() # prepare data frame for alpha diversity calculation 

write.csv(otu_raw_t3_df, "processed_data/otutable_wide.csv")  # write out the data frame 

alpha_rarefy <- vegan::rarefy(otu_raw_t3_df[,-1], sample=1814) # calculate richness with rarefaction   

alpha_rarefy_df <- data.frame(sample_id = otu_raw_t3_df[[1]],
                              richness_rarefied =  alpha_rarefy) # convert to data frame   



# Calculate also different types of measurement for alpha diversity 
alpha_multiple <- otu_raw_t3 %>% 
        group_by(sample_id) %>%
        summarise(sobs= specnumber(read_count),
                  shannon = diversity(read_count, index = "shannon"),
                  simpson = diversity(read_count, index = "simpson"),
                  invsimpson = 1/simpson,
                  n = sum(read_count)) %>% 
        inner_join(alpha_rarefy_df, by ="sample_id") %>% 
        select(sample_id:invsimpson, richness_rarefied, read_count =n) %>%
        mutate(richness_rarefied = round(richness_rarefied, 0))

# write out the alpha-diversity data including rarefied richness
write.csv(alpha_multiple, "processed_data/alpha_multiple.csv")  




# Beta diversity with rarefaction 
otu_raw_t3_matrix <- as.matrix(otu_raw_t3_df[,-1])
rownames(otu_raw_t3_matrix) <- otu_raw_t3_df[[1]]

# use vegan::avgdist to calculate rarefied Bray-Curtis distance matrix with 
beta_dist_rarefied <- avgdist(otu_raw_t3_matrix, 
                              dmethod = "bray", 
                              sample=min_librarysize)  

beta_dist_rarefied_df <- beta_dist_rarefied %>% 
        as.matrix() %>%
        as_tibble(rownames = "sample_id") %>%
        pivot_longer(-sample_id) %>% 
        filter(name < sample_id) # keep only the lower half 

# write out the lower half of the distance matrix 
write.csv(beta_dist_rarefied_df,
          "C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/beta_dist-matrix_rarefied.csv")

beta_dist_rarefied2 <- beta_dist_rarefied %>%
        as.matrix() %>% 
        as_tibble(rownames = "sample_id") # convert to tibble (alternative data frame)   

# write out rarefied Bray-Curtis distance
write.csv(beta_dist_rarefied2, "processed_data/beta_dist_rarefied_wide.csv")  





