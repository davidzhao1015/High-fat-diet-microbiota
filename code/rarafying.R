library(vegan)
library(tidyverse)


# CC188 using the vegan r package to generate ecological distances avgdist() taking consideration of rarefaction  

# CC189 (to watch) rarefying ecological distances with r: should you? 
# CC190 is normalization an acceptable alternative to rarefaction? Nope. 
# CC191 differences in sampling effort impact bray-Curtis distances and rarefaction can minimize it 
# CC198 generate a rarefaction curve from collector's curves in r within the tidyverse (to watch) 
# CC200 how to rarefy community data with tidyverse and vegan packages 
# CC201 is richness estimation an alternative to rarefaction? trying breakaway and chao1 
# CC202 how to find the best sampling depth for rarefaction 

set.seed(19861015)

add_genus_data <- read_csv("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/agg_genus_data.csv")


otu_raw <- read_delim("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/raw_data/refseq-based_otutable_wu2011.txt")

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


# library size per sample 
min_librarysize <- otu_raw_t3 %>%
        mutate(sample_id =  factor(sample_id)) %>% # convert to factor variable 
        mutate(read_count = as.numeric(read_count)) %>% # convert to numeric variable 
        group_by(sample_id) %>% 
        summarize(n= sum(read_count))%>%
        summarize(min = min(n),
                  max = max(n)) %>% # max and min library size, min = 1892 

        pull(min)  # obtain the lowest library size for rarefaction  
        
        
# apply vegan functions to calculate alpha diversity with rarefaction  
otu_raw_t3_df <- otu_raw_t3 %>% 
        spread(key= taxa, value = read_count, fill=0) %>% 
        as.data.frame() # prepare data frame for alpha diversity calculation 

write.csv(otu_raw_t3_df,
          "C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/otutable_wide.csv")



alpha_rarefy <- vegan::rarefy(otu_raw_t3_df[,-1], sample=1814) 

alpha_rarefy_df <- data.frame(sample_id = otu_raw_t3_df[[1]],
                              richness_rarefied =  alpha_rarefy) # richness with rarefaction based on vegan package 


# calculate also different types of measurement for alpha diversity 
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

write.csv(alpha_multiple, 
          "C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/alpha_multiple.csv")


# beta diversity with rarefaction 
otu_raw_t3_matrix <- as.matrix(otu_raw_t3_df[,-1])
rownames(otu_raw_t3_matrix) <- otu_raw_t3_df[[1]]


beta_dist_rarefied <- avgdist(otu_raw_t3_matrix, 
                              dmethod = "bray", 
                              sample=min_librarysize)  

beta_dist_rarefied_df <- beta_dist_rarefied %>% # bray-curtis with rarefaction 
        as.matrix() %>%
        as_tibble(rownames = "sample_id") %>%
        pivot_longer(-sample_id) %>% 
        filter(name < sample_id) # keep only the lower half 

write.csv(beta_dist_rarefied_df,
          "C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/beta_dist-matrix_rarefied.csv")


# NMDS - CC187 by pat schloss 
beta_dist_rarefied2 <- beta_dist_rarefied %>%
        as.matrix() %>% 
        as_tibble(rownames = "sample_id")

nmds_bc <- metaMDS(beta_dist_rarefied) 

stress_bc <- nmds_bc$stress 

metadata2 <- read_csv("~/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data/metadata2.csv")

scores(nmds_bc) %>%
        as_tibble(rownames = "sample.id") %>%
        inner_join(., metadata2, by="sample.id") %>%
        ggplot(aes(x=NMDS1, y=NMDS2, color=DIET)) +
        geom_point()




