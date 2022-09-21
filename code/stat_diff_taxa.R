library(tidyverse)
library(broom)
library(purrr)

setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data")


meta <- read_csv("./metadata2.csv")  # meta data 
agg_genus <- read_csv("./agg_genus_data.csv")  # aggregate genus after being trimmed 


genus_meta <- inner_join(agg_genus, meta, by="sample.id")  # merge genus and meta data 

genus_meta$DIET <- as.factor(genus_meta$DIET)

# hypothesis testing - wilcox.test 
genus_test <- genus_meta %>% 
        nest(sample_data = c(-genus)) %>%
        mutate(test=map(sample_data, ~tidy(wilcox.test(prop~DIET, data=., exact=FALSE))))%>% 
        unnest(test) %>% 
        mutate(p.value.adj = p.adjust(p.value, method="BH")) %>%  # multiple test correction  
        arrange(p.value.adj) 

# pull genus those are significantly different by diets 
sig_genus <- genus_test %>% 
        filter(p.value.adj <= 0.05) %>% 
        pull(genus)

# box plot of significant genus 
boxplot_genus <- genus_meta %>% 
        filter(genus %in% sig_genus) %>% 
        mutate(genus = factor(genus, levels = sig_genus)) %>%
        mutate(prop= prop + 1/21000) %>% 
        ggplot(aes(x= genus, y=prop, color=DIET)) +
        geom_hline(yintercept = 1/10530, color = "gray") +
        geom_boxplot() +
        scale_color_manual(name=NULL,
                           values = c("blue", "red"),
                           breaks = c("HighFat", "LowFat"),
                           labels = c("High Fat", "Low Fat"))+
        labs(title = "Genus significantly associated with diets",
             x= NULL, 
             y="Relative Abundance (%)")+
        scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1),
                      labels=c(1e-2, 1e-1, 1, 10, 100))+
        theme_classic() +
        coord_flip()

ggsave("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/documentation/boxplot_genus.tiff",
       width = 5, height=6)








