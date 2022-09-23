library(tidyverse)
library(broom)
library(purrr)



setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/processed_data")


meta <- read_csv("./metadata2.csv")  # meta data 
agg_genus <- read_csv("./agg_genus_data.csv")  # aggregate genus after being trimmed 

genus_meta <- inner_join(agg_genus, meta, by="sample.id")  # merge genus and meta data 

genus_meta$DIET <- as.factor(genus_meta$DIET) 

genus_meta2 <- genus_meta %>% 
        select(sample.id, DIET, genus, prop) # select necessary variables 


genus_pool <- genus_meta2 %>% 
        group_by(DIET, genus) %>% 
        summarize(median = median(prop), .groups = "drop") %>% 
        group_by(genus) %>%
        summarize(pool = max(median) <0.01,
                  median = max(median),
                  .groups ="drop") %>% 
        arrange(desc(median))


genus_diet_rel_abund <- inner_join(genus_meta2, genus_pool, by="genus") %>%
        mutate(genus = if_else(pool, "Other", as.character(genus))) %>% 
        group_by(sample.id, DIET, genus) %>%
        summarise(prop = sum(prop),
                  median=max(median),
                  .groups = "drop") %>% 
        mutate(genus = factor(genus),
               genus = fct_reorder(genus, median, .desc = FALSE)) 


# hypothesis testing - wilcox.test 
genus_significance <- genus_diet_rel_abund %>% 
        nest(sample_data = c(-genus)) %>%
        mutate(test=map(sample_data, ~tidy(wilcox.test(prop~DIET, data=., exact=FALSE))))%>% 
        unnest(test) %>% 
        mutate(p.value.adj = p.adjust(p.value, method="BH")) %>%  # multiple test correction  
        arrange(p.value.adj) %>%
        select(genus, sample_data, p.value.adj) %>% 
        filter(p.value.adj < 0.05) %>%
        filter(is.na(genus) != 1)

# box plot of significant genus 
boxplot_genus_stat <- genus_significance %>% 
        unnest(sample_data) %>% 
        select(genus, prop, DIET, p.value.adj) %>% 
        mutate(prop= prop + 1/21000) %>% 
        ggplot(aes(x= prop, y=genus, color=DIET)) +
        stat_summary(fun.data=median_hilow, geom = "pointrange",
                     fun.args=list(conf.int=0.5),
                     position = position_dodge(width=0.6)) +  # stat_summary is worth learning more 
        coord_trans(x="log10") +
        scale_x_continuous(limits=c(NA, 100),
                           breaks=c(0.01, 0.1, 1, 10, 100),
                           labels=c(0.01, 0.1, 1, 10, 100)) +
        scale_color_manual(name=NULL,
                           values = c("blue", "red"),
                           breaks = c("HighFat", "LowFat"),
                           labels = c("High Fat", "Low Fat"))+
        labs(title = "Genus significantly associated with diets",
             x= "Relative Abundance (%)", 
             y= NULL)+ 
        theme_classic() 

ggsave("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/documentation/boxplot_genus.tiff",
       width = 5, height=6)


# add stars and bars to the above figure 




























