library(ape) # load library for reading FASTA files 
library(tidyverse)
library(readr) # package to read in any tabular data  
library(vegan)
library(ggtext)
library(ggthemes)
library(RColorBrewer)
library(broom)
library(RVAideMemoire) # pair-wise comparison between groups in multivariate analysis 


# Permutation Multivariate Analysis of Variance (PERMANOVA) with adonis() [adonis2() may be the upgraded version?] 
# in vegan package 

# Ref
# 1) minimalR by Dr. Patrick Schloss
# 2) Paper: Application of multivariate statistical techniques in microbial ecology 
# 3) Chap 9 - Book: Statistical analysis of microbiome data with R  



# Objective of codes: test significance of difference in Bray-Curtis distance between two diet groups; 
# in addition, test significance of difference in Bray-Curtis distance by diet and sex. 





# Ensure analysis is reproducible
set.seed(19861015)

# Set permutation parameter 
permutation <- 1000 

# Change working directory 
setwd("C:/Users/17803/Documents/RStuodio-link-GitHub/Wu_2011_MLRepo/High-fat-diet-microbiota/") 




# Read in rarefied Bray-Curtis distance matrix 
bc_dist_rare <- read_csv("processed_data/beta_dist_rarefied_long.csv") %>%
        select(-1) # drop the first column (index)

# Read in the meta data 
meta <- read_csv("processed_data/metadata2.csv")  %>% 
        rename(sample_id = sample.id) %>% 
        select(sample_id, diet = DIET, sex = SEX) 



# join distance matrix and metadata 
dist_meta <- inner_join(bc_dist_rare, 
                        meta, 
                        by="sample_id") %>% 
        mutate(diet_sex = interaction(diet, sex))  # add a new interaction grouping diet * sex 
        
        

# convert to dist data 
all_dist <- dist_meta %>%
        select(all_of(.[["sample_id"]])) %>% 
        as.dist() 



# test significance for diet alone 
set.seed(19861015) 
adonis_diet <- adonis(all_dist~diet,
                      data = dist_meta,
                      permutations = permutation)

p.val_diet <- adonis_diet[["aov.tab"]][["Pr(>F)"]][1]  # significant different at 0.05 level  


# sequential test for [diet x sex] OR [sex X diet]  
set.seed(19861015)
adonis_diet_sex <- adonis(all_dist~diet*sex, data = dist_meta, permutations = 1000) 
adonis_diet_sex$aov.tab  # view output 
# diet explains 8.27% variance with p=0.001 while sex explains 8.54% variance 
# diet and sex have interaction effect on community structure with p=0.001 

set.seed(190861015)
adonis_sex_diet <- adonis(all_dist~sex*diet, data = dist_meta, permutations = 1000) 
adonis_sex_diet$aov.tab # result seems similar to the above although orders of two explanatory variables matters 



# test for newly-created four groups diet_sex 
set.seed(19861015)
adonis_group4 <- adonis(all_dist~diet_sex, data = dist_meta, permutations = 1000) 
adonis_group4$aov.tab
# based on Bray-Curtis distance, the group explains about 26.6% variance with significance at p=0.001 


# change the default dummy contrasts in R to the default "Sum" contrasts in vegan package 
set.seed(19861015) 
adonis_group4_contr.sum <- adonis(all_dist~diet_sex, 
                                  data = dist_meta, 
                                  permutations = 1000,
                                  contr.unordered = "contr.sum")

adonis_group4_contr.sum$aov.tab  # groups explain about 26.6% with p=0.001 


# add the contrasts for ordered factors 
set.seed(19861015)
adonis_group4_contr.poly <- adonis(all_dist~diet_sex, 
                                   data = dist_meta, 
                                   permutations = 1000,
                                   contr.unordered = "contr.sum",  # not quite sure this argument ... 
                                   contr.ordered = "contr.poly") # not quite sure this argument ... 

adonis_group4_contr.poly$aov.tab # groups explain about 26.6% with p=0.0010 




# post hoc analysis: pairwise permutation MANOVA using RVAideMemoire package [see p.297-301 Ref(3)]
set.seed(19861015) 

pairwise.perm.manova(all_dist, dist_meta$diet_sex,
                     test = c("Pillai","Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                     nperm = 1000,
                     progress = TRUE,
                     p.method = "fdr")  # adjust p-values using "FDR" method 




# evaluate group homogeneity using betadisper() [see p.301-304 Ref(3)]
group4 <- dist_meta$diet_sex

homo_to.centroid <- betadisper(all_dist, group4, type = "centroid") 
homo_to.centroid 

# percentage PCoA1 and PCoA2
homo_to.centroid$eig[1]/sum(homo_to.centroid$eig)  # 25.2% 
homo_to.centroid$eig[2]/sum(homo_to.centroid$eig) # 12.9% 

# alternative homogeneity testing 
homo_to.median <- betadisper(all_dist, group4, type = "median") 
homo_to.median 


# display betadisper object 
plot(homo_to.centroid)  # PCoA 

boxplot(homo_to.centroid)  # box plot illustrates association btw distance to centroid against four groups 


# customize PCoA illustrating dispersion homogeneity 
## extract scores of centers and samples in PCoA ordination 
homo_score_centroids <- scores(homo_to.centroid,
                         choices = c(1,2),
                         display = c("centroids")) %>% 
        as_tibble(rownames = "group4")  


homo_score_sample <- scores(homo_to.centroid,
                      choices = c(1,2),
                      display = c("sites")) %>% 
        as_tibble(rownames = "sample_id") %>%  
        inner_join(dist_meta, by="sample_id") %>%
        select(sample_id, PCoA1, PCoA2, diet_sex) %>% 
        full_join(homo_score_centroids, 
                  by = c("diet_sex" = "group4"), 
                  suffix = c("_sample", "_centroid")) 


# calculate the hulls for each group 
hull_group4 <- homo_score_sample %>% 
        group_by(diet_sex) %>%
        slice(chull(PCoA1_sample, PCoA2_sample)) %>% 
        ungroup() 

label <- homo_score_centroids 


ggplot(data = homo_score_sample) + 
        geom_point(aes(x= PCoA1_sample, 
                       y= PCoA2_sample, 
                       color = diet_sex))+
        geom_point(aes(x= PCoA1_centroid, 
                       y= PCoA2_centroid, 
                       color=diet_sex),
                   size = 6,
                   show.legend = F)+
        geom_segment(aes(x=PCoA1_centroid, xend =PCoA1_sample,
                         y=PCoA2_centroid, yend =PCoA2_sample), 
                     color = "grey") +
        geom_polygon(data = hull_group4, 
                     aes(x=PCoA1_sample, y=PCoA2_sample, color=diet_sex), 
                     fill = NA,
                     show.legend = FALSE)+
        # geom_text(aes(x= PCoA1_centroid,
        #               y= PCoA2_centroid,
        #               label= diet_sex,
        #               color = diet_sex),
        #           vjust = "top",
        #           hjust = "left") +
        geom_label(data = label, aes(x= PCoA1, y= PCoA2, label = group4, fill = group4),
                   nudge_x = 0.09,
                   fontface = "bold",
                   color="grey30",
                   show.legend = F,
                   alpha= 0.6) +
        coord_equal() +
        scale_color_manual(values = c("#B2182B", "#2166AC", "#D6604D", "#4393C3"),
                           breaks = c("HighFat.female", "HighFat.male", "LowFat.female", "LowFat.male"),
                           labels = c("High-fat, Female",
                                      "High-fat, Male",
                                      "Low-fat, Female",
                                      "Low-fat, Male"),
                           aesthetics = c("colour", "fill"))+
        # scale_fill_manual(values = c("#B2182B", "#2166AC", "#D6604D", "#4393C3"),
        #                   breaks = c("HighFat.female", "HighFat.male", "LowFat.female", "LowFat.male"))
        xlab("PCoA1 (25.2%)")+
        ylab("PCoA2 (12.9%)")+
        labs(title = "Homogeneity of dispersion in diet X sex groups",
             subtitle = "PCoA ordination on rarefied Bray-Curtis distance matrix",
             caption = "Bacterial community structure in four distinct groups were pairwisely different (p<0.01, PERMANOVA).\nDispersoin of the group Low-fat, female was significantly different from other three groups (p<0.05, betadisper() function)\nwhereas dispersion of remaining three groups were statistically equal (p>0.05, betadisper() function)")+ 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = NA,
                                              colour = "grey10",
                                              size = 0.5, linetype = "solid"),
              legend.position = "bottom",
              legend.title = element_blank(),
              plot.caption = element_text(hjust = 0),
              plot.caption.position= "panel")

ggsave("documentation/Homogeneityplot.tiff", width = 10, height =6)





# use either parametric ANOVA or permutation tests (permutest) to analyze the sig of the fitted model 
anova(homo_to.centroid) # p<0.001 

permutest(homo_to.centroid)  # p=0.001 
# both ANOVA and non-parametric test indicate that multivariate dispersion are significantly different 
# at 0.05 significance level 




# parametric Tukey's HSD test to compare distance to centroid in pairwise groups 
TukeyHSD(homo_to.centroid, conf.level = 0.95) 

# results show that low-fat female group is significantly different from other three groups whereas, 
# the remaining groups do not appear significantly different 
























