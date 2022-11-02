library(edgeR)  # model on over-dispersed count data 
library(DESeq2) # model on over-dispersed count data 
library(GUniFrac)
library(tidyverse)
library(statmod)



# Modeling over-dispersed microbiome data - Chap 11 in Book: Statistical analysis of microbiome data with R 


# Principles: 
# 1. high-throughput sequencing data sets were advised to be treated as count data; and the methods based on 
# count distribution or correspondence analysis were often used in the literature 
# 2. biological variations: biological replicates are used to test the biological variations between samples 
# within group, thus providing information that is necessary to make inferences between groups, and to generalize 
# conclusions. 
# 3. technical variations: technical replicates aim to test the variation in the testing protocol itself 
# 4. the count models: advantage for separating biological from technical variation
# 5. RNA-seq data 

# models for over-dispersed microbiome data: 
# poisson model (not suitable to deal with over-dispersion existing in microbiome count data)
# Negative Binomial (NB) model 
# Wald test 
# NB-based exact test 


# implement edgeR - chap 11.3.2 
# steps: 1) normalization, 2) dispersion estimation, 3) test for differential expression 

require(GUniFrac)  # it provides microbiome data set of the study by Charlson et al 2010 

data("throat.otu.tab")

head(throat.otu.tab)

throat <- t(throat.otu.tab)  # transpose to taxa-sample format 

counts <- throat 

head(counts)



require(edgeR) # load package for modeling over-dispersed microbiome data 

data("throat.meta")

group <- throat.meta$SmokingStatus  # grouping variables 

head(group)

dim(counts) # 856 taxa, 60 samples 

length(group)  # 60 samples 

y <- DGEList(counts = counts, group = group) # build the edgeR object 

names(y)  # elements that object contains 

head(y$counts)

y$samples # contains a summary of samples 

sum(y$all.zeros) # how many genes have 0 counts across all samples 





# step 3. filter taxa with too low counts in any of the exp conditions 
# based on counts per million (CPM), remove taxa with a CPM value less than the cutoff from the analysis 

dim(y)  # 856 taxa, 60 samples 

y_full <- y # keep the old one in case we mess up  

head(y$counts)

apply(y$counts, 2, sum)  # total OTU counts per sample 

keep <- rowSums(cpm(y) > 100) >=2 # keep OTU with a cpm of 100 or greater at least two samples 

y <- y[keep, ]

dim(y) # reduce the data set from 856 to 616 OTUs. 

y$samples$lib.size  <- colSums(y$counts) # reset the library sizes 



# normalize the data 
y <- calcNormFactors(y)
y
y$samples # without the replacement, the default value is 1 for all values in y$samples$norm.factors 

# effective lib size is the product of the original lib size and the scaling factor 
y$samples$lib.size*y$samples$norm.factors 




# step 5. explore the data by multi-dimensional scaling (MDS) plot 

pdf("MDS_plot.pdf", width = 7, height = 7) # in inches 

plotMDS(y,
        method = "bcv",
        main = "MDS plot for throat count data",
        col = as.numeric(y$samples$group),
        labels = colnames(y$counts)) 

legend("topright",
       as.character(unique(y$samples$group)),
       col=1:2,
       cex=0.8,
       pch=16)

dev.off() # tells R to turn off device and writing to the pdf 




# step 6. estimate the dispersion parameter for each taxon/ OTU 
# dispersion estimation is critical for NB modeling 

# the dispersion measures the biological variability of within-group variability i.e. variability 
# between replicates (or called inter-library variation) for that OTU 

# for strongly abundant OTU, the dispersion = (variation between samples of the same treatment)^2 



# estimate common dispersion - overall variability across data
y1 <- estimateCommonDisp(y, verbose = T) 

names(y1)

# fit a trended model to get a tag/ taxon wise dispersion 
y1 <- estimateTagwiseDisp(y1)  

names(y1)

plotBCV(y1)  # plot the taxon-wise biological coefficient of variation (sq root of dispersion) vs log2-CPM 



# fit a generalized linear model to estimate the taxon-wise dispersion 
# design matrix 
design <- model.matrix(~group) 
rownames(design) <- colnames(y)
design

require(statmod)

y2 <- estimateDisp(y, design, robust = TRUE)
y2$common.dispersion 

plotBCV(y2) # plot the taxon-wise biological coefficient of variation (sq root of dispersion) vs log2-CPM
#the plot shows that the trended dispersion deceases with taxon count. at low logCPM, the dispersion are 
#very large indeed. 

fit <- glmQLFit(y2, design, robust = TRUE)  # quasi-likelihood GLM 
plotQLDisp(fit)



# step 7. test the differential abundance 
# the exactTest() approach - make the pairwise comparisons between the groups 
# null hypothesis: the observed difference between Smoker and NonSmoker was merely caused by 
# experimental variability i.e., the type of variability that we can just as well expect between 
# different samples in the same group 

et <- exactTest(y1, pair = c("NonSmoker", "Smoker")) 
topTags(et)  # tabulate the top differential abundant OTU/ taxon 

y3 <- y
y3$samples$group <- relevel(y3$samples$group, ref="Smoker") # re-level "Smoker" as the control or ref level 
levels(y3$samples$group)



# GLM approach 
# GLM method is similar to but more feasible than exactTest() approach 

# create design matrix  to describe the treatment conditions 

design <- model.matrix(~group)
rownames(design) <- colnames(y)
design 

fit <- glmQLFit(y1, design)  #first conduct GLM F test and likelihood ratio test using this design as one argument 
qlf <- glmQLFTest(fit, contrast = c(-1,1))  # use glmQLFTest() with contrast argument to compare smokers vs non-smoker
topTags(qlf)

FDR <- p.adjust(qlf$table$PValue, method = "BH")

sum(FDR < 0.05)

topTags(qlf, n=15)



# step 8. interpret the results of differential expression analysis with diagnostic plots 

# MA-plot using plotSmear() 

da = decideTestsDGE(et, p.value = 0.1)
da_OTUs <- rownames(y1)[as.logical(da)]
plotSmear(et, de.tags = da_OTUs, cex = .5)
abline(h = c(-2,2), col = "blue")



# volcano plot - summarize both fold-change and a measure of statistical test (ie. p-val)

tab <- data.frame(logFC = et$table[,1],
                  negLogPval = - log10(et$table[,3]))

head(tab)

par(mar = c(5,4,4,4)) 
plot(tab, 
     pch = 16, 
     cex=0.6,
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue)) 

# log2 fold change and p-value cutoffs 
lfc = 2
pval = 0.1 

# selecting interest OTUs 
sig_OTUs <- (abs(tab$logFC) >lfc & tab$negLogPval > -log10(pval)) 

# identify the selected OTUs 
points(tab[sig_OTUs, ],
       pch=16,
       cex=0.8,
       col="red") 
abline(h=-log10(pval),
       col = "green3",
       lty=2)
abline(v=c(-lfc, lfc),
       col="blue",
       lty=2)
mtext(paste("pval=", pval),
      side=4,
      at=-log10(pval),
      cex=0.8,
      line=0.5,
      las=1)
mtext(c(paste("-", lfc, "fold"),
        paste("+", lfc,"fold")),
      side=3,
      at=c(-lfc, lfc),
      cex = 0.8,
      line = 0.5)






