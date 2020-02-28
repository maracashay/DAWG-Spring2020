#### Setup ####

### Clear workspace ###
rm(list=ls())

### Install Packages ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
BiocManager::install("DirichletMultinomial")

### load packages ###
### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("reshape2","magrittr","vegan", "dendextend")
ipak(packages)


library(microbiome)
library(DirichletMultinomial)
library(dplyr)
library(phyloseq)
library(ggplot2)

### Prep Phyloseq Objects ###
#genus level
#ps.gen<-taxa_level(ps.adj,"Genus") #microbiomeSeq package

#Separate by sample type
ps.gen.stool<-subset_samples(ps.gen, Sample =="stool")
ps.gen.trach<-subset_samples(ps.gen, Sample =="tracheal")

#Extract the OTU count matrix and convert it into right format
otu.stool.count <- as.matrix(t(abundances(ps.gen.stool)))
otu.trach.count <- as.matrix(t(abundances(ps.gen.trach)))

## Fit the DMM model *tracheal ##
fit <- mclapply(1:8, dmn, count=otu.trach.count, verbose=TRUE)

lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

best <- fit[[2]]
mixturewt(best)
predict <- apply(mixture(best), 1, which.max)
heatmapdmn(otu.trach.count, fit[[1]], fit[[2]], 30)

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.9))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

predict <- as.factor(predict)
predict <- as.data.frame(predict)
sample_data_new <- cbind(sample_data(ps.gen.trach), predict)
sample_data_new <- sample_data(sample_data_new)
ps.gen.trach.new <- phyloseq(otu_table(ps.gen.trach), sample_data_new)
sample_data(ps.gen.trach.new)$month <- as.factor(sample_data(ps.gen.trach.new)$month)

PCoA_trach <- ordinate(ps.gen.trach.new, "PCoA", "bray")
PCoA_trach.plot <- plot_ordination(ps.gen.trach.new, PCoA_trach, color="predict", shape="month")
PCoA_trach.plot + geom_point(size=3) + 
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,panel.background = element_blank(), axis.line = element_line(colour = "black"))

sample_trach <- sample_data(ps.gen.trach.new)

#Divide samples into stool and tracheal#
#and calculate alpha diversity#
ps.adj.stool <- subset_samples(ps.adj, Sample=="stool")
ps.adj.trach <- subset_samples(ps.adj, Sample=="tracheal")

tab1 <- richness(ps.gen.trach)
tab1

alpha <- cbind(sample_data(ps.gen.trach.new), predict, tab1)
alpha$predict <- as.factor(alpha$predict)
boxplot(chao1 ~ predict, data=alpha)
t.test(chao1 ~ predict, data=alpha)

##Grouping the model
month <- sample_data(ps.gen.trach)$month
names(month) <- rownames(otu.trach.count)

bestgrp <- dmngroup(otu.trach.count, month, k=1:5, 
                    verbose=TRUE, mc.preschedule=FALSE)


xtabs(~month + predict(bestgrp, otu.trach.count, assign=TRUE))
roc <- roc(month[rownames(otu.trach.count)]=="1", predict(bestgrp, otu.trach.count)[,"1"])
roc

xyplot(TruePostive ~ FalsePositive, roc, tyep="l")

test <- as.data.frame(otu_table(ps.gen))
write.csv(test, "test.csv")
write.csv(month, "month.csv")




#==========================================#

## Fit the DMM model *stool ##
fit_stool <- mclapply(1:8, dmn, count=otu.stool.count, verbose=TRUE)

lplc_stool <- sapply(fit_stool, laplace) # AIC / BIC / Laplace
aic_stool  <- sapply(fit_stool, AIC) # AIC / BIC / Laplace
bic_stool  <- sapply(fit_stool, BIC) # AIC / BIC / Laplace
plot(lplc_stool, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic_stool, type="b", lty = 2)
lines(bic_stool, type="b", lty = 3)

best_stool <- fit_stool[[2]]
mixturewt(best_stool)
predict_stool <- apply(mixture(best_stool), 1, which.max)
heatmapdmn(otu.stool.count, fit_stool[[1]], fit_stool[[2]], 30)

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.9))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

predict_stool <- as.factor(predict_stool)
predict_stool <- as.data.frame(predict_stool)
sample_data_new <- cbind(sample_data(ps.gen.stool), predict_stool)
sample_data_new <- sample_data(sample_data_new)
ps.gen.stool.new <- phyloseq(otu_table(ps.gen.stool), sample_data_new)


PCoA_stool <- ordinate(ps.gen.stool.new, "PCoA", "bray")
PCoA_stool.plot <- plot_ordination(ps.gen.stool.new, PCoA_stool, color="predict_stool")
PCoA_stool.plot + geom_point(size=3) + 
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,panel.background = element_blank(), axis.line = element_line(colour = "black"))



#### Heirarchacal clustering #####

dist <- vegdist(ps.pruned@otu_table, "bray")
clusters <- hclust(dist, method = "average")

plot(clusters)

group <- cutree(clusters, k=5)

groups <- as.data.frame(group)

groups$SampleID <- rownames(groups)

rownames(groups) <- NULL

groups <- groups[order(groups$group), ]
head(groups)




clust.dend <- as.dendrogram(clusters)
Day <- factor(ps.pruned@sam_data$Day)
n_Day <- length(unique(Day))
cols_2 <- colorspace::rainbow_hcl(n_Day, c = 70, l  = 50)
col_Day <- cols_2[Day]

k234 <- cutree(clust.dend, k = 2:4)

# color labels by Day:
labels_colors(clust.dend) <- col_Day[order.dendrogram(clust.dend)]
# color branches based on cutting the tree into 4 clusters:
clust.dend <- color_branches(clust.dend, k = 4)

### plots
par(mar = c(12,4,1,1))
plot(clust.dend)

# comparing methods

hc_complete <- hclust(dist, method = "complete")
hc_ave <- hclust(dist, method = "average")

hc_complete.dend <- as.dendrogram(hc_complete)
hc_ave.dend <- as.dendrogram(hc_ave)

tanglegram(hc_complete.dend, hc_ave.dend)
# similar subtrees are connected with lines of same color
# different subtrees are marked with a dashed line
