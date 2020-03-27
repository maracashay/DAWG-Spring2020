### March 2020 DAWG Script ###
# Breakaway package maintained by Amy Willis

setwd("C:/Users/Nina/OneDrive - The Pennsylvania State University/DAWG")

#### Setup ####

### Clear workspace ###
rm(list=ls())

### Install Packages ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("...")

### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("devtools","magrittr","phyloseq","tidyverse","tibble","corncob","ggplot2")
ipak(packages)

library(phyloseq)
library(magrittr)
library(tidyverse)
library(tibble)
library(devtools)
library(ggplot2)

devtools::install_github("adw96/breakaway")
library(breakaway)

devtools::install_github("adw96/DivNet")
library(DivNet)


### Prep Phyloseq Objects ###

#Separate by sample type
ps.stool<-subset_samples(ps.adj, Sample =="stool")
ps.stool <- filter_taxa(ps.stool, function(x) sum(x) > 0, TRUE)
ps.trach<-subset_samples(ps.adj, Sample =="tracheal")
ps.trach <- filter_taxa(ps.trach, function(x) sum(x) > 0, TRUE)

#Extract the OTU tables, make samples columns
otu.stool <- as.data.frame(t(otu_table(ps.stool))) # 7342 ASVs
otu.trach <- as.data.frame(t(otu_table(ps.trach))) # 6834 ASVs

#### Run breakaway ####

#Collapse OTU tables into frequency tables
#A frequency table is a list with one entry for each sample.
#Each entry indicates how many species [freq] had a certain abundance [Var1]
frequencytable.stool <- build_frequency_count_tables(otu.stool)
frequencytable.trach <- build_frequency_count_tables(otu.trach)
head(frequencytable.stool[[65]])
plot(frequencytable.stool[[65]][,1],frequencytable.stool[[65]][,2])

head(frequencytable.stool[[2]])
plot(frequencytable.stool[[2]][,1],frequencytable.stool[[2]][,2])

#not enough data in most samples for Breakaway to be used
frequencytable.stool[[1]]
View(frequencytable.stool)

# run breakaway on individual sample
sample_richness(frequencytable.stool[[65]]) # 473
ba65 <- breakaway(frequencytable.stool[[65]]) # 473, same with useAll = T
plot(ba65)
sample_richness(frequencytable.stool[[2]]) # 19
breakaway_nof1(frequencytable.stool[[2]]) # 19 (exploratory function)

# run breakaway on all samples in ps object
richness.all <- breakaway(ps.adj)
richness.all$model
plot(richness.all, physeq=ps.adj, color="Day", shape = "Sample")
richness.stool <- breakaway(ps.stool) # 11 warnings, same with useAll = T
plot(richness.stool, physeq=ps.stool, color="month") # 8 rows with missing values
richness.trach <- breakaway(ps.trach) # 4 warnings, same with useAll = T
plot(richness.trach, physeq=ps.trach, color="Death") # 4 rows with missing values
summary(richness.stool)

#### Run DivNet ####

# doesn't work at ASV level because no ASV is observed in all samples
dv.stool <- DivNet::divnet(ps.stool, X = NULL)
dv.trach <- DivNet::divnet(ps.trach, X = NULL)

## make ps obejcts at genus level ##
#ps.gen <- tax_glom(ps.adj, taxrank="Phylum")
#ps.gen.stool<-subset_samples(ps.gen, Sample =="stool")
#ps.gen.trach<-subset_samples(ps.gen, Sample =="tracheal")

# takes a couple minutes to run
dv.stool <- divnet(ps.gen.stool, X = NULL)
#dv.trach <- divnet(ps.gen.trach, X = NULL)
#dv.all <- divnet(ps.gen, X = NULL, ncores = 4)

names(dv.stool) # gives estimates for these alpha diversity metrics

head(dv.stool$shannon)

dv.all %>% 
  plot(physeq=ps.gen, color = "Sample") +
  xlab("Sample Names") +
  ylab("Shannon diversity estimate\n(phylum level)") +
  coord_cartesian(ylim = c(0,2))

plot(dv.all)
