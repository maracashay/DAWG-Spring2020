#Network analyses are a way of evaluating correlations between all objects of interest in a study
#This is a good review on how network analyses can be applied on omics studies:
#https://doi.org/10.1093/bib/bbaa005

#network analyses can be compelted with several tools and methods of analyzing correlations (approaches are included in vegan, for example).
#P-values are calculated from the Total Information Coefficient (TIC) to answer the question: are they correlated/co-occuring?
#Maximal Information Coefficient (MIC) is becoming the preferred statistic to evaluate the strength of the association (see: Reshef et al. 2011 in Science )
#My understanding is that the "mictools" pipeline is a good option for balancing computation, power, and statistical approach

# Mictools pipeline is available on python, C++, and Docker
# These pipelines can break it down step-by-step, but in R, there are just a few lines of code to complete the entire pipeline
# For R: https://github.com/rsamantha/minerva
# Most recent publication/version of the MICtools pipeline: https://academic.oup.com/gigascience/article/7/4/giy032/4958979


### 1. Initial Setup ###

### Set the working directory ###

setwd("/Users/Fleishman/Desktop/DAWGFall2019/DAWG-Spring2020-master/April")

### Clear workspace ###

rm(list=ls())

### Install Packages ###
#install.packages("minerva") 
#devtools::install_github("briatte/ggnet")
#install.packages("fastDummies")


###Load packages ###
#ipak function for quickly loading/downloading packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# apply ipak to package list
packages <- c("minerva","phyloseq","microbiomeSeq","ggnet","dplyr","sna","fastDummies","vegan")
ipak(packages)

### Load PS object ###
load("/Users/Fleishman/Desktop/DAWGFall2019/DAWG-Spring2020-master/Feb28/DAWG 022820.rdata")

#### Prep input file from phyloseq object ####
# Minerva uses ONE matrix
# You can think of it as a taxa/count table with metadata file tacked on the end of it
# Rows are samples and each column is *binary* metadata or counts of a taxa

#subset to stool and tracheal for separate processing, we will start with "trachael"
ps.gen.stool<-subset_samples(ps.gen, Sample =="stool")
ps.gen.trach<-subset_samples(ps.gen, Sample =="tracheal")

genTable<-as.data.frame(ps.gen.stool@otu_table)
meta<-as.data.frame(ps.gen.stool@sam_data)

#check metadata to ensure it is all numeric and that categorical variables are binary
head(meta)
#it is all numeric, but not binary.
#we also note that some columns are not what we would want to include in the analysis

#remove extraneous info, such as "sampleID", sample, and Paired
meta$sampleID<-NULL
meta$Basespace_ID<-NULL
meta$Sample<-NULL
meta$Paired<-NULL

#convert categorical variables to dummy
#ID *could* (and maybe should) be dummy but it would make it complicated to interpret, so we will leave for our purposes
#similarly, "year", "month", etc probably should be dummy, but we can leave for now

#location
head(dummy_cols(meta$Location)) #produces dataframe with the original column followed by dummy columns
meta$home<-dummy_cols(meta$Location)[,2] #assign each location to its own new dummy column
meta$nursing<-dummy_cols(meta$Location)[,3]
meta$rehab<-dummy_cols(meta$Location)[,4]
meta$outsidehosp<-dummy_cols(meta$Location)[,5]
meta$otherlocal<-dummy_cols(meta$Location)[,6]
meta$Location<-NULL # remove the original location column

#gender
head(dummy_cols(meta$Gender)) #see that the second column is for zeros, or females
meta$female<-dummy_cols(meta$Gender)[,2]
meta$male<-dummy_cols(meta$Gender)[,3]
meta$Gender<-NULL

#Race
head(dummy_cols(meta$Race))
meta$white<-dummy_cols(meta$Race)[,2]
meta$black<-dummy_cols(meta$Race)[,3]
meta$Hispanic<-dummy_cols(meta$Race)[,3]
meta$Race<-NULL


#merge genera table and metadata
#Typically the input table is "short and fat" rather than "long and skinny" data
genera_meta<-cbind(genTable,meta)

#verify that the column names are the genera then metadata and that sample names carried through
colnames(genera_meta)
rownames(genera_meta)

netmatrix <- as.matrix(genera_meta)

#### MICtools pipeline ####

#1. Run mictools command to obtain inputs for a MIC evaluation
ticnull <- mictools(netmatrix, nperm=10000, seed=1234)


##Examine output:
names(ticnull) 
#our output includes the TIC, null distribution, observation distribution, and a table of p-values

head(ticnull$nulldist) #Null distribution of the data, this will be used to determine the strength (MIC)
head(ticnull$obsdist)#observed distribution of the data, this is what p-val are calculated from

hist(ticnull$tic) # observe distribution of TIC
hist(ticnull$pval$adj.P.Val, breaks=50, freq=FALSE) # observe distribution of adjusted p-values


#2. Use columns of indicies and FDR adjusted pvalue to calculate the strength of association (MIC)
#6 is the adjusted columns for pval.col. The 2,3 means the comparative values are in column 2 and 3 in the output
#default p-val cutoff is P<0.05
micres <- mic_strength(netmatrix, ticnull$pval, pval.col=c(6, 4, 5))  #pval.col =c(pval to use,variable column1, variable column 2) from the Ticnull matrix

#examine output:
head(micres)
#TicePval is the adjusted p-val 
#MIC is the strength of the relationship, greater numbers mean stronger relationship
#I1 abd I2 are the two variables being compared for each pval

hist(micres$MIC)
# It is common to have few observations with a strong (high MIC) relationship


#### Combine output for downstream visualization ####
#networks are usually only visualized with *significant* correlations, at a particular P-value
#Nodes (e.g. a specific taxa) are connected by lines (called edges)
#Nodes and/or edges can be colored, shaped, or sized based on other data: e.g. MIC, spearmans, metadata

#Network analyses are most commonly in cytoscape, which has many options: https://cytoscape.org/ 
#Here we will create a csv for upload to cytoscape, then use r visualization software for quick visualization

#1. replicate matrix
micres2<-micres


#2. Add in any other statistics that measure strength of relationship (that you may be interested in)

#To get CORRELATIONS you use the R cor function
netmatrix2<-as.data.frame(netmatrix)
#pearcorr = cor(netmatrix, method = c("pearson")) #linear, may not be appropriate for microbiome data
spearcorr = cor(netmatrix2, method = c("spearman")) #rank-based

#If you want more MINE statistics you can use the mine function
#mine(netmatrix)

##Filter for only significant relationships and add to MIC table, must do for each attribute calculated
#spearmans
micres2$spearcorr<-0
for (i in 1:nrow(micres2)) {
  micres2$spearcorr[i]<-spearcorr[micres2$I1[i], micres2$I2[i]]
}

#pearsons
#micres2$pearcorr<-0
#for (i in 1:nrow(micres2)) {
#  micres2$pearcorr[i]<-pearcorr[micres2$I1[i], micres2$I2[i]]
#}


#1. Create network file for quick visualization with ggnet2
micres2[is.na(micres2)]<-0 #NAs in the table cause errors, for our purposes we will just turn them to zeros
micres3<-select(micres2,I1,I2,MIC,spearcorr,TicePval)#reorder columns
net<-network(micres3, matrix.type = "edgelist", ignore.eval = F) #create network object
network.vertex.names(net) = colnames(netmatrix)#add node names to network object
quartz()
ggnet2(net,label=TRUE, edge.size = abs(micres3$MIC)) #visualize network object with ggnet2

#Very busy! This is easier to manage in cytoscape.
# Here, we can filter to only use lower value MIC (e.g. strong relationships) to do a quick visual
micres_filt<-micres3[abs(micres3$MIC)>0.5,] #filter MIC value
net_filt<-network(micres_filt, matrix.type = "edgelist", ignore.eval = F)
network.vertex.names(net_filt) = colnames(netmatrix)
quartz()
ggnet2(net_filt,label=TRUE, edge.size = abs(micres_filt$MIC))

#2. Observe in more detail if any metadata of interest have very strong relationships
micres_filt[micres_filt$I1=="Death"|micres_filt$I2=="Death",]
micres_filt[micres_filt$I1=="ICU"|micres_filt$I2=="ICU",]
micres_filt[micres_filt$I1=="CAP"|micres_filt$I2=="CAP",]
micres_filt[micres_filt$I1=="Sepsis"|micres_filt$I2=="Sepsis",]
#seems that while metadata are related to each other, there are not many strong relationships to the genera


#3. Output network csv for upload to cytoscape
#Your first import to Cytoscape will only submit the node and node columns, plus "edge attributes" (e.g. MIC, correlations, the p value, etc.) 
#You may import, in a second stage, the attribute information about your nodes (e.g. overarching categories for dummy variables or family assignments for genera). 
write.csv(micres2,"stoolnetworkp05.csv")

#4. create nodes attribute table
#this can be any manner that you want to categorize your nodes
#Here we will just categorize them as: genera, demographic variables, design variables, and health variables
length(genTable)#554 genera
colnames(meta)#1-3 experimental, 4 is demographic, 5-23 are health, 24 is experimental, 25-34 are demographic

nodedata<-as.data.frame(colnames(netmatrix))
colnames(nodedata)<-"nodes"
nodedata$type<-"genera"
nodedata$type[(length(genTable)+1):(length(genTable)+3)]<- "exp"
nodedata$type[(length(genTable)+4)]<- "dem"
nodedata$type[(length(genTable)+5):(length(genTable)+23)]<- "health"
nodedata$type[(length(genTable)+24)]<- "exp"
nodedata$type[(length(genTable)+25):(length(genTable)+34)]<- "exp"
nodedata[550:570,] #check

#filter for significant nodes
micres_nodedata<-subset(nodedata, nodedata$nodes %in% c(micres3$I1,micres3$I2))

#export
write.csv(micres_nodedata,"stoolnodedata.csv")



####Move into cytoscape for downstream analyses and improved visualization ####
#Network Analyzer: http://manual.cytoscape.org/en/latest/Network_Analyzer.html
#Notes for getting started:
#I have to open then save my network CSV file to have it upload correctly
#Their manual and tutorials are overly specific for importing: simply select I1 as the green circle and I2 as the red target
#when importing node data select the node names as the key



