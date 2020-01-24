ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("phyloseq", "DESeq2", "ALDEx2", "phyloseq", "ggplot2", "dplyr")
ipak(packages)

# BiocManager::install("ALDEx2")


#### Run PICRUSt2  Session 1, Jan 2020#####
#### Import PICRUSt2 tables ####

# import sample data
map <- read.table(file = "sample-metadata-pruned.txt", row.names =1, header=TRUE)

# import EC table
EC <- read.table(file = 'EC_pred_metagenome_unstrat.txt', sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)

# import pathways table
path <- read.table(file = "MetaCyc-path_abun_unstrat.txt", sep = '\t', header=TRUE, row.names=1, check.names = FALSE)

# make phyloseq objects for EC and Pathways
ps.EC <- merge_phyloseq(otu_table(t(EC), taxa_are_rows = FALSE), sample_data(map))
ps.path <- merge_phyloseq(otu_table(t(path), taxa_are_rows = FALSE), sample_data(map))
ps.EC@sam_data$month <- as.factor(ps.EC@sam_data$month)
ps.path@sam_data$month <- as.factor(ps.path@sam_data$month)

ps.EC.D1 <- subset_samples(ps.EC, Day =="1")
EC.D1 <- data.frame(otu_table(ps.EC.D1))
map.d1 <- data.frame(sample_data(ps.EC.D1))


EC.int <- mutate_all(EC.D1, function(x) as.integer(as.numeric(x)))
EC.int <- t(EC.int)
aldex2_EC <- ALDEx2::aldex(EC.int, map.d1$Sample, test="t", effect = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_EC, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

ps.path.D1 <- subset_samples(ps.path, Day =="1")
path.D1 <- data.frame(otu_table(ps.path.D1))

path.int <- mutate_all(path.D1, function(x) as.integer(as.numeric(x)))
path.int <- t(path.int)
aldex2_path <- ALDEx2::aldex(path.int, map.d1$Sample, test="t", effect = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_path, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# EC Tracheal ~ month
ps.EC.D1.trach <- subset_samples(ps.EC.D1, Sample == "tracheal")
map.D1.trach <- data.frame(sample_data(ps.EC.D1.trach))
path.D1.trach <- data.frame(otu_table(ps.EC.D1.trach))

path.int <- mutate_all(path.D1.trach, function(x) as.integer(as.numeric(x)))
path.int <- t(path.int)
aldex2_path <- ALDEx2::aldex(path.int, map.D1.trach$month, test="t", effect = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_path, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# EC Stool ~ month
ps.EC.D1.stool <- subset_samples(ps.EC.D1, Sample == "stool")
map.D1.stool <- data.frame(sample_data(ps.EC.D1.stool))
path.D1.stool <- data.frame(otu_table(ps.EC.D1.stool))

path.int <- mutate_all(path.D1.stool, function(x) as.integer(as.numeric(x)))
path.int <- t(path.int)
aldex2_path <- ALDEx2::aldex(path.int, map.D1.stool$month, test="t", effect = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_path, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# path Tracheal ~ month
ps.path.D1.trach <- subset_samples(ps.path.D1, Sample == "tracheal")
map.D1.trach <- data.frame(sample_data(ps.path.D1.trach))
path.D1.trach <- data.frame(otu_table(ps.path.D1.trach))

path.int <- mutate_all(path.D1.trach, function(x) as.integer(as.numeric(x)))
path.int <- t(path.int)
aldex2_path <- ALDEx2::aldex(path.int, map.D1.trach$month, test="t", effpatht = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_path, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# path Stool ~ month
ps.path.D1.stool <- subset_samples(ps.path.D1, Sample == "stool")
map.D1.stool <- data.frame(sample_data(ps.path.D1.stool))
path.D1.stool <- data.frame(otu_table(ps.path.D1.stool))

path.int <- mutate_all(path.D1.stool, function(x) as.integer(as.numeric(x)))
path.int <- t(path.int)
aldex2_path <- ALDEx2::aldex(path.int, map.D1.stool$month, test="t", effpatht = TRUE, denom="zero")
ALDEx2::aldex.plot(aldex2_path, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

# Nayaran et al 2020 in BMC Genomics used DESEq2 for analying picrust2
des.samp <- phyloseq_to_deseq2(ps.EC.D1, ~ Sample)
samp.counts <- counts(des.samp)
geoMeans = apply(samp.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.sample <-estimateSizeFactors(des.samp, geoMeans=geoMeans) 
dds.sample <- DESeq(dds.sample, test="Wald", fitType="parametric")
resultsNames(dds.sample)
alpha=0.05

res <- results(dds.sample, cooksCutoff = FALSE, contrast = c("Sample", "tracheal", "stool"))
sigtab.res = res[which(res$padj < alpha), ]
sigtab.res

# DESEq ~ Month EC #
des.month <- phyloseq_to_deseq2(ps.EC.D1, ~ month)
month.counts <- counts(des.month)
geoMeans = apply(month.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.month <-estimateSizeFactors(des.month, geoMeans=geoMeans) 
dds.month <- DESeq(dds.month, test="Wald", fitType="parametric")
resultsNames(dds.month)

res <- results(dds.month, cooksCutoff = FALSE, contrast = c("month", "0", "1"))
sigtab.res = res[which(res$padj < alpha), ]
sigtab.res

# EC tracheal ~ month
des.month <- phyloseq_to_deseq2(ps.EC.D1.trach, ~ month)
month.counts <- counts(des.month)
geoMeans = apply(month.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.month <-estimateSizeFactors(des.month, geoMeans=geoMeans) 
dds.month <- DESeq(dds.month, test="Wald", fitType="parametric")
resultsNames(dds.month)

res <- results(dds.month, cooksCutoff = FALSE, contrast = c("month", "0", "1"))
sigtab.res = res[which(res$padj < alpha), ]
sigtab.res

# EC stool ~ month
ps.EC.D1.stool <- subset_samples(ps.EC.D1, Sample == "stool")
des.month <- phyloseq_to_deseq2(ps.EC.D1.stool, ~ month)
month.counts <- counts(des.month)
geoMeans = apply(month.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.month <-estimateSizeFactors(des.month, geoMeans=geoMeans) 
dds.month <- DESeq(dds.month, test="Wald", fitType="parametric")
resultsNames(dds.month)

res <- results(dds.month, cooksCutoff = FALSE, contrast = c("month", "0", "1"))
sigtab.res = res[which(res$padj < alpha), ]
sigtab.res
sigtab.res.df = data.frame(sigtab.res) # convers the output to a data frame
sigtab.res.df = cbind(EC=row.names(sigtab.res.df), sigtab.res.df)
sigtab.res.df = sigtab.res.df[order(sigtab.res.df$log2FoldChange),]
# orders the ESVs by the log 2 fold change
sigtab.res.df


stool.month <- ggplot(sigtab.res.df, aes(x=log2FoldChange, y=reorder(EC, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="EC") + 
  theme(axis.text.x = element_text(hjust = 0, size=10),axis.text.y= element_text(size=10)) +
  theme(axis.title= element_text(size=10))
stool.month = stool.month + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
stool.month
# EC 3.6.3.4. is a Cu 2+ exporting ATPase
# EC 4.1.1.49 is a phosphoenolpyruvate carboxykinase enzyme
# EC 3.4.11.4 is a tripeptide aminopeptidase enzyme


EC.D1.Stool <- data.frame(otu_table(ps.EC.D1.stool))
map.D1.stool <- data.frame(sample_data(ps.EC.D1.stool))


plot(log(EC.D1.Stool$EC.4.3.1.1) ~ map.D1.stool$month)
# 0 = no mortality in 1 month
# EC 4.3.1.1 is an aspartate ammonia-lyase 

plot(EC.D1.Stool$EC.3.1.1.61 ~ map.D1.stool$month)
# EC 3.1.1.61 is a protein glutamate methylesterase also important for chemotaxis

