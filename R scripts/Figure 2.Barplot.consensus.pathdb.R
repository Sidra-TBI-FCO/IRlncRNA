# Consensus db Enrichment visualization 

rm(list=ls())
setwd("~/")

# Load the required packages 
source("~/R scripts/ipak.function.R")
required.packages = c("func2vis")
ibiopak(required.packages)

#load data
tab1.exp = read.csv("./DEG/Overlap.Walk.diff_mRNA.TCGA.BRCA.ICR.H.L.normalized.127.log2.csv")
tab1.exp$X.1 = NULL

pathways <- read.table(("~/Downloads/ORA_results.tab"),header=TRUE,sep="\t")

#Rename "Gene" and "LogFC" columns to "gene" and "fc" respectively     (clean_pathways() function accepts "gene" and "fc" colnames only) 

tab1.exp$gene = tab1.exp$X
names(tab1.exp)[names(tab1.exp) == "logFC"] <- "fc"

cleaned.pathways = clean_pathways(tab1.exp, pathways)

dev.new()

#Stacked_barplot
svg(paste0("./Figures/Consensus.pathdb/TCGA.127.DEG.walk.ICR.HL.ConsensusDB.pathways.svg"), height = 4, width = 8, pointsize = 4)

plot_pathways_stacked_barplot(cleaned.pathways)
dev.off()

#Bubble plot
png(("./R figures/CPBD_enriched_pathwyas_bubble.png"), res = 600, width = 13, height = 10, units = "in")
plot_pathways(cleaned.pathways, 570, 13)
dev.off()