# Calculate the estimate stromal, immune scores in correlation with the 3Lnc. signature 
#clean global environment
rm(list=ls())

#Set working directory
load("~/R.Config.Rdata")
setwd("~/")

# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
#Define parameters
required.packages <- c("corrplot", "stringr","utils","estimate")
ipak(required.packages)

# Set parameters
approach = "estimate" #choices are: "consensusTME" or "consensusTME_TBI" or "standard_ssGSEA"
Cancertypes = c("BRCA","RAQA","UCEC", "OV",  "READ", "COAD","HNSC","SKCM", "LIHC", "STAD", "BLCA",  "CESC", "KICH","LUSC", "LUAD", "GBM","KIRP","KIRC","LGG")

N.Cancertypes = length(Cancertypes)

i=1
for (i in 1:N.Cancertypes)
  {
Cancer = Cancertypes[i]

# Load Corresponding data
load(paste0("./Processed.data/TCGA-", Cancer, "/003_TCGA-", Cancer, "_EDAseq_normalized.gene.info.2017.Rdata"))

colnames(RNASeq.NORM.quantiles) = substr(colnames(RNASeq.NORM.quantiles),1,12)
colnames(RNASeq.NORM.quantiles) = gsub("\\.","-",colnames(RNASeq.NORM.quantiles))
RNASeq.NORM.quantiles = RNASeq.NORM.quantiles[,!duplicated(colnames(RNASeq.NORM.quantiles))]
dim(RNASeq.NORM.quantiles)

# Log2 transformation
RNASeq.QN.counts.filtered = log(RNASeq.NORM.quantiles+1,2)

dir.create("./Validation.dataset/Analysis/signatures/Estimate/",showWarnings = FALSE)  
write.table(RNASeq.QN.counts.filtered, sep = "\t", file = paste0("./Validation.dataset/Analysis/signatures/Estimate/003_",Cancer,"_normalized_TCGAbiolinks_TP_filtered.txt"), quote = FALSE)

filterCommonGenes(input.f=paste0("./Validation.dataset/Analysis/signatures/Estimate/003_",Cancer,"_normalized_TCGAbiolinks_TP_filtered.txt"),
                  output.f=paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,".filtered.ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))

estimateScore(input.ds = paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,".filtered.ESTIMATE.input.gct"),
              output.ds = paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,".filtered.ESTIMATE.score.gct"),
              platform= "illumina")

estimate_.gct<-read.table(paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,".filtered.ESTIMATE.score.gct"), skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate_.gct) = estimate_.gct$NAME
estimate_.gct$NAME = NULL
estimate_.gct$Description = NULL

ESTIMATE = t(estimate_.gct)
save(ESTIMATE, file = paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,"_ESTIMATE_scores.Rdata"))

}

######################################### Scatter plot of the Estimate score and the 3Lnc signature score ######################################

## Setup environment
rm(list=ls())
setwd("~/")

load("~/R.Config.Rdata")
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages = c("ggplot2", "ggpubr")
ipak(required.packages)

# Set Parameters
var2 = "ESTIMATE Score"
var1 = "3 Lnc Sig score"

 Cancertypes = c("BRCA","RAQA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")

# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=1
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]


# Load data
load(paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,"_ESTIMATE_scores.Rdata"))
load(paste0("./Validation.dataset/data.file.ICR.ssgsea.scaled.Rdata"))

rownames(ESTIMATE) = gsub("\\.","-",rownames(ESTIMATE))

all_scores_df = all_scores_df[which(all_scores_df$X %in% rownames(ESTIMATE)),]

ESTIMATE = as.data.frame(ESTIMATE)

# Analysis
ESTIMATE$sig.3.lnc.Score = all_scores_df$SSGSEA_Score[match(rownames(ESTIMATE), all_scores_df$X)]

ESTIMATE = ESTIMATE[which(rownames(ESTIMATE) %in% all_scores_df$X),] 
  
plot_df = data.frame(Sample = row.names(ESTIMATE), var1 = NA, var2 = NA)

plot_df$var1 = ESTIMATE$sig.3.lnc.Score
plot_df$var2 = ESTIMATE$ESTIMATEScore

score = "ESTIMATEScore"   ## "StromalScore"    "ImmuneScore"     "ESTIMATEScore"   "sig.3.lnc.Score"  

#dev.new()

#Plotting 
dir.create("./Validation.dataset/Figures/Correlations/Scattered_plots",showWarnings = FALSE)
dir.create(paste0("./Validation.dataset/Figures/Correlations/Scattered_plots/",score),showWarnings = FALSE)
svg(filename = paste0("./Validation.dataset/Figures/Correlations/Scattered_plots/",score,"/",var2,"_",var1,"_stattered_plot_",Cancer,".svg"), height =4 ,width = 5)


plot = ggplot(plot_df, aes(var1, var2)) +
  geom_point(size = 0.4) +
  stat_cor(method = "pearson", size = 5) +
  geom_smooth(method="lm") +
  ylab(var2) +
  xlab(var1) +
  theme_bw() +
  ggtitle(Cancer) +
  theme(axis.title.x = element_text(size = 17, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        axis.text.x = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 17, colour = "black"))
plot(plot)
dev.off()

}

#### Scatter plot of the Estimate score and the immune score

# Set Parameters
var2 = "ESTIMATE Score"
var1 = "Immune Score"

# Need to add RAQA 
Cancertypes = c("BRCA","RAQA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")

# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=1
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  
  # Load data
  load(paste0("./Validation.dataset/Analysis/signatures/Estimate/TCGA.",Cancer,"_ESTIMATE_scores.Rdata"))
  ESTIMATE = as.data.frame(ESTIMATE)
  
  plot_df = data.frame(Sample = row.names(ESTIMATE), var1 = NA, var2 = NA)
  
  plot_df$var1 = ESTIMATE$ImmuneScore
  plot_df$var2 = ESTIMATE$ESTIMATEScore
  
  score = "ESTIMATEScore"   ## "StromalScore"    "ImmuneScore"     "ESTIMATEScore"   "sig.3.lnc.Score"  
  
  #Plotting 
  dir.create("./Validation.dataset/Figures/Correlations/Scattered_plots",showWarnings = FALSE)
  dir.create(paste0("./Validation.dataset/Figures/Correlations/Scattered_plots/",score),showWarnings = FALSE)
  svg(filename = paste0("./Validation.dataset/Figures/Correlations/Scattered_plots/",score,"/",var2,"_",var1,"_stattered_plot_",Cancer,".svg"), height =4 ,width = 5)
  
  
  plot = ggplot(plot_df, aes(var1, var2)) +
    geom_point(size = 0.4) +
    stat_cor(method = "pearson", size = 5) +
    geom_smooth(method="lm") +
    ylab(var2) +
    xlab(var1) +
    theme_bw() +
    ggtitle(Cancer) +
    theme(axis.title.x = element_text(size = 17, colour = "black"),
          axis.title.y = element_text(size = 17, colour = "black"),
          axis.text.x = element_text(size = 17, colour = "black"),
          axis.text.y = element_text(size = 17, colour = "black"))
  
  
  #dev.new()  
  plot(plot)
  dev.off()
}