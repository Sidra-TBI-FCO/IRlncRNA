# Calculate the ConsensusTME 

#clean global environment
rm(list=ls())

#Set working directory
load("~/R.Config.Rdata")
setwd("~/")
source(paste0(toolbox.path,"./R scripts/ipak.function.R"))

library(ConsensusTME)

# Set parameters
approach = "consensusTME_TBI"

Cancertypes = c("BRCA","UCEC", "OV",  "READ", "COAD","HNSC","SKCM", "LIHC", "STAD", "BLCA",  "CESC", "KICH", "LUSC",  "LUAD", "GBM","KIRP","KIRC","LGG")

# Parameters 
method = "ssgsea" 

N.Cancertypes = length(Cancertypes)

i=1
for (i in 1:N.Cancertypes) {
  cancer_type = Cancertypes[i]
  
  # Load Corresponding data
  load(paste0("./Processed.data/TCGA-", cancer_type, "/003_TCGA-", cancer_type, "_EDAseq_normalized.gene.info.2017.Rdata"))
  
  colnames(RNASeq.NORM.quantiles) = substr(colnames(RNASeq.NORM.quantiles),1,12)
  colnames(RNASeq.NORM.quantiles) = gsub("\\.","-",colnames(RNASeq.NORM.quantiles))
  RNASeq.NORM.quantiles = RNASeq.NORM.quantiles[,!duplicated(colnames(RNASeq.NORM.quantiles))]
  dim(RNASeq.NORM.quantiles)
  
  # Log2 transformation
  RNASeq.QN.counts.filtered = log(RNASeq.NORM.quantiles+1,2)
 
  dir.create(paste0("./Validation.dataset/Analysis/signatures/consensusTME_TBI/"),showWarnings = FALSE)
  
  if(approach == "consensusTME_TBI"){
    load(paste0("./GSEA list/immune.gene.lists.v4.Rdata"))
    source(paste0("./R scripts/consensusTMEAnalysis_WJed.R"))
    ES = consensusTMEAnalysis_WJed(RNASeq.QN.counts.filtered, cancer = cancer_type, statMethod = method, include_source = TRUE) # include_source = TRUE (when you want to see the source of the code behind it)
    save(ES, file = paste0("./Validation.dataset/Analysis/signatures/consensusTME_TBI/", approach, "_", cancer_type, "_", method, ".Rdata"))
  }
  
}

### Calculate the ConsensusTME Enrichment scores in RAQA

# Parameters
cancer_type = "BRCA"
method = "ssgsea"
  
  # Loading the expression matriX of Arabs 
  load("~/Dropbox (SMCP)/NNN-BRCA/RNAseq-RA-QA/Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.annotated.EDAseq.QN.Rdata")
  
  filtered.norm.RNAseqData = RNASeq.QN.LOG2
  
  dim(filtered.norm.RNAseqData)

  # Log2 transformation
  RNASeq.QN.counts.filtered = log(filtered.norm.RNAseqData+1,2)
  
  dir.create(paste0("./Validation.dataset/Analysis/signatures/consensusTME_TBI/"),showWarnings = FALSE)
  
  rownames(RNASeq.QN.counts.filtered) = substr(rownames(RNASeq.QN.counts.filtered),3,8)
  
  rownames(RNASeq.QN.counts.filtered) <- sub("^0+", "", rownames(RNASeq.QN.counts.filtered))    
  
  if(approach == "consensusTME_TBI"){
    load(paste0("~/GSEA list/immune.gene.lists.v4.Rdata"))
    source(paste0("~/R scripts/consensusTMEAnalysis_WJed.R"))
    ES = consensusTMEAnalysis_WJed(RNASeq.QN.counts.filtered, cancer = cancer_type, statMethod = method, include_source = TRUE)
    save(ES, file = paste0("./Validation.dataset/Analysis/signatures/consensusTME_TBI/", approach, "_RAQA_", method, ".Rdata"))
  }
  