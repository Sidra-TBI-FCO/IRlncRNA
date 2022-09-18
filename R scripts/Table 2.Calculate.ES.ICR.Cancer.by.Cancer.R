
# Prepare environment 
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/")

# Install required packages 
source(paste0("~/R scripts/ipak.function.R"))
required.packages = c("scales","GSVA","gclus","readxl","tidyr","corrplot","ComplexHeatmap","reshape2","dplyr")
ibiopak(required.packages)

# Loading the data 
load(paste0("./GSEA list/immune.gene.lists.v4.2.Rdata"))


Gene.set = "ICR_list"
Cancertypes = c("BRCA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")

# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=3
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  
  # Loading the expression matrix of TCGA 
  
  load(paste0("~/RNAseq-TCGA_BRCA/Processed.data/TCGA-",Cancer,"/003_TCGA-",Cancer,"_EDAseq_normalized.gene.info.2017.Rdata"))
  RNASeq.NORM.quantiles["CCL5",]
  
  dim(RNASeq.NORM.quantiles)
  
  Expression.data = log(RNASeq.NORM.quantiles +1, 2)
  available_genes = rownames(Expression.data)
  
  Expression.data = as.matrix(Expression.data)
  Gene.list = ICR_list # DElncRNA
  unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]
  
  cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
             " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
             " genes are available in expression data."), append = TRUE, sep = "\n")
  
  ## ssGSEA
  ES = gsva(Expression.data,Gene.list,method="ssgsea")
  
  colnames(ES) = substr(colnames(ES),1,12)
  colnames(ES) = gsub("\\.","-",colnames(ES))
  ES = ES[,!duplicated(colnames(ES))]
  dim(ES)
  
  ESz = ES 
  for(j in 1: nrow(ESz))  {
    ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
  }
  
  
  ESs = as.matrix(ES)
  
  # Sum them together if we want to use the collective score 
  data = data.frame(ICR_ES = ESs[1,], ICR_ES_scaled = rescale(ESs[1,],to = c(1,10)))
  
  dir.create("./Validation.dataset/Analysis/signatures",showWarnings = FALSE)
  dir.create(paste0("./Validation.dataset/Analysis/signatures/",Gene.set),showWarnings = FALSE)
  
  save(ES,data, 
       file = paste0("./Validation.dataset/Analysis/signatures/",Gene.set,"/",Gene.set,"_",Cancer,".Rdata"))
}


## Calculate the ICR Enrichment score in RAQA

# Parameters 
Cancer = "RAQA"
Gene.set = "ICR_list"

# Loading the expression matric of Arabs 
load("./Data/expression matrix/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.tissue.annotated.EDAseq.QN.Rdata")


Expression.data = RNASeq.QN.LOG2
available_genes = rownames(Expression.data)

Gene.list = ICR_list
unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]

cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
           " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
           " genes are available in expression data."), append = TRUE, sep = "\n")


Expression.data = as.matrix(Expression.data)

## ssGSEA
ES = gsva(Expression.data,Gene.list,method="ssgsea")
ESz = ES 
for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
}

ESs = as.matrix(ES)

# Sum them together if we want to use the collective score 
data = data.frame(ICR_ES = ESs[1,], ICR_ES_scaled = rescale(ESs[1,],to = c(1,10)))
rownames(data) = substr(rownames(data),3,8)

rownames(data) <- sub("^0+", "", rownames(data))        # Apply sub function

dir.create("./Validation.dataset/Analysis/signatures",showWarnings = FALSE)
dir.create(paste0("./Validation.dataset/Analysis/signatures/",Gene.set),showWarnings = FALSE)

save(ES,data, 
     file = paste0("./Validation.dataset/Analysis/signatures/",Gene.set,"/",Gene.set,"_",Cancer,".Rdata"))
