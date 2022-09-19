## Calculate the enrichment scores of different signatures in TCGA and RAQA 

# Prepare environment 
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/Dropbox (SMCP)/NNN-BRCA/LncRNA_Analysis_Raghvendra")

# Install required packages 
source(paste0("./R scripts/ipak.function.R"))
required.packages = c("scales","GSVA","gclus","scales","readxl","tidyr","corrplot","ComplexHeatmap","reshape2","dplyr")
ibiopak(required.packages)

#### Calculate the Enrichment score in TCGA 

common.3.pos.up.corr = list(c("ENSG00000247774.2", "ENSG00000256039.1", "ENSG00000235576.1"))
Gene.set = "common.3.pos.up.corr"

Cancertypes = c("BRCA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")
# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=3
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  
  load(paste0("./Validation.dataset/Proceesed_data/",Cancer,".quantile.normalized.matrix.Rdata"))
  
  Expression.data = expression.matrix
  available_genes = rownames(Expression.data)
  
  Expression.data = as.matrix(Expression.data)
  Gene.list = common.3.pos.up.corr # DElncRNA
  unavailable_genes_RNAseq = unlist(Gene.list)[-which(unlist(Gene.list) %in% rownames(Expression.data))]
  
  cat(paste0(Gene.set," ssGSEA ", ". Total number of genes is ", length(unlist(Gene.list)), ".",
             " Of which ", length(unlist(Gene.list)[unlist(Gene.list) %in% available_genes]), 
             " genes are available in expression data."), append = TRUE, sep = "\n")
  
  ## ssGSEA
  ES = gsva(Expression.data,Gene.list,method="ssgsea")
  ESz = ES 
  for(j in 1: nrow(ESz))  {
    ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
  }
  
  colnames(ES) = gsub("\\.","-",colnames(ES))
  
  ESs = as.matrix(ES)
 
  # Sum them together if we want to use the collective score 
  data = data.frame(common.3.pos.up.corr = ESs[1,])
  
  dir.create("./Validation.dataset/Analysis/signatures/",showWarnings = FALSE)
  
  save(ES,data, 
       file = paste0("./Validation.dataset/Analysis/signatures/ssGSEA",Cancer,".",Gene.set,".in.TCGA.Rdata"))
}


### Calculate the prognostic value of this signature in TCGA

#Setup environment
rm(list = ls())
setwd("~/")

# Install required packages 
source(paste0("./R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

Gene.set = "common.3.pos.up.corr"
Cancertypes = c("BRCA","BRCA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")

# Empty object to store HR data 
results_df = data.frame(Cancertype = Cancertypes, HR = NA, p_val = NA, CI_lower = NA, CI_upper = NA)

# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=4
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  
  load(paste0("./Validation.dataset/Analysis/signatures/ssGSEA",Cancer,".",Gene.set,".in.TCGA.Rdata"))
  # Set parameters
  Surv.cutoff.years = 10
  Gene.set = "common.3.pos.up.corr"  #  "High.corr"     #Inv.corr      ICR.score     Sum.High.Inv
  
  # Load clinical data
  load("./Clinical Data/clinical_data.Rdata")
  
  clinical_data = clinical_data[which(clinical_data$bcr_patient_barcode %in% rownames(data)),]
  
  clinical_data = clinical_data[!duplicated(clinical_data$bcr_patient_barcode),]
  
  rownames(clinical_data) = clinical_data$bcr_patient_barcode
  
  ES = data
  
  clinical_data = merge(clinical_data, ES, by = "row.names")
  
  HR_table = data.frame(Signature = colnames(ES), p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)
  
  i=2
  for (i in 1:ncol(ES)){
    Group.of.interest = colnames(ES)[i]
    Y = Surv.cutoff.years * 365
    # time / event object creation
    clinical_data = clinical_data[!is.na(clinical_data$OS.time),]
    TS.Alive = clinical_data[clinical_data$vital_status == "Alive", c("vital_status", "OS.time", Group.of.interest)]
    colnames(TS.Alive) = c("Status","Time", Group.of.interest)
    TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
    TS.Alive$Time[TS.Alive$Time > Y] = Y
    
    TS.Dead = clinical_data[clinical_data$vital_status == "Dead", c("vital_status", "OS.time", Group.of.interest)]
    colnames(TS.Dead) = c("Status","Time", Group.of.interest)
    TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
    TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
    TS.Dead$Time[TS.Dead$Time > Y] = Y
    
    TS.Surv = rbind (TS.Dead,TS.Alive)
    TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
    TS.Surv$Status <- TS.Surv$Status == "Dead"
    TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
    
    uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
    summary = summary(uni_variate)
    HR = summary$conf.int[1]
    CI_lower = summary$conf.int[3]
    CI_upper = summary$conf.int[4]
    p_value = summary$coefficients[5]
    HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
    HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
    HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
    HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
    
    # fill out results_df object 
    results_df$HR[which(results_df$Cancertype == Cancer)] = HR
    results_df$p_val[which(results_df$Cancertype == Cancer)] = p_value
    results_df$CI_lower[which(results_df$Cancertype == Cancer)] = CI_lower
    results_df$CI_upper[which(results_df$Cancertype == Cancer)] = CI_upper
    
  }
}

colnames(results_df)[which(colnames(results_df) %in% "Cancertype")] = "Cancer"
colnames(results_df)[which(colnames(results_df) %in% "p_val")] = "p_value"

write.csv(results_df, file = paste0("./Validation.dataset/Analysis/survival_analysis/",Gene.set,"_HR_all.cancer_results.csv"),row.names = FALSE)


### RAQA cohort 
rm(list = ls())

Cancer = "RAQA"
Gene.set = "common.3.pos.up.corr"

# Check the ENS ID and take corresponding Lnc ID
diff.genes = read.csv("../Analysis/DEG/diff_lncRNAs_TCGA.Shimaa.reversed.csv")

# Parameters 
common.3.pos.up.corr = list(c("PCED1B-AS1","RP11-291B21.2","AC092580.4")) # common.3.pos.up.corr = list(c("ENSG00000247774.2", "ENSG00000256039.1", "ENSG00000235576.1"))
Gene.set = "common.3.pos.up.corr"

# Loading the expression matric of Arabs 
load("./data/lncRNA.QN.patients.Rdata")

Expression.data = expression.matrix.QN.LOG2

available_genes = rownames(Expression.data)
Gene.list = common.3.pos.up.corr # DElncRNA
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

ESs = ES

# Sum them together if we want to use the collective score 
data = data.frame(common.3.pos.up.corr = ESs[1,])
save(data, file = paste0("./Validation.dataset/Analysis/signatures/ssGSEA_",Cancer,".",Gene.set,".in.TCGA.spearman.matching.Rdata"))


### Calculate the prognostic value of this signature in RAQA
#Setup environment
load("~/R.Config.Rdata")
setwd(master.location)
setwd("/")
# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
colnames(data)
Surv.cutoff.years = 10
Gene.set = "common.3.pos.up.corr" #  "High.corr"     #Inv.corr      ICR.score     Sum.High.Inv

# Load clinical data
load("./data/lncRNA.QN.patients.Rdata")
annotation$days_to_last_follow_up <- as.Date(as.character(annotation$Last.seen.date..updated.), format="%m/%d/%Y")-
  as.Date(as.character(annotation$Diagnosis.Date), format="%m/%d/%Y")

clinical_data = annotation
rownames(clinical_data) = clinical_data$Coding

clinical_data = merge(clinical_data, ES, by = "row.names")

HR_table = data.frame(Signature = colnames(ES), p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)

i=2
for (i in 1:ncol(ES)){
  Group.of.interest = colnames(ES)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  clinical_data = clinical_data[!is.na(clinical_data$days_to_last_follow_up),]
  TS.Alive = clinical_data[clinical_data$Status == "Alive", c("Status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$Status == "Dead", c("Status", "days_to_last_follow_up", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
  HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
  HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
  HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
}