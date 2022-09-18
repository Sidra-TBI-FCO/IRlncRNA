# Prepare environment 
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/")

# Install required packages 
source(paste0("./R scripts/ipak.function.R"))
required.packages = c("scales","GSVA","gclus","scales","readxl","tidyr","corrplot","ComplexHeatmap","reshape2","dplyr")
ibiopak(required.packages)

## Calculate the Enrichment scores of TCGA 

common.3.pos.up.corr = c("ENSG00000247774.2", "ENSG00000256039.1", "ENSG00000235576.1")

Gene.set = "mutivar.3.pos.up.corr"

Cancertypes = c("BRCA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")


# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=4
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
 
  load(paste0("./Validation.dataset/Proceesed_data/",Cancer,".quantile.normalized.matrix.Rdata"))
  
  Expression.data = expression.matrix
  available_genes = rownames(Expression.data)
  
  Expression.data = as.matrix(Expression.data)
  Gene.list = mutivar.3.pos.up.corr # DElncRNA
 
  Expression.data = Expression.data[which(rownames(Expression.data) %in% Gene.list),]
  
  ES = Expression.data
  colnames(ES) = gsub("\\.","-",colnames(ES))
  
  ESs = as.matrix(ES)
  
  ## Match the samples with the samples with mRNA data from the ICR clustering table 
  clustering =  read.csv("../RNAseq-TCGA_BRCA/Analysis/ICRscores/Table_ICR_scores.all.TCGA.csv")
  
  clustering$TCGA_Barcode = substr(clustering$TCGA_Barcode,1,12)
  
  data = as.data.frame(t(ES))
  
  data = data[which(rownames(data) %in% clustering$TCGA_Barcode),]
  
  dir.create("./Validation.dataset/Analysis/signatures/",showWarnings = FALSE)
  
  save(ES,data, 
       file = paste0("./Validation.dataset/Analysis/signatures/ssGSEA",Cancer,".",Gene.set,".in.TCGA.Rdata"))
}


### Calculate the prognostic value of this signature 

# Install required packages 
source(paste0("~/Sidra Medicine - Research Division/TBI-LAB - General/Bioinformatics tools/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

# Parameters
Gene.set = "mutivar.3.pos.up.corr"

Cancertypes = c(
  "BRCA","HNSC","SKCM", "UCEC",  "LIHC", "STAD", "BLCA",  "CESC", "KICH", "OV", "LUSC", "READ", "COAD", "LUAD", "GBM","KIRP","KIRC","LGG")

# Empty object to store HR data 
results_df = data.frame(Cancertype = NA, HR = NA, p_val = NA, CI_lower = NA, CI_upper = NA,Gene = NA, log.rank = NA)

# Lenght of cancertypes 
N.Cancertypes = length(Cancertypes)

i=1
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  load(paste0("./Validation.dataset/Analysis/signatures/ssGSEA",Cancer,".",Gene.set,".in.TCGA.Rdata"))

    # Set parameters
  Surv.cutoff.years = 10
  Gene.set = "mutivar.3.pos.up.corr"
  
  # Load clinical data
  load("./Clinical Data/clinical_data.Rdata")
  clinical_data = clinical_data[which(clinical_data$bcr_patient_barcode %in% rownames(data)),]
  clinical_data = clinical_data[!duplicated(clinical_data$bcr_patient_barcode),]
  rownames(clinical_data) = clinical_data$bcr_patient_barcode
  ES = data
  clinical_data = merge(clinical_data, ES, by = "row.names")
  colnames(clinical_data)[which(colnames(clinical_data) %in% "ENSG00000235576.1")] = "Gene1"    #   AC092580.4
  colnames(clinical_data)[which(colnames(clinical_data) %in% "ENSG00000247774.2")] = "Gene2"    #   PCED1B-AS1
  colnames(clinical_data)[which(colnames(clinical_data) %in% "ENSG00000256039.1")] = "Gene3"    #   RP11-291B21.2
  
  Y = Surv.cutoff.years * 365
  clinical_data = clinical_data[!is.na(clinical_data$OS.time),]
  TS.Alive = clinical_data[clinical_data$vital_status == "Alive", c("vital_status", "OS.time", "Gene1", "Gene2" ,"Gene3")]
  colnames(TS.Alive) = c("Status","Time", "Gene1", "Gene2" ,"Gene3")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$vital_status == "Dead", c("vital_status", "OS.time", "Gene1", "Gene2" ,"Gene3")]
  colnames(TS.Dead) = c("Status","Time", "Gene1", "Gene2" ,"Gene3")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  
  multivariate_rev = coxph(formula = Surv(Time, Status) ~ Gene1 + Gene2 + Gene3, data = TS.Surv)
  summary(multivariate_rev)
  
  # To extract the coxp table 
  res.cox.sum <- summary(multivariate_rev)$coefficients
  res.cox.sum = as.data.frame(res.cox.sum)
  summary = summary(multivariate_rev)$conf.int
  summary = as.data.frame(summary)
  res.cox.sum$CI_lower = summary$`lower .95`
  res.cox.sum$CI_upper = summary$`upper .95`
  res.cox.sum$log.rank = summary(multivariate_rev)$sctest[3]
  
  HR_table = data.frame(HR = res.cox.sum$`exp(coef)`,p_val = res.cox.sum$`Pr(>|z|)`,CI_lower = res.cox.sum$CI_lower, CI_upper = res.cox.sum$CI_upper,log.rank = res.cox.sum$log.rank)
  HR_table$Gene = rownames(res.cox.sum)
  HR_table$Cancertype = Cancer
  
  results_df = rbind(results_df, HR_table)
}


colnames(results_df)[which(colnames(results_df) %in% "Cancertype")] = "Cancer"
colnames(results_df)[which(colnames(results_df) %in% "log.rank")] = "p_value"
results_df$gene.name = NA
results_df$gene.name[which(results_df$Gene %in% "Gene1")] = "AC092580.4"
results_df$gene.name[which(results_df$Gene %in% "Gene2")] = "PCED1B-AS1"
results_df$gene.name[which(results_df$Gene %in% "Gene3")] = "RP11-291B21.2"


write.csv(results_df, file = paste0("./Validation.dataset/Analysis/survival_analysis/",Gene.set,"_HR_all.cancer_results.csv"),row.names = FALSE)

### Calculate the Enrichment score of RAQA cohort

# Check the ENS ID and take corresponding Lnc ID
diff.genes = read.csv("./Analysis/DEG/diff_lncRNAs_TCGA.csv")


# Parameters 
Cancer = "RAQA"
Gene.set = "mutivar.3.pos.up.corr"
mutivar.3.pos.up.corr = c("AC092580.4","PCED1B-AS1","RP11-291B21.2") # common.3.pos.up.corr = list(c("ENSG00000247774.2", "ENSG00000256039.1", "ENSG00000235576.1"))
Gene.set = "mutivar.3.pos.up.corr"

# Loading the expression matric of Arabs 
load("./data/lncRNA.QN.patients.Rdata")

Expression.data = expression.matrix.QN.LOG2
available_genes = rownames(Expression.data)

Expression.data = as.matrix(Expression.data)
Gene.list = mutivar.3.pos.up.corr # DElncRNA

Expression.data = Expression.data[which(rownames(Expression.data) %in% Gene.list),]

ES = Expression.data

load("./data/CTA_TNBC_BRCA_Qatar_RNASeq.tumor.Clustering.RData")

clustering$Sample_ID = substr(clustering$Sample_ID,3,7)
ES = ES[,which(colnames(ES) %in% clustering$Sample_ID)]

data = as.data.frame(t(ES))

save(data, file = paste0("./Validation.dataset/Analysis/signatures/ssGSEA_",Cancer,".",Gene.set,".in.TCGA.spearman.matching.Rdata"))


### Survival analysis 

# Install packages and load
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
colnames(data)
Surv.cutoff.years = 10
Gene.set = "mutivar.3.pos.up.corr"
Cancer = "RAQA"

# Load clinical data
load("./data/lncRNA.QN.patients.Rdata")
load(paste0("./Validation.dataset/Analysis/signatures/ssGSEA_",Cancer,".",Gene.set,".in.TCGA.spearman.matching.Rdata"))

# Load Enrichment scores
load(paste0("./Analysis/ssGSEA/Scaled.ssGSEA.TACGA.",Gene.set,".in.RAQA.Rdata"))

annotation$days_to_last_follow_up <- as.Date(as.character(annotation$Last.seen.date..updated.), format="%m/%d/%Y")-
  as.Date(as.character(annotation$Diagnosis.Date), format="%m/%d/%Y")

clinical_data = annotation
rownames(clinical_data) = clinical_data$Coding
ES = data

clinical_data = merge(clinical_data, ES, by = "row.names")

colnames(clinical_data)[which(colnames(clinical_data) %in% "AC092580.4")] = "Gene1"    #   AC092580.4
colnames(clinical_data)[which(colnames(clinical_data) %in% "PCED1B-AS1")] = "Gene2"    #   PCED1B-AS1
colnames(clinical_data)[which(colnames(clinical_data) %in% "RP11-291B21.2")] = "Gene3"    #   RP11-291B21.2

# Empty object to store HR data 
results_df = data.frame(Cancer = NA, HR = NA, p_val = NA, CI_lower = NA, CI_upper = NA,Gene = NA, log.rank = NA)

Y = Surv.cutoff.years * 365
# time / event object creation
clinical_data = clinical_data[!is.na(clinical_data$days_to_last_follow_up),]
TS.Alive = clinical_data[clinical_data$Status == "Alive", c("Status", "days_to_last_follow_up", "Gene1", "Gene2" ,"Gene3")]
colnames(TS.Alive) = c("Status","Time", "Gene1", "Gene2" ,"Gene3")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = clinical_data[clinical_data$Status == "Dead", c("Status", "days_to_last_follow_up", "Gene1", "Gene2" ,"Gene3")]
colnames(TS.Dead) = c("Status","Time", "Gene1", "Gene2" ,"Gene3")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

multivariate_rev = coxph(formula = Surv(Time, Status) ~ Gene1 + Gene2 + Gene3, data = TS.Surv)
summary(multivariate_rev)

# To extract the coxp table 
res.cox.sum <- summary(multivariate_rev)$coefficients
res.cox.sum = as.data.frame(res.cox.sum)
summary = summary(multivariate_rev)$conf.int
summary = as.data.frame(summary)
res.cox.sum$CI_lower = summary$`lower .95`
res.cox.sum$CI_upper = summary$`upper .95`
res.cox.sum$log.rank = summary(multivariate_rev)$sctest[3]

HR_table = data.frame(HR = res.cox.sum$`exp(coef)`,p_val = res.cox.sum$`Pr(>|z|)`,CI_lower = res.cox.sum$CI_lower, CI_upper = res.cox.sum$CI_upper,log.rank = res.cox.sum$log.rank)
HR_table$Gene = rownames(res.cox.sum)
HR_table$Cancer = Cancer

results_df = rbind(results_df, HR_table)

colnames(results_df)[which(colnames(results_df) %in% "log.rank")] = "p_value"
results_df$gene.name = NA
results_df$gene.name[which(results_df$Gene %in% "Gene1")] = "AC092580.4"
results_df$gene.name[which(results_df$Gene %in% "Gene2")] = "PCED1B-AS1"
results_df$gene.name[which(results_df$Gene %in% "Gene3")] = "RP11-291B21.2"


write.csv(results_df, file = paste0("./Validation.dataset/Analysis/survival_analysis/",Gene.set,"_HR_",Cancer,"_results.csv"),row.names = FALSE)
