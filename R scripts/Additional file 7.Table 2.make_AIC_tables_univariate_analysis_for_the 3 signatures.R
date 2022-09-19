# Script to create the Supplementary AIC table

#Setup environment
rm(list = ls())
setwd("~/")

# Install required packages 
source(paste0("~/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival","RColorBrewer", "forestplot")
ipak(required.packages)

library(survivalROC)
library(risksetROC)
library(ROCR)
library(data.table)
library(ggplot2)
library(survminer)
library(tidyverse)
library(timeROC)
library(pracma)
library(pROC)
library(survivalAnalysis)


#Gene.set = "common.3.pos.up.corr"
load("./data.file.ICR.ssgsea.scaled.all.signatures.Rdata")
all_scores_df <- all_scores_df[!is.na(all_scores_df$OS.status) & !is.na(all_scores_df$OS.time),]
Cancertypes <- unique(all_scores_df$Cancer)

colnames(all_scores_df)

N.Cancertypes = length(Cancertypes)

# Set parameters
Surv.cutoff.years = 10

# Analysis
rownames(all_scores_df) = all_scores_df$X
colnames(all_scores_df)

output_df <- NULL

i=2
for (i in 1:N.Cancertypes) {
  Cancer = Cancertypes[i]
  
  subset_all_scores_df <- all_scores_df[all_scores_df$Cancer==Cancer,]
   
  Y = Surv.cutoff.years * 365
  # time / event object creation
  subset_all_scores_df = subset_all_scores_df[!is.na(subset_all_scores_df$OS.time),]
  TS.Alive = subset_all_scores_df[subset_all_scores_df$OS.status == "Alive", c("OS.status", "OS.time", "ICR_ES_scaled","scaled_ssGSEA","Combin_ES_ICR_3LNC")]
  colnames(TS.Alive) = c("Status","Time", "ICR_ES_scaled","scaled_ssGSEA","Combin_ES_ICR_3LNC")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = subset_all_scores_df[subset_all_scores_df$OS.status == "Dead", c("OS.status", "OS.time","ICR_ES_scaled","scaled_ssGSEA","Combin_ES_ICR_3LNC" )]
  colnames(TS.Dead) = c("Status","Time", "ICR_ES_scaled","scaled_ssGSEA","Combin_ES_ICR_3LNC")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  

  icr_fit = coxph(formula = Surv(Time, Status) ~ ICR_ES_scaled, data = TS.Surv)
  lncrna_fit = coxph(formula = Surv(Time, Status) ~ scaled_ssGSEA, data = TS.Surv)
  combin_fit = coxph(formula = Surv(Time, Status) ~ Combin_ES_ICR_3LNC, data = TS.Surv)
  

#Convert survival models into dataframe
icr_fit_df <- cox_as_data_frame(icr_fit)
lncrna_fit_df <- cox_as_data_frame(lncrna_fit)
combin_fit_df <- cox_as_data_frame(combin_fit)

#Get HR for ICR (continuous) in appropriate format
icr_hr <- signif(icr_fit_df$HR,digits=4)
icr_pval <- signif(icr_fit_df$p, digits=4)
icr_low <- signif(icr_fit_df$Lower_CI,digits=4)
icr_high <- signif(icr_fit_df$Upper_CI,digits = 4)
icr_hr_info <- paste0(icr_hr," [",icr_low,"-",icr_high,"]")

#Get HR for lncRNA (continuous) in appropriate format
lncrna_hr <- signif(lncrna_fit_df$HR,digits=4)
lncrna_pval <- signif(lncrna_fit_df$p, digits = 4)
lncrna_low <- signif(lncrna_fit_df$Lower_CI,digits=4)
lncrna_high <- signif(lncrna_fit_df$Upper_CI,digits = 4)
lncrna_hr_info <- paste0(lncrna_hr," [",lncrna_low,"-",lncrna_high,"]")


#Get HR for lncRNA (continuous) in appropriate format
combin_hr <- signif(combin_fit_df$HR,digits=4)
combin_pval <- signif(combin_fit_df$p, digits = 4)
combin_low <- signif(combin_fit_df$Lower_CI,digits=4)
combin_high <- signif(combin_fit_df$Upper_CI,digits = 4)
combin_hr_info <- paste0(combin_hr," [",combin_low,"-",combin_high,"]")

aic_icr <- extractAIC(icr_fit)
aic_lncrna <- extractAIC(lncrna_fit)
aic_combin <- extractAIC(combin_fit)


temp <- cbind(Cancer, nrow(TS.Surv),icr_low,icr_high , icr_pval,icr_hr,icr_hr_info, lncrna_low,lncrna_high, lncrna_pval,lncrna_hr,lncrna_hr_info, combin_low,combin_high, combin_pval,combin_hr,combin_hr_info ,signif(aic_icr[2],4), signif(aic_lncrna[2],4),  signif(aic_combin[2],4), signif(aic_icr[2]-aic_lncrna[2],4))
output_df <- rbind(output_df, temp)

}

output_df <- as.data.frame(output_df)

colnames(output_df)

colnames(output_df) <- c("Cancer","N_Samples","icr_low","icr_high", "ICR_Pval","ICR_HR","icr_hr_info","lncrna_low","lncrna_high","LncRNA_Pval", "LncRNA_HR","lncrna_hr_info","combin_low",  "combin_high", "combin_pval", "combin_hr" ,"combin_hr_info","ICR_AIC","LncRNA_AIC","Combine_AIC","dAIC")
output_df$N_Samples <- as.numeric(as.vector(output_df$N_Samples))
output_df$ICR_Pval <- as.numeric(as.vector(output_df$ICR_Pval))
output_df$LncRNA_Pval <- as.numeric(as.vector(output_df$LncRNA_Pval))
output_df$ICR_AIC <- as.numeric(as.vector(output_df$ICR_AIC))
output_df$LncRNA_AIC <- as.numeric(as.vector(output_df$LncRNA_AIC))
output_df$Combine_AIC <- as.numeric(as.vector(output_df$Combine_AIC))
output_df$dAIC <- as.numeric(as.vector(output_df$dAIC))

write.csv(output_df, file = paste0("./survival_analysis/ssGSEA.ICRscaled.Combine_HR_all.cancer_results.csv"),row.names = FALSE)
