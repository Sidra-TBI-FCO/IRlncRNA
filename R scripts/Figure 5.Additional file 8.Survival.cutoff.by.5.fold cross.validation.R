rm(list = ls())

# Install required packages 
source(paste0("./R scripts/ipak.function.R"))

library(survivalROC)
library(risksetROC)
library(ROCR)
library(data.table)
library(ggplot2)
library(survminer)
library(tidyverse)
library(glmnet)
library(MLSurvival)
library(doParallel)
library(pracma)

ipak("caret")
registerDoParallel(cores=6)
setwd("./Validation.dataset")

load("./data.file.ICR.ssgsea.scaled.all.signatures.Rdata")
lncrna_df = all_scores_df

lncrna_df <- lncrna_df[!is.na(lncrna_df$OS.status) & !is.na(lncrna_df$OS.time),]

cancers <- unique(lncrna_df$Cancer)

Y = 3650

TS.Alive = lncrna_df[lncrna_df$OS.status=="Alive",]
TS.Alive$OS.time[TS.Alive$OS.time > Y] = Y

TS.Dead = lncrna_df[lncrna_df$OS.status=="Dead",]
TS.Dead$OS.status[which(TS.Dead$OS.time> Y)] = "Alive"
TS.Dead$OS.time[TS.Dead$OS.time > Y] = Y

lncrna_df = rbind(TS.Alive,TS.Dead)

lncrna_df$status <- 1
lncrna_df[lncrna_df$OS.status=="Alive",]$status <- 0

lncrna_df = subset(lncrna_df,lncrna_df$OS.time > 1) 

subset_lncrna_df_BRCA <- lncrna_df[lncrna_df$Cancer=="BRCA",]

#Have to remove as in READ when performing cross-validation, there are training set with all samples belonging to high/low categories and hence model will not run
cancers <- setdiff(cancers,c("READ","COAD"))

lncrna_df$time <- lncrna_df$OS.time/365
for (cancer in cancers)
{
  subset_lncrna_df <- lncrna_df[lncrna_df$Cancer==cancer,]
  
  #Get total samples
  set.seed(123)
  N <- nrow(subset_lncrna_df)
  subset_lncrna_df <- subset_lncrna_df[randperm(N),]
  samples <- c(1:N)
  
  #Get the fold sample ids for 5-fold cross validation
  flds <- createFolds(y=subset_lncrna_df$SSGSEA_Score, k = 5, list = TRUE, returnTrain = T)
  
  #Choose threshold to be between 25-75 percentile of the score distribution
  summary_scores <- as.numeric(quantile(subset_lncrna_df$SSGSEA_Score))
  min_t <- signif(summary_scores[2],digits=1)
  max_t <- signif(summary_scores[4],digits=1)
  
  opt_cutoffs <- seq(min_t,max_t,0.01)
  cv_mean_list <- list()
  
  #Find optimal cutoff based on concordance index
  for (i in 1:length(opt_cutoffs))
  {
    opt_cutoff <- opt_cutoffs[i]
    
    res_cat <- subset_lncrna_df[,c("time","status","SSGSEA_Score")]
    res_cat$SSGSEA_Category <- "low"
    if (sum(res_cat$SSGSEA_Score>opt_cutoff)>0)
    {
      res_cat[res_cat$SSGSEA_Score>opt_cutoff,]$SSGSEA_Category <- "high"
      
      #Get cv info for this cut-off through 5-fold cross-validation
      cv_info <- NULL
      for (j in 1:length(flds))
      {
        train_ids <- flds[[j]]
        test_ids <- setdiff(samples, train_ids)
        train_set <- res_cat[train_ids,]
        test_set <- res_cat[test_ids,]
        
        fit <- coxph(formula = Surv(time,status)~SSGSEA_Category, data = train_set)
        
        #Calculate CI
        y_predict <- predict(fit, newdata=test_set, type="lp")
        y_test <- as.matrix(dplyr::select(test_set, all_of(c("time","status"))))
        temp_ci <- Cindex(y_predict,y_test) 
        cv_info <- c(cv_info,temp_ci)
      }
      mean_cv <- mean(cv_info, na.rm=T)
    }else{
      mean_cv <- 0
    }
    cv_mean_list[i] <- mean_cv
  }
  cv_mean_list <- unlist(cv_mean_list)
  # plot(opt_cutoffs, cv_mean_list)
  
  #Final optimal cutoff
  order_ids <- order(cv_mean_list,decreasing = T)
  cv_mean_list < cv_mean_list[order_ids]
  opt_cutoffs <- opt_cutoffs[order_ids]
  
  #Select that cutoff for which the cv mean is maximum and the number of samples with high and low categories is at least 20% of population each
  for (i in 1:length(opt_cutoffs))
  {
    #Make the final survival model
    res.cat <- subset_lncrna_df
    res.cat$SSGSEA_Category <- "low"
    if (sum(res.cat$SSGSEA_Score>opt_cutoffs[i])>0)
    {
      res.cat[res.cat$SSGSEA_Score>opt_cutoffs[i],]$SSGSEA_Category <- "high"
    }
    N_low <- sum(res.cat$SSGSEA_Category=="low")
    N_high <- sum(res.cat$SSGSEA_Category=="high")
    if (N_low > 0.2*N & N_high > 0.2*N)
    {
      best_cutoff <- opt_cutoffs[i]
      break;
    }
  }
  
  #When we break the res.cat is corresponding to optimal cutoff
  fit <- survfit(Surv(time, status) ~ SSGSEA_Category, data = res.cat)
  dir.create("./Figures/Survival/",showWarnings = FALSE)
  svg(paste0("./Figures/Survival/",cancer,"_KM_best_cutoff.svg"),height=6, width=7)
 # dev.new()
  
  survp <- ggsurvplot(fit, data = res.cat, conf.int = F, pval = T,
                      pval.size = 8,
                      #     surv.median.line = "v",
                      break.time.by = 1,
                      fontsize = 7,
                      censor.shape = 3,
                      censor.size = 1.5,
                      censor = TRUE,
                      risk.table = TRUE,
                      tables.y.text.col = TRUE,
                      tables.y.text = FALSE,
                      tables.height = 0.2,
                      tables.theme = theme_cleantable(),
                      font.x =  22,font.y = 22,font.tickslab = 20,
                      title=paste0("Survival Analysis for Cancer: ",cancer)
  )
  print(survp)
  dev.off()
}

