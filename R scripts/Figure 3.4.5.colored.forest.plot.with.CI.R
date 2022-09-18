#clean global environment
rm(list=ls()) 

#Set working directory 
setwd(paste0("~/Validation.dataset/"))
source(paste0("~/R scripts/ipak.function.R"))

ipak("forestplot")

Gene.set = "ICRscore"    #   "scaled_ssGSEA" ,  "ICR_ES_scaled" , "Combin_ES_ICR_3LNC"  "SSGSEA

# Load data
t.test_results = read.csv(paste0("./Analysis/survival_analysis/",Gene.set,"_HR_all.cancer_results.csv"), stringsAsFactors = FALSE,sep=",")

colnames(t.test_results)

#remove digits from the end of value and meand_difference using (signif)
t.test_results$p_value = signif(t.test_results$p_value, digits=4)
t.test_results$HR = signif(t.test_results$HR, digits=4)
t.test_results$CI_lower = signif(t.test_results$CI_lower, digits=4)
t.test_results$CI_upper = signif(t.test_results$CI_upper, digits=4)

row.names(t.test_results) = t.test_results$Cancertype

# set t.test_results$Cancer as length 
N.Cancer = length(t.test_results$Cancertype)

# N.Cancer + 2
x = N.Cancer + 1

class(t.test_results)
t.test_results$CI = NA
t.test_results$CI = paste0(t.test_results$HR," [",t.test_results$CI_lower,"-",t.test_results$CI_upper,"]")[c(1:N.Cancer)]

cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,t.test_results$HR[1:N.Cancer]), NA),
    lower = c(NA,t.test_results$CI_lower[c(1:N.Cancer)], NA),
    upper = c(NA,t.test_results$CI_upper[c(1:N.Cancer)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

#labels and columns to appear in the plot
colnames(t.test_results)
table.text = cbind(
  c("Cancer", as.character(t.test_results$Cancer)[c(1:N.Cancer)]),
  c("N samples", as.character(t.test_results$`N_Samples`)[c(1:N.Cancer)]),
  c("p value", t.test_results$p_value[c(1:N.Cancer)]),
  c("HR", t.test_results$CI[c(1:N.Cancer)]))

col_lines = list(gpar(col = "black"))
col_boxes = list(gpar(fill = "black"))

tabletext = table.text

pdf(file = paste0("./Figures/Survival/Forest.",Gene.set,".signature.tumors.with.CI.N.samples.pdf"), height = 6, width = 10)

i = 1
for (i in c(2:x)) {
  if (is.na(tabletext[i,3])){
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "black",col = "black"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "black"))))
    next
  }
  if (tabletext[i,3] == ""){
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "black",col = "black"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "black"))))
    next
  }
  if (as.numeric(tabletext[i,3]) < 0.05 & cochrane_from_rmeta$mean[i] > 1) {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#05AEE9",col = "#05AEE9"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#05AEE9"))))
  } else if (as.numeric(tabletext[i,3]) < 0.05 & cochrane_from_rmeta$mean[i] < 1) {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#F83C5D",col = "#F83C5D"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#F83C5D"))))
  } else {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#8C9599",col = "#8C9599"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#8C9599"))))
  }
} 

styles <- fpShapesGp(lines = col_lines,box = col_boxes)

forestplot(tabletext, 
           cochrane_from_rmeta,
           is.summary=c(TRUE,rep(FALSE,x-2),FALSE),
           #  hrzl_lines = gpar(col="#444444"),
           clip=c(0.7,10), 
           xlim=c(-8,20),
           xlog=FALSE, 
           zero = 1,
           #  ci.vertices = TRUE,
           new_page = FALSE,
           col=fpColors(box="darkblue",line="black"),
           shapes_gp = styles,
           align = c("l", "l", "l"), 
           boxsize = .2,
           vertices = FALSE,
           xticks = c(0,1,2,3),
           xlab = paste0("HR"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 13), xlab = gpar(fontsize = 17),
                            ticks = gpar(fontsize = 17)))

dev.off()

### BRCA only forest plot 

t.test_results = read.csv(paste0("./Analysis/BRCA signatures.csv"), stringsAsFactors = FALSE,sep=",")
t.test_results$N_Samples = NULL

colnames(t.test_results)[which(colnames(t.test_results) %in% "Cancertype")] = "Signature"

#remove digits from the end of value and meand_difference using (signif)
t.test_results$p_val = signif(t.test_results$p_val, digits=4)
t.test_results$HR = signif(t.test_results$HR, digits=4)
t.test_results$CI_lower = signif(t.test_results$CI_lower, digits=4)
t.test_results$CI_upper = signif(t.test_results$CI_upper, digits=4)

row.names(t.test_results) = t.test_results$Signature

# set t.test_results$Cancer as length 
N.Signature = length(t.test_results$Signature)

# N.Cancer + 2
x = N.Signature + 1

class(t.test_results)
t.test_results$CI = NA
t.test_results$CI = paste0(t.test_results$HR," [",t.test_results$CI_lower,"-",t.test_results$CI_upper,"]")

#x = N.Cancer
# Cochrane data from the 'rmeta'-package, (must use cochrane??)
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,t.test_results$HR[1:N.Signature]), NA),
    lower = c(NA,t.test_results$CI_lower[c(1:N.Signature)], NA),
    upper = c(NA,t.test_results$CI_upper[c(1:N.Signature)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

#labels and columns to appear in the plot

colnames(t.test_results)
table.text = cbind(
  c("Signature", as.character(t.test_results$Signature)[c(1:N.Signature)]),
  c("p value", t.test_results$p_val[c(1:N.Signature)]),
  c("HR", t.test_results$CI[c(1:N.Signature)]))

col_lines = list(gpar(col = "black"))
col_boxes = list(gpar(fill = "black"))


tabletext = table.text

pdf(file = paste0("./Figures/Survival/Forest.signatures.BRCA.with.CI.N.samples.pdf"), height = 2, width = 6)

i = 1
for (i in c(2:x)) {
  if (is.na(tabletext[i,2])){
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "black",col = "black"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "black"))))
    next
  }
  if (tabletext[i,2] == ""){
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "black",col = "black"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "black"))))
    next
  }
  if (as.numeric(tabletext[i,2]) < 0.05 & cochrane_from_rmeta$mean[i] > 1) {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#05AEE9",col = "#05AEE9"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#05AEE9"))))
  } else if (as.numeric(tabletext[i,2]) < 0.05 & cochrane_from_rmeta$mean[i] < 1) {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#F83C5D",col = "#F83C5D"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#F83C5D"))))
  } else {
    col_boxes = do.call(c, list(col_boxes, list(gpar(fill = "#8C9599",col = "#8C9599"))))
    col_lines = do.call(c, list(col_lines, list(gpar(col = "#8C9599"))))
  }
} 


styles <- fpShapesGp(lines = col_lines,box = col_boxes)

dev.new()
forestplot(tabletext, 
           cochrane_from_rmeta,
           is.summary=c(TRUE,rep(FALSE,x-2),FALSE),
           clip=c(0.7,10), 
           xlim=c(-8,20),
           xlog=FALSE, 
           zero = 1,
           new_page = FALSE,
           col=fpColors(box="darkblue",line="black"),
           shapes_gp = styles,
           align = c("l", "l", "l"), 
           boxsize = .1,
           vertices = FALSE,
           xticks = c(0,0.2,0.4,0.6,0.8,1),
           xlab = paste0("HR"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 8), xlab = gpar(fontsize = 15),
                            ticks = gpar(fontsize = 15)))

dev.off()
