# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(master.location)
setwd("~/")
source(paste0(toolbox.path,"/R scripts/ipak.function.R"))
required.packages <- c("corrplot", "stringr")
ipak(required.packages)    

# Set Parameters
colpattern = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 410)
selected_genes = "ConsensusTME"                                                                                                   # Specify which genes will be correlated
test = "pearson"
score = "SSGSEA_Score"

load("./review.data.file.ICR.ssgsea.scaled.all.signatures.Rdata")

all_scores_df = all_scores_df[-which(all_scores_df$Cancer %in% "RAQA"),]

cancers <- unique(all_scores_df$Cancer)


for (cancer in cancers)
{  
load(paste0("./Validation.dataset/Analysis/signatures/consensusTME_TBI/consensusTME_TBI_",cancer,"_ssgsea.Rdata"))

  annotation = all_scores_df
  
  All_cancers_data = ES

All_cancers_data = All_cancers_data[,!duplicated(colnames(All_cancers_data))]

dim(All_cancers_data)

colnames(All_cancers_data) = substr(colnames(All_cancers_data),1,12)

All_cancers_data = All_cancers_data[,which(colnames(All_cancers_data) %in% all_scores_df$X)]
All_cancers_data = as.data.frame(t(All_cancers_data))

All_cancers_data$lnc.3.sig = all_scores_df$SSGSEA_Score[match(rownames(All_cancers_data),all_scores_df$X)]

annotation = annotation[which(annotation$Cancer %in% cancer),]


Hallmark.enrichment.score.df = All_cancers_data
score = "SSGSEA_Score"
signature1 = "SSGSEA_Score"

Hallmark.enrichment.score.df[,score] = all_scores_df[,score][match(row.names(Hallmark.enrichment.score.df), all_scores_df$X)]

Hallmark_GSEA_cor <- cor (Hallmark.enrichment.score.df,method=test)

mean_correlation_table = data.frame(Cancers = cancer, Mean.correlation = 0)

cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
Hallmark_GSEA_cor_sign <- cor.mtest(Hallmark_GSEA_cor, 0.95)

dir.create(paste0("./Validation.dataset/Analysis/Correlations/"),showWarnings = FALSE)

save(Hallmark_GSEA_cor, Hallmark_GSEA_cor_sign, file = paste0("./Validation.dataset/Analysis/Correlations/Correlation_",score,"_",selected_genes,"_",cancer,".Rdata"))

}

              #######################     Part 2     #########################

required.packages <- c("corrplot", "stringr","ComplexHeatmap","ggplot2")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
my.palette = colorRampPalette(c("#161C52", "white", "#800000"))(n = 444)

selected_genes = "Selected.Pathways"                                                                                                   # Specify which genes will be correlated
test = "pearson"
Cancer = "Pancancer"
selected_genes = "ConsensusTME"                                                                                                   # Specify which genes will be correlated
display_correlations = "only_significant"
score = "SSGSEA_Score"

# Load data
load("./Validation.dataset/For.Ragh/review.data.file.ICR.ssgsea.scaled.all.signatures.Rdata")
load(paste0("./Validation.dataset/Analysis/Correlations/Correlation_",score,"_",selected_genes,"_LGG.Rdata"))
all_scores_df = all_scores_df[-which(all_scores_df$Cancer %in% "RAQA"),]

N.sets = length(unique(all_scores_df$Cancer))

cancersets = data.frame(unique(all_scores_df$Cancer))
pancancer_Geneset_cor_table = t(cancersets)
colnames(pancancer_Geneset_cor_table) = pancancer_Geneset_cor_table[c(1),]

pancancer_Geneset_cor_table = rbind(pancancer_Geneset_cor_table,matrix(nrow = nrow(Hallmark_GSEA_cor),ncol=N.sets))
pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[-c(1),]
rownames(pancancer_Geneset_cor_table) = rownames(Hallmark_GSEA_cor)
cancers = unique(all_scores_df$Cancer)
  
i=1
for (cancer in cancers) 
  {
  load(paste0("./Validation.dataset/Analysis/Correlations/Correlation_",score,"_",selected_genes,"_",cancer,".Rdata"))
  pancancer_Geneset_cor_table[, cancer] = as.numeric(Hallmark_GSEA_cor[,score])
}

# convert to numeric matrix
mode(pancancer_Geneset_cor_table) = "numeric"
pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[-which(rownames(pancancer_Geneset_cor_table) == score),]

#Correlation complex heatmap
pdf(paste0("./Validation.dataset/Figures/Correlations/Correlation_lnc.3.sig_",selected_genes,"_all_cancers.red.green.pdf"), width = 11, height = 12)

Geneset_cor = pancancer_Geneset_cor_table
cex.before <- par("cex")
par(cex = 0.35)
lims=c(-1,1)
if (length(Geneset_cor[Geneset_cor<0]) == 0) {lims=c(0,1)}
annotation = data.frame (Cancer = colnames(Geneset_cor),color = c("brown" , "#1ABC9C" ,"#4B0082","#CD6600","green","light blue","orange","cyan","#CD3278","purple", "black","red","pink","violet","dark blue","yellow","#8B6508","#40E0D0"),stringsAsFactors = FALSE)
Geneset_cor = as.data.frame(Geneset_cor)
annotation$Cancer = factor(annotation$Cancer, levels = c("BRCA" ,"HNSC", "SKCM" ,"UCEC", "LIHC" ,"STAD" ,"BLCA" ,"CESC", "KICH" ,"OV", "LUSC", "READ" ,"COAD","LUAD" ,"GBM"  ,"KIRP" ,"KIRC", "LGG"))
ha = HeatmapAnnotation(`Cancer` = annotation$Cancer,
                       col = list(`Cancer` = c("BRCA"="brown" ,"HNSC"="#1ABC9C", "SKCM"="#4B0082" ,"UCEC"="#CD6600", "LIHC"="green" ,"STAD"="light blue" ,"BLCA"="orange" ,"CESC"="cyan", "KICH"="#CD3278" ,"OV"="purple", "LUSC"="black", "READ" ="red","COAD"="pink","LUAD"="violet" ,"GBM"="dark blue"  ,"KIRP"="yellow" ,"KIRC"= "#8B6508", "LGG" = "#40E0D0")),annotation_name_gp = gpar(fontsize = 22, fontface = "bold"))

Geneset_cor = Geneset_cor[-c(35),]
Geneset_cor = Geneset_cor[-c(27),]

Geneset_cor = as.matrix(Geneset_cor)
HM = Heatmap(Geneset_cor,
             column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
             row_title_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 22),
             cluster_rows = TRUE ,cluster_columns = FALSE ,
             show_column_names = FALSE,top_annotation = ha,name = "value",col=my.palette,
             row_names_max_width = unit(9, "in")
)
draw(HM, heatmap_legend_side = "left", annotation_legend_side = "left")

par(cex = cex.before)
dev.off()

