## Top 20 from the 2nd heatmap in TCGA 

# Correlation of IC and DELncRNAs 
# Prepare environment 
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/")

# Install required packages 
source(paste0("./R scripts/ipak.function.R"))
required.packages = c("scales","GSVA","gclus","scales","readxl","tidyr","ComplexHeatmap","corrplot","reshape2","dplyr","circlize")
ibiopak(required.packages)

# Loading the data 
load(("./Processed Data/BRCA/RNASeqData/BRCA_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))
# Check which patients included in the LNcRNA analysis 
colnames(filtered.norm.RNAseqData) = substr(colnames(filtered.norm.RNAseqData),1,12)
filtered.norm.RNAseqData = log(filtered.norm.RNAseqData+1,2)

# Load the 30 IC list 
IC_list = read_excel("~/Analysis/TCGA_differential_lncRNA_and_ICG.xlsx")
IC_list = unique(IC_list$gene)
filtered.norm.RNAseqData = filtered.norm.RNAseqData[which(rownames(filtered.norm.RNAseqData) %in% IC_list),]

# Load the diff expressed LNC  
matrix = read.csv("./data/lncRNA_TCGA_patients_QN.csv")
diff.genes = read.csv("./Analysis/DEG/diff_lncRNAs_TCGA.csv")

diff.genes = diff.genes[which(diff.genes$adj.P.Val < 0.05),]
diff.genes= diff.genes[order(diff.genes$logFC,decreasing = TRUE),]   # Ordr by FC 

#diff.genes = diff.genes[which(diff.genes$logFC > 1),]
matrix = matrix[which(rownames(matrix) %in% diff.genes$X.1),]
colnames(matrix) = gsub("\\.","-",colnames(matrix))

# Check names in both matrices 
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(colnames(filtered.norm.RNAseqData) %in% colnames(matrix))]
matrix = matrix[,which(colnames(matrix) %in% colnames(filtered.norm.RNAseqData))] 
matrix = as.matrix(matrix)
filtered.norm.RNAseqData = filtered.norm.RNAseqData[,colnames(matrix)]
Correlation = cor(t(matrix),t(filtered.norm.RNAseqData),method = "spearman") # With compliments to @user20650  

#do a little loop
dim(Correlation)
cor.matrix=Correlation
cor.matrix[cor.matrix<10000]=NA
p.matrix=cor.matrix
for (i in rownames(filtered.norm.RNAseqData)) {
  # print (i)
  for (j in rownames(matrix)){
    # print(j)
    test = cor.test(filtered.norm.RNAseqData[i,],matrix[j,],method = "spearman")
    cor.matrix[j,i] = test$estimate
    p.matrix[j,i] = test$p.value
  }
}

### Prepare for the heatmap

# Convert the matric to 3 columns dataframe 
data.cor  = setNames(melt(cor.matrix), c('Lnc', 'ICG', 'Cor'))
data.p  = setNames(melt(p.matrix), c('Lnc', 'ICG', 'p.value'))

# Check the order 
data.cor$Lnc == data.p$Lnc
data.cor$ICG == data.p$ICG

# Match the 2 matrices 
data.cor$ass = NA
data.cor$ass = paste0(data.cor$Lnc,"_",data.cor$ICG)
data.p$ass = NA
data.p$ass = paste0(data.p$Lnc,"_",data.p$ICG)
res <- left_join(data.cor,data.p,by="ass")
res$Cor[res$p.value > 0.05] = 0

TCGA_df = data.frame(Lnc = res$Lnc.x, gene = res$ICG.x, cor = res$Cor)
TCGA_matrix = reshape(TCGA_df, idvar="gene", timevar="Lnc", direction="wide")

rownames(TCGA_matrix) = TCGA_matrix$gene
TCGA_matrix$gene = NULL
colnames(TCGA_matrix) = substring(colnames(TCGA_matrix),5)

# Order the columns by col sums 
X1 = TCGA_matrix[,order(colSums(-TCGA_matrix),colSums(-TCGA_matrix,na.rm=TRUE))] ## tiebreaker orders NA columns properly

# Order the rows by absolute row sums 
X2 = X1[order(-rowSums(abs(X1),-rowSums(abs(X1),na.rm=TRUE))),] 

X2 = as.matrix(X2)

# Draw the heatmap 
svg(file = paste0("./Figures/Correlations/TCGA.IC.heatmap.cut.0.05.grey.ordered.with.NA.grey.spearman.svg"), 
    width = 20, height = 12,pointsize = 25)

colors = colorRamp2(c(min(X2),0, max(X2)), c("#2B71AB","white","#D41524"))
Heatmap(X2,cluster_rows = FALSE, cluster_columns = FALSE,row_names_gp = gpar(fontsize = 30),show_column_names = FALSE,col = colors)

dev.off()
