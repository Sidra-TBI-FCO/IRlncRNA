
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/")
source(paste0("~/R scripts/ipak.function.R"))

required.packages = c("plyr", "ggplot2", "ggrepel", "limma", "AnnotationDbi","stringr","EnsDb.Hsapiens.v79")
ipak(required.packages)
required.bioconductor.packages = c("limma", "ComplexHeatmap","EnhancedVolcano")
ibiopak(required.bioconductor.packages)

Cancer = "BRCA"
#Quantile normalize lnc-RNASeq for TCGA
clustering <- read.table("./data/TCGA_Clustering.csv",header=TRUE,sep="\t")
expression.matrix <- read.table("./data/lncRNA_TCGA_patients.csv",header=TRUE,sep="\t")
expression.matrix <- as.matrix(expression.matrix)
expression.matrix[expression.matrix>quantile(expression.matrix,0.998)] <- as.numeric(quantile(expression.matrix,0.998))
write.table(expression.matrix,"./data/lncRNA_TCGA_patients_QN.csv",row.names=T,col.names=T,quote=F,sep=",")
save(expression.matrix,file = "./Validation.dataset/Proceesed_data/BRCA.quantile.normalized.matrix.Rdata")

#expression.matrix = expression.matrix[which(rowSums(expression.matrix) > 0),]

Expression.matrix = as.data.frame(t(expression.matrix))

rownames(Expression.matrix)= gsub("\\.","-",rownames(Expression.matrix))
     
Expression.matrix$ICR.Clusters = clustering$HL.ICR.Cluster[match(rownames(Expression.matrix), clustering$Sample_ID)]

table(Expression.matrix$ICR.Clusters)
Expression.matrix = Expression.matrix[order(Expression.matrix$ICR.Clusters,decreasing = TRUE),]

Expression.matrix$ICR.Clusters = factor(Expression.matrix$ICR.Clusters,levels = c("ICR-Low", "ICR-High"))

design=model.matrix(~ Expression.matrix$ICR.Clusters)
Expression.matrix$ICR.Clusters = NULL

Expression.matrix = as.matrix(Expression.matrix)
mode(Expression.matrix) = "numeric"

fit = lmFit(t(Expression.matrix) , design=design)
fit = eBayes(fit)
diff.stats = topTable(fit, coef=2,adjust.method = "BH", number = ncol(Expression.matrix))

#Get the mapping
Gtf_df = read.table("./data/gencode.v29.long_noncoding_RNAs.gtf", header = FALSE, sep = "\t")  # gtf downloaded from: https://www.gencodegenes.org/human/release_29.html
Gtf_df = Gtf_df[which(Gtf_df$V3 == "gene"),]
gene_info_list = strsplit(as.character(Gtf_df$V9),";",fixed = TRUE)
Gtf_sub = NULL
for (i in 1:length(gene_info_list))
{
  temp <- cbind(gene_info_list[[i]][1],gene_info_list[[i]][3])
  Gtf_sub <- rbind(Gtf_sub,temp)
}
Gtf_sub <- as.data.frame(Gtf_sub)
colnames(Gtf_sub) = c("gene_id", "gene_name")
Gtf_sub$gene_id = gsub("gene_id ", "", Gtf_sub$gene_id)
Gtf_sub$gene_name = gsub("gene_name ", "", Gtf_sub$gene_name)
Gtf_sub$gene_name = gsub(" ", "", Gtf_sub$gene_name)

# strip version number of gene_id
Gtf_sub$gene_id = str_replace(Gtf_sub$gene_id,
                              pattern = ".[0-9]+$",
                              replacement = "")

diff.stats$X <- str_replace(rownames(diff.stats), pattern = ".[0-9]+$",replacement="")
diff.stats$ENSEMBL_ID <- diff.stats$X
mapping_df <- select(EnsDb.Hsapiens.v79, key=diff.stats$X, 
                     columns=c("ENTREZID", "SYMBOL"), 
                     keytype="GENEID")

colnames(mapping_df) <- c("Ensemble_ID","ENTREZID","X")
gene_names <- NULL
for (i in 1:nrow(diff.stats))
{
  gene_names <- c(gene_names,mapping_df[mapping_df$Ensemble_ID==diff.stats$X[i],]$X[1])
}
diff.stats$gene_name <- gene_names

diff.genes = diff.stats[diff.stats$adj.P.Val < 0.05,]

dir.create(paste0("./Analysis/DEG"), showWarnings = FALSE)

write.csv(diff.stats, file = paste0("./Analysis/DEG/diff_lncRNAs_TCGA.csv"))
