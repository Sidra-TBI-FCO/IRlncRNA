# final script for the Normalization of the RNA Seq data together , (NBL =161)

rm(list=ls())

load("~/R.Config.Rdata")
setwd(paste0(master.location.db,"/"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

load("./GeneInfo/geneInfo.July2017.RData")

geneInfo["CCL5","ensembl_gene_id"][which(geneInfo["CCL5","ensembl_gene_id"] %in% "ENSG00000161570")] = "ENSG00000271503"


required.bioconductor.packages = c("SummarizedExperiment","EDASeq", "preprocessCore","limma","TCGAbiolinks")
required.packages = c("base64enc", "HGNChelper","RCurl","httr","stringr","digest","bitops",
                      "rjson","Matrix","latticeExtra","matrixStats")
ipak(required.packages)
ibiopak(required.bioconductor.packages)

 # Parameters
Cancer= "TCGA-BRCA"  # "OV" , "COAD" , "READ"

load(paste0("./Processed.data/",Cancer,"/001_Exp_",Cancer,"_Processed_Data.Rdata"))

data = data[-c(60484:60488),]
data = as.data.frame(data)
data$X1 = substr(data$X1,1,15)
data$X2 = rownames(geneInfo)[match(data$X1, geneInfo$ensembl_gene_id)]
data = data[!is.na(data$X2),]
data = data[!duplicated(data$X2),]
rownames(data) = data$X2
data$X1 = NULL
data$X2 = NULL

RNA_To_Normalize_all = data
cat(paste0("Number of samples in Data is ", ncol(RNA_To_Normalize_all), ".\n", "number of genes in Data is ", 
           nrow(RNA_To_Normalize_all), ".", "\n", "Number of patients is ", 
           length(unique(substring(colnames(RNA_To_Normalize_all),1,16))), "."))

RNAseqData.to.normalize = RNA_To_Normalize_all
RNAseqData.to.normalize = as.matrix(RNAseqData.to.normalize)

# drop rows(genes) without approved gene-symbol
RNAseq.genes = rownames(RNAseqData.to.normalize)
info.genes = rownames(geneInfo)
available.genes = unique(RNAseq.genes[which(RNAseq.genes %in% info.genes)])
geneInfo = geneInfo[which(rownames(geneInfo) %in% available.genes),]                                                      # drop the genes without RNAseq.DATA
RNAseqData.to.normalize = RNAseqData.to.normalize[which(rownames(RNAseqData.to.normalize) %in% available.genes),]         # drop the genes without info
mode(RNAseqData.to.normalize) <- "numeric"

geneInfo <- geneInfo[ order(row.names(geneInfo)), ]
RNAseqData.to.normalize <- floor(RNAseqData.to.normalize[order(row.names(RNAseqData.to.normalize)),])

cat(paste0(length(RNAseqData.to.normalize[,1]), " genes were available for normalization."))
# get information on number of genes included for normalization to logfile

# remove the duplicated column names
colnames(RNAseqData.to.normalize) = substring(colnames(RNAseqData.to.normalize),1,12)
RNA_To_Normalize_all=as.matrix(RNA_To_Normalize_all)
RNAseqData.to.normalize <- RNAseqData.to.normalize[, !duplicated(colnames(RNAseqData.to.normalize))]

RNASeq.expr.set = newSeqExpressionSet(RNAseqData.to.normalize, featureData = geneInfo)                              # Create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "gcContent"] = as.numeric(geneInfo[, "gcContent"])                                               # Make sure gcContenet is numeric
RNASeq.expr.set = withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE)             # Removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set = betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)                         # Removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM = log(RNAseqData.to.normalize + .1) + offst(RNASeq.expr.set)                                      # Apply the Edaseq Ofset
RNASeq.NORM = floor(exp(RNASeq.NORM) - .1)                                                                                # Return non decimal values

#Quantile normalization RNA
RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                                                                 # Quantile normalize
RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                                                                     # Return non decimal values
rownames(RNASeq.NORM.quantiles) <- rownames(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)

dev.new()
plotDensities(log(RNASeq.NORM.quantiles[,1:1091],2))

save(RNASeq.NORM.quantiles,geneInfo,file=paste0("./Processed.data/",Cancer,"/003_",Cancer,"_EDAseq_normalized.gene.info.2017.Rdata"))

