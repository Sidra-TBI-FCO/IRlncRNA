####################################################################
###
### LncRNAs2Pathways script to be performed on HPC
### 
####################################################################

rm(list=ls())
setwd("/export/cse/rmall/Network_Analysis/Julie_Related/")                                                                    # Setwd to location were output files have to be saved.

#trace(utils:::unpackPkgZip, edit= TRUE)

#required.packages = c("e1071", "doSNOW", "ipred", "xgboost", "caret", "foreach", "ggplot2")
#install.packages(required.packages)
#install.packages("ddalpha", repos = "http://172.32.75.19")
library(ggplot2)
library(plyr)
library(LncPath)
library(stringr)

# Set parameters
version = "GRCh38"

## Load Data
diff.genes <- read.csv("../results/diff_lncRNAs_Arab.csv",stringsAsFactors = FALSE,sep=" ")
#diff.genes = read.csv("../results/diff_lncRNAs_TCGA.csv",stringsAsFactors = FALSE,sep=" ")

if(version == "GRCh38"){
  Gtf_df = read.table("../data/gencode.v29.long_noncoding_RNAs.gtf", header = FALSE, sep = "\t")  # gtf downloaded from: https://www.gencodegenes.org/human/release_29.html
  Gtf_df = Gtf_df[which(Gtf_df$V3 == "gene"),]
  gene_info_list = strsplit(as.character(Gtf_df$V9),";",fixed = TRUE)
  Gtf_sub = NULL
  for (i in 1:length(gene_info_list))
  {
    temp <- cbind(gene_info_list[[i]][1],gene_info_list[[i]][3])
    Gtf_sub <- rbind(Gtf_sub,temp)
  }
  Gtf_sub <- as.data.frame(Gtf_sub)
  #Gtf_sub = data.frame(do.call('rbind', strsplit(as.character(Gtf_df$V9),';',fixed=TRUE)))[c(1,5)]
}
if(version == "GRCh37"){
  Gtf_df = read.table("../Tools/gencode.v19.long_noncoding_RNAs.gtf", header = FALSE, sep = "\t")  # gtf downloaded from: https://www.gencodegenes.org/human/release_29.html
  Gtf_df = Gtf_df[which(Gtf_df$V3 == "gene"),]
  Gtf_sub = data.frame(do.call('rbind', strsplit(as.character(Gtf_df$V9),';',fixed=TRUE)))[c(1,4)]
}

colnames(Gtf_sub) = c("gene_id", "gene_name")
Gtf_sub$gene_id = gsub("gene_id ", "", Gtf_sub$gene_id)
Gtf_sub$gene_name = gsub("gene_name ", "", Gtf_sub$gene_name)
Gtf_sub$gene_name = gsub(" ", "", Gtf_sub$gene_name)

# strip version number of gene_id
Gtf_sub$gene_id = str_replace(Gtf_sub$gene_id,
                              pattern = ".[0-9]+$",
                              replacement = "")

# convert to ensemble IDs for Arabs
diff.genes$gene_name = Gtf_sub$gene_name[match(diff.genes$Gene, Gtf_sub$gene_name)]
diff.genes.annotated = diff.genes[which(diff.genes$Gene %in% Gtf_sub$gene_name),]   # not all of the lncRNA's are annotated
diff.genes.annotated$Ensemble_ID = Gtf_sub$gene_id[match(diff.genes.annotated$gene_name, Gtf_sub$gene_name)]
write.table(diff.genes.annotated,"../data/lncRNA_Id_to_Ensemble_Id_mapping.csv",col.names=T,row.names=F,quote=F,sep=",")

#Use it for TCGA
#diff.genes$Gene <- str_replace(diff.genes$Gene, pattern = ".[0-9]+$",replacement="")
#diff.genes$gene_name <- diff.genes$Gene
#diff.genes.annotated <- diff.genes[diff.genes$adj.P.Val<=0.05,]
#diff.genes.annotated$Ensemble_ID <- diff.genes.annotated$gene_name

# Perform analysis
NetLncPath <- getNet() #Get the background lncRNA-mRNA interaction network, 
#it was built by intergrating an lncRNAmRNA co-expression network and the protein-protein interaction network.

length(unique(NetLncPath$V1)) # 21107 unique values in column 1
length(unique(NetLncPath$V2)) # 27375 unique values in column 2

#get example lncRNA sets 
SigLncs <- diff.genes.annotated$Ensemble_ID
print(head(SigLncs), row.names = FALSE)

#Perform the Walkscore to get score for each gene (pcg/mRNA)
if (length(SigLncs) == 0) 
  stop("The list is empty.")
cat("Now start the random walking...\n")
if (!exists("LncPathEnvir")) 
  LncPathEnvir <- initializeLncPathEnvir()
Network <- NetLncPath
NetLncPath <- graph.edgelist(as.matrix(Network), directed = FALSE)
VertexWeight <- rep(0, length(V(NetLncPath)))
names(VertexWeight) <- V(NetLncPath)$name

# Can even consider a weighted start set
# for (i in 1:length(SigLncs))
# {
#   geneid <- SigLncs[i]
#   pval_info <- -log10(diff.genes.annotated[diff.genes.annotated$Ensemble_ID == geneid,]$"P.Value")
#   VertexWeight[V(NetLncPath)$name == geneid] <- pval_info
# }

VertexWeight[V(NetLncPath)$name %in% SigLncs] <- 1
print(names(which(VertexWeight > 0)))
cat(paste(length(which(VertexWeight > 0)), "of", length(SigLncs), 
          "were mapped to the huge net with", length(VertexWeight), 
          "Nodes.", "\n"))

WalkRes <- RandomWalk2igraph(NetLncPath, VertexWeight, EdgeWeight = FALSE)
WalkScore <- as.data.frame(cbind(V(NetLncPath)$name, WalkRes))
colnames(WalkScore) <- c("Gene","WalkRes")
WalkScore$Gene <- as.character(as.vector(WalkScore$Gene))
WalkScore$WalkRes <- as.numeric(as.vector(WalkScore$WalkRes))

if (length(which(as.numeric(WalkScore[[2]]) == 0)) > 0) 
  WalkScore <- WalkScore[-(which(as.numeric(WalkScore[[2]]) == 
                                   0)), ]
WalkScore[[2]] <- sqrt(as.numeric(WalkScore[[2]]))
WalkScore <- WalkScore[order(WalkScore[[2]], decreasing = TRUE), 
                       ]
PCEnsem2Sym <- get("PCEnsem2Sym", envir = LncPathEnvir)
PCWalkScore <- merge(PCEnsem2Sym, WalkScore, by.x = 2, by.y = 1)
PCWalkScore <- PCWalkScore[order(PCWalkScore[[3]], decreasing = TRUE), 
                           ]

differntial_PCWalkScore <- PCWalkScore[PCWalkScore$WalkRes>=0.01,]
write.csv(differntial_PCWalkScore,"../results/Paper_Text_Results/prop_genes_WalkScore_Arab.csv", row.names = F, quote=F)
write.csv(PCWalkScore,"../results/Paper_Text_Results/all_genes_WalkScore_Arab.csv",row.names=F,quote=F)

# NetLncPath <- getNet() #Get the background lncRNA-mRNA interaction network 
# Result = lncPath(SigLncs, NetLncPath, Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 1000, minPathSize = 15, maxPathSize = 500)
# Table = lncPath2Table(Result)
# save(list = c("Table"), file = "results/TCGA.Pathways.KEGG.completed.Rdata")
# 
# #Load Walkscore for RAQA
# df_walkscore <- read.table("results/20200604_all_genes_WalkScore_TCGA.csv",header=TRUE,sep=",")
