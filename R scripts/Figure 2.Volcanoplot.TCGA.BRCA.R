
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd("~/")
source(paste0("~/R scripts/ipak.function.R"))
required.bioconductor.packages = c("EnhancedVolcano")
ibiopak(required.bioconductor.packages)

## Create the volcano plot 

#load data
DIFF = read.csv("./Analysis/DEG/diff_lncRNAs_TCGA.csv")

DIFF$gene_name = str_trim(DIFF$gene_name)

dev.new()

svg(filename = paste0("./Figures/Volcano.DELncRNAs.limma.TCGA.labelled.svg"), width = 7, height = 7,pointsize = 12)

# Creating the Volcano plot
EnhancedVolcano(DIFF,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,FCcutoff = 1,
                selectLab = DIFF$gene_name[DIFF$adj.P.Val<= 1e-12 & abs(DIFF$logFC) >= 1],
                labSize = 6.0,
                pointSize = 2.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
)

dev.off()
