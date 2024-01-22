# 1. Identification of differentially activated TFs in each CCS
# maligant.seurat: seurat object of malignant cells

library(Matrix)
library(Seurat)
library(dplyr)
library(foreach)
library(SCENIC) 
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCopeLoomR)

# Initialize the SCENIC object
exprMat <- as.matrix(maligant.seurat@assays$RNA@counts)
hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org="hgnc", dbDir="scenic_db", dbs = hg38_dbs, nCores=1)

# Filter out genes
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

# Infer Co-expression
runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)

# Infer TF-target relationship
runGenie3(exprMat_filtered_log, scenicOptions)

# Run SCENIC
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat_filtered_log)

# Get the nonextended (high confidence) regulons
aucell_matrix <-getAUC(loadInt(scenicOptions,"aucell_regulonAUC"))      
aucell_matrix <- aucell_matrix[ -grep("extended", rownames(aucell_matrix)), ]
 


# 2. Differential regulon analysis
# clusters: a list of cells in CCSs
# aucell_matrix: a matrix, activity of high confidence regulons in cells

# Calculate log2FC, p.val and p.adj
library(stats)
log2FC <- c()
p.val <- c()
for(i in 1: length(clusters)){
  a <- aucell_matrix[,clusters[[i]]]
  b <- aucell_matrix[,-which(colnames(aucell_matrix) %in% clusters[[i]])]
  
  for(j in 1:dim(aucell_matrix)[1]){
    fc <- mean(a[j,])/mean(b[j,])
    log2FC <- c(log2FC,log2(fc))
    
    wil <- wilcox.test(a[j,],b[j,])
    
    p.val <- c(p.val,wil$p.value)
  }
}
p.adj <- p.adjust(p.val,method = "fdr")

diff.regulon <- data.frame(CCS = rep(names(clusters), each = dim(aucell_matrix)[1]), 
                           regulon = rep(rownames(aucell_matrix), times = length(clusters)),
                           log2FC = log2FC, p.val = p.val, p.adj = p.adj)
diff.regulon <- split(diff.regulon, diff.regulon$CCS)

# log2FC > 1 and p.adj < 0.05
# Filter out top 5 regulons in each CCS

library(tibble)
library(dplyr)

temp <- list()
top5 <- list()
for(i in length(clusters)){
  temp[[i]] <- subset(diff.regulon[[i]],log2FC > 1 & p.adj < 0.05)
  top5[[i]] <- arrange(temp[[i]], desc(log2FC))[1:5,]
}














