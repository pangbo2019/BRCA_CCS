# 1. Association analysis: spearman
#' @param maligant.seurat seurat object of malignant cells
#' @param deg.list the list of differentially expressed genes

library(Seurat)
library(stats)
# calculate the average expression of differential genes  
gene <- AverageExpression(maligant.seurat, group.by = "harmony_clusters", assays = "RNA")
gene <- gene[[1]]
avg <- list()
for(i in 1: 11){
  avg [[i]] <- apply(gene[deg.list[[i]],], 2, mean)
}
avg.mat <- matrix(unlist(avg), nrow = 11)
rownames(avg.mat) <- names(avg[[1]])

# spearman
corr <- cor(t(avg.mat),use = "all.obs",method = "spearman")



